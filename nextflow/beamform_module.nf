nextflow.preview.dsl = 2


params.obsid = null
params.pointings = null
params.calid = null

params.begin = 0
params.end = 0
params.all = false

params.summed = false
params.channels = null
params.vcstools_version = 'master'
params.mwa_search_version = 'master'

params.basedir = '/group/mwavcs/vcs'
params.scratch_basedir = '/astro/mwavcs/vcs'
params.didir = "${params.scratch_basedir}/${params.obsid}/cal/${params.calid}/rts"
params.publish_fits = false
params.publish_fits_scratch = false

params.no_combined_check = false

//Calculate the max pointings used in the launched jobs
if ( params.pointings ) {
    max_job_pointings = params.pointings.split(",").size()
    if ( max_job_pointings > params.max_pointings ) {
        max_job_pointings = params.max_pointings
    }
}
else {
    // No input pointings (this happens in beamform_fov_sources.nf) so assuming max
    max_job_pointings = params.max_pointings
}

//Work out total obs time
if ( params.all ) {
    // an estimation since there's no easy way to make this work
    obs_length = 4805
}
else {
    obs_length = params.end - params.begin + 1
}

//Calculate expected number of fits files
n_fits = (int) (obs_length/200)
if ( obs_length % 200 != 0 ) {
    n_fits = n_fits + 1
}


//Beamforming ipfb duration calc
mb_ipfb_dur = ( obs_length * (params.bm_read + 3 * (params.bm_cal + params.bm_beam) + params.bm_write) + 200 ) * 1.2

//Beamforming duration calc
mb_dur = ( obs_length * (params.bm_read + params.bm_cal + max_job_pointings * (params.bm_beam +params.bm_write)) + 200 ) * 1.2

//Required temp SSD mem required for gpu jobs
temp_mem = (int) (0.0012 * obs_length * max_job_pointings + 1)
temp_mem_single = (int) (0.0024 * obs_length + 2)
if ( ! params.summed ) {
    temp_mem = temp_mem * 4
    temp_mem_single = temp_mem_single *4
}


if ( params.summed ) {
    bf_out = " -p -s "
}
else {
    bf_out = " -p "
}

range = Channel.from( ['001', '002', '003', '004', '005', '006',\
                       '007', '008', '009', '010', '011', '012',\
                       '013', '014', '015', '016', '017', '018',\
                       '019', '020', '021', '022', '023', '024'] )


process beamform_setup {
    output:
    file "${params.obsid}_beg_end.txt"
    file "${params.obsid}_channels.txt"
    file "${params.obsid}_utc.txt"

    """
    #!/usr/bin/env python
    import csv

    from mwa_metadb_utils import obs_max_min, get_channels
    from process_vcs import ensure_metafits, gps_to_utc, create_link
    from mdir import mdir
    import csv

    # Work out begin and end time of obs
    if "${params.all}" == "true":
        beg, end = obs_max_min(${params.obsid})
    else:
        beg = $params.begin
        end = $params.end
    with open("${params.obsid}_beg_end.txt", "w") as outfile:
        spamwriter = csv.writer(outfile, delimiter=',')
        spamwriter.writerow([beg, end])

    # Find the channels 
    if "$params.channels" is "null":
        channels = get_channels($params.obsid)
    else:
        channels = [$params.channels]
    with open("${params.obsid}_channels.txt", "w") as outfile:
        spamwriter = csv.writer(outfile, delimiter=',')
        spamwriter.writerow(channels)
        
    # Ensure the metafits files is there
    ensure_metafits("${params.basedir}/${params.obsid}", "${params.obsid}",\
                    "${params.scratch_basedir}/${params.obsid}/${params.obsid}_metafits_ppds.fits")

    # Covert gps time to utc
    with open("${params.obsid}_utc.txt", "w") as outfile:
        spamwriter = csv.writer(outfile, delimiter=',')
        spamwriter.writerow([gps_to_utc(beg)])

    # Make sure all the required directories are made
    mdir("${params.scratch_basedir}/${params.obsid}", "Data")
    mdir("${params.scratch_basedir}/${params.obsid}", "Products")
    mdir("${params.scratch_basedir}/batch", "Batch")
    mdir("${params.scratch_basedir}/${params.obsid}/pointings", "Pointings")
    """
}

process combined_data_check {
    when:
    params.no_combined_check == false

    input:
    tuple val(begin), val(end)

    """
    #!/usr/bin/env python

    import sys
    from check_known_pulsars import check_data

    #Perform data checks
    dur = $end-$begin + 1
    check = check_data("$params.obsid", beg=$begin, dur=dur)
    if not check:
        print("ERROR: Recombined check has failed. Cannot continue.")
        sys.exit(1)
    else:
        print("Recombined check passed, all files present.")
    """
}


process make_beam {
    label 'gpu'
    //time '2h'
    time "${mb_dur}s"
    errorStrategy 'retry'
    maxRetries 1
    if ( "$HOSTNAME".startsWith("garrawarla") ) {
        maxForks 70
    }
    else {
        maxForks 120
    }

    if ( "$HOSTNAME".startsWith("farnarkle") ) {
        clusterOptions = "--gres=gpu:1  --tmp=${temp_mem}GB"
        scratch '$JOBFS'
        beforeScript "module use ${params.module_dir}; module load vcstools/${params.vcstools_version}"
    }
    else if ( "$HOSTNAME".startsWith("x86") ) {
        clusterOptions = "--gres=gpu:1"
        scratch '/ssd'
        //container = "file:///${params.containerDir}/vcstools/vcstools_${params.vcstools_version}.sif"
        beforeScript "module use ${params.module_dir}; module load vcstools/${params.vcstools_version}"
    }
    else if ( "$HOSTNAME".startsWith("garrawarla") ) {
        clusterOptions = "--gres=gpu:1  --tmp=${temp_mem}GB"
        scratch '/nvmetmp'
        container = "file:///${params.containerDir}/vcstools/vcstools_${params.vcstools_version}.sif"
    }
    else if ( "$HOSTNAME".startsWith("galaxy") ) {
        beforeScript "module use ${params.module_dir}; module load vcstools/${params.vcstools_version}"
    }
    else {
        container = "cirapulsarsandtransients/vcstools:${params.vcstools_version}"
    }

    input:
    each channel_pair
    val utc
    each point
    tuple val(begin), val(end)

    output:
    file "*fits"


    //TODO add other beamform options and flags -F
    """
    make_beam -o $params.obsid -b $begin -e $end -a 128 -n 128 \
-f ${channel_pair[0]} -J ${params.didir}/DI_JonesMatrices_node${channel_pair[1]}.dat \
-d ${params.scratch_basedir}/${params.obsid}/combined -P ${point.join(",")} \
-r 10000 -m ${params.scratch_basedir}/${params.obsid}/${params.obsid}_metafits_ppds.fits \
${bf_out} -t 6000 -z $utc
    mv */*fits .
    """
}


process make_beam_ipfb {
    publishDir "${params.basedir}/${params.obsid}/pointings/${point}", mode: 'copy', enabled: params.publish_fits, pattern: "*hdr"
    publishDir "${params.basedir}/${params.obsid}/pointings/${point}", mode: 'copy', enabled: params.publish_fits, pattern: "*vdif"
    publishDir "${params.scratch_basedir}/${params.obsid}/pointings/${point}", mode: 'copy', enabled: params.publish_fits_scratch, pattern: "*hdr"
    publishDir "${params.scratch_basedir}/${params.obsid}/pointings/${point}", mode: 'copy', enabled: params.publish_fits_scratch, pattern: "*vdif"

    label 'gpu'
    //time '2h'
    time "${mb_ipfb_dur}s"
    errorStrategy 'retry'
    maxRetries 1
    if ( "$HOSTNAME".startsWith("garrawarla") ) {
        maxForks 70
    }
    else {
        maxForks 120
    }
    
    if ( "$HOSTNAME".startsWith("farnarkle") ) {
        clusterOptions = "--gres=gpu:1  --tmp=${temp_mem_single}GB"
        scratch '$JOBFS'
        beforeScript "module use ${params.module_dir}; module load vcstools/${params.vcstools_version}"
    }
    else if ( "$HOSTNAME".startsWith("x86") ) {
        clusterOptions = "--gres=gpu:1"
        scratch '/ssd'
        container = "file:///${params.containerDir}/vcstools/vcstools_${params.vcstools_version}.sif"
    }
    else if ( "$HOSTNAME".startsWith("garrawarla") ) {
    clusterOptions = "--gres=gpu:1  --tmp=${temp_mem_single}GB"
        scratch '/nvmetmp'
        container = "file:///${params.containerDir}/vcstools/vcstools_${params.vcstools_version}.sif"
    }
    else if ( "$HOSTNAME".startsWith("galaxy") ) {
        beforeScript "module use ${params.module_dir}; module load vcstools/${params.vcstools_version}"
    }
    else {
        container = "cirapulsarsandtransients/vcstools:${params.vcstools_version}"
    }

    when:
    point != " " //Don't run if blank pointing given

    input:
    each channel_pair
    val utc
    each point
    tuple val(begin), val(end)

    output:
    file "*fits"
    //file "*hdr"
    //file "*vdif"

    //TODO add other beamform options and flags -F
    """
    if $params.publish_fits; then
        mkdir -p -m 771 ${params.basedir}/${params.obsid}/pointings/${point}
    fi
    if $params.publish_fits_scratch; then
        mkdir -p -m 771 ${params.scratch_basedir}/${params.obsid}/dpp_pointings/${point}
    fi

    make_beam -o $params.obsid -b $begin -e $end -a 128 -n 128 \
-f ${channel_pair[0]} -J ${params.didir}/DI_JonesMatrices_node${channel_pair[1]}.dat \
-d ${params.scratch_basedir}/${params.obsid}/combined -P ${point} \
-r 10000 -m ${params.scratch_basedir}/${params.obsid}/${params.obsid}_metafits_ppds.fits \
-p -v -t 6000 -z $utc
    ls *
    mv */*fits .
    """
}

process splice {
    publishDir "${params.basedir}/${params.obsid}/pointings/${unspliced[0].baseName.split("_")[2]}_${unspliced[0].baseName.split("_")[3]}", mode: 'copy', enabled: params.publish_fits
    publishDir "${params.scratch_basedir}/${params.obsid}/pointings/${unspliced[0].baseName.split("_")[2]}_${unspliced[0].baseName.split("_")[3]}", mode: 'copy', enabled: params.publish_fits_scratch
    label 'cpu'
    time '2h'
    maxForks 300
    errorStrategy 'retry'
    maxRetries 1

    input:
    val chan
    each file(unspliced)

    output:
    file "${params.obsid}*fits"
    val "${unspliced[0].baseName.split("_")[2]}_${unspliced[0].baseName.split("_")[3]}"

    if ( "$HOSTNAME".startsWith("farnarkle") ) {
        beforeScript "module use ${params.module_dir}; module load vcstools/${params.vcstools_version}"
    }
    else if ( "$HOSTNAME".startsWith("x86") ) {
        beforeScript "module use ${params.module_dir}; module load vcstools/${params.vcstools_version}"
    }
    else if ( "$HOSTNAME".startsWith("garrawarla") ) {
        container = "file:///${params.containerDir}/vcstools/vcstools_${params.vcstools_version}.sif"
    }
    else if ( "$HOSTNAME".startsWith("galaxy") ) {
        beforeScript "module use ${params.module_dir}; module load vcstools/${params.vcstools_version}"
    }
    else {
        container = "cirapulsarsandtransients/vcstools:${params.vcstools_version}"
    }

    """
    splice_wrapper.py -o ${params.obsid} -c ${chan.join(" ")}
    """
}


workflow pre_beamform {
    main:
        beamform_setup()
        combined_data_check(beamform_setup.out[0].splitCsv())
    emit:
        beamform_setup.out[0].splitCsv()
        beamform_setup.out[1].splitCsv()
        beamform_setup.out[2].splitCsv().flatten()
}


workflow beamform {
    take: 
        obs_beg_end
        channels
        utc
        pointings
    main:
        make_beam( channels.flatten().merge(range),\
                   utc,\
                   pointings,\
                   obs_beg_end )
        splice( channels,\
                make_beam.out.flatten().map { it -> [it.baseName.split("ch")[0], it ] }.\
                groupTuple( size: 24 ).map { it -> it[1] } )
    emit:
        make_beam.out.flatten().map{ it -> [it.baseName.split("ch")[0], it ] }.groupTuple().map{ it -> it[1] }
        splice.out[0].flatten().map{ it -> [it.baseName.split("ch")[0], it ] }.groupTuple().map{ it -> it[1] }
        splice.out[1]
        splice.out[0] | flatten() | map { it -> [it.baseName.split("_ch")[0].split("${params.obsid}_")[-1], it ] } | groupTuple()
}

workflow beamform_ipfb {
    take:
        obs_beg_end
        channels
        utc
        pointings
    main:
        make_beam_ipfb( channels.flatten().merge(range),\
                        utc,\
                        pointings.flatten(),\
                        obs_beg_end )
        splice( channels,\
                make_beam_ipfb.out[0].flatten().map { it -> [it.baseName.split("ch")[0], it ] }.\
                groupTuple( size: 24 ).map { it -> it[1] } )
    emit:
        make_beam_ipfb.out.flatten().map{ it -> [it.baseName.split("ch")[0], it ] }.groupTuple().map{ it -> it[1] }
        splice.out[0].flatten().map{ it -> [it.baseName.split("ch")[0], it ] }.groupTuple().map{ it -> it[1] }
        splice.out[1]
}
