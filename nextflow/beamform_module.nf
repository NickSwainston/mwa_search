nextflow.preview.dsl = 2


params.obsid = null
params.pointings = null
params.calid = null

params.begin = 0
params.end = 0
params.all = false

params.summed = false
params.vcstools_version = 'master'
params.mwa_search_version = 'master'

params.basedir = '/group/mwaops/vcs'
params.scratch_basedir = '/astro/mwaops/vcs'
params.didir = "${params.basedir}/${params.obsid}/cal/${params.calid}/rts"
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
mb_ipfb_dur = ( obs_length * (params.bm_read + 3 * (params.bm_cal + params.bm_beam) + params.bm_write) + 20 ) * 2

//Beamforming duration calc
mb_dur = ( obs_length * (params.bm_read + params.bm_cal + max_job_pointings * (params.bm_beam +params.bm_write)) + 20 ) * 2

//Required temp SSD mem required for gpu jobs
temp_mem = (int) (0.0012 * obs_length * max_job_pointings + 1)
if ( ! params.summed ) {
    temp_mem = temp_mem * 4
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



// Handling begin and end times
process get_beg_end {
    script:
    if ( params.all )
        """
        #!/usr/bin/env python3

        from mwa_metadb_utils import obs_max_min

        beg, end = obs_max_min(${params.obsid})
        print("{},{}".format(beg, end), end="")
        """
    else
        """
        #!/usr/bin/env python3

        beg = "$params.begin"
        end = "$params.end"
        print("{},{}".format(beg, end), end="")
        """
}


process get_channels {
    output:
    file "${params.obsid}_channels.txt"

    """
    #!/usr/bin/env python3

    from mwa_metadb_utils import get_channels
    import csv

    channels = get_channels($params.obsid)
    with open("${params.obsid}_channels.txt", "w") as outfile:
        spamwriter = csv.writer(outfile, delimiter=',')
        spamwriter.writerow(channels)
    """
}


process ensure_metafits {

    """
    #!/usr/bin/env python3

    from process_vcs import ensure_metafits
    
    ensure_metafits("${params.basedir}/${params.obsid}", "${params.obsid}",\
                    "${params.basedir}/${params.obsid}/${params.obsid}_metafits_ppds.fits")
    """
}


process gps_to_utc {
    input:
    tuple val(begin), val(end)

    """
    #!/usr/bin/env python3

    from process_vcs import gps_to_utc

    print(gps_to_utc(${begin})),
    """
}


process make_directories {
    """
    #!/usr/bin/env python3

    from mdir import mdir
    from process_vcs import create_link

    mdir("${params.basedir}/${params.obsid}", "Data")
    mdir("${params.scratch_basedir}/${params.obsid}", "Products")
    mdir("${params.basedir}/batch", "Batch")
    mdir("${params.basedir}/${params.obsid}/pointings", "Pointings")
    mdir("${params.scratch_basedir}/${params.obsid}/dpp_pointings", "DPP Products")
    create_link("${params.basedir}/${params.obsid}", "dpp_pointings",
                "${params.scratch_basedir}/${params.obsid}", "dpp_pointings")
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
    maxRetries 3
    if ( "$HOSTNAME".startsWith("farnarkle") ) {
        clusterOptions = "--gres=gpu:1  --tmp=${temp_mem}GB"
        scratch '$JOBFS'
    }

    input:
    each channel_pair
    val utc
    each point
    tuple val(begin), val(end)

    output:
    file "*/*fits"

    beforeScript "module use $params.module_dir; module load vcstools/$params.vcstools_version"

    //TODO add other beamform options and flags -F
    """
    make_beam -o $params.obsid -b $begin -e $end -a 128 -n 128 \
-f ${channel_pair[0]} -J ${params.didir}/DI_JonesMatrices_node${channel_pair[1]}.dat \
-d ${params.basedir}/${params.obsid}/combined -P ${point.join(",")} \
-r 10000 -m ${params.basedir}/${params.obsid}/${params.obsid}_metafits_ppds.fits \
${bf_out} -z $utc
    """
}


process make_beam_ipfb {
    publishDir "${params.basedir}/${params.obsid}/pointings/${point}", mode: 'move', enabled: params.publish_fits, pattern: "*hdr"
    publishDir "${params.basedir}/${params.obsid}/pointings/${point}", mode: 'move', enabled: params.publish_fits, pattern: "*vdif"
    publishDir "${params.scratch_basedir}/${params.obsid}/dpp_pointings/${point}", mode: 'move', enabled: params.publish_fits_scratch, pattern: "*hdr"
    publishDir "${params.scratch_basedir}/${params.obsid}/dpp_pointings/${point}", mode: 'move', enabled: params.publish_fits_scratch, pattern: "*vdif"

    label 'gpu'
    //time '2h'
    time "${mb_ipfb_dur}s"
    errorStrategy 'retry'
    maxRetries 3
    if ( "$HOSTNAME".startsWith("farnarkle") ) {
        clusterOptions = "--gres=gpu:1  --tmp=${temp_mem}GB"
        scratch '$JOBFS'
    }

    input:
    each channel_pair
    val utc
    each point
    tuple val(begin), val(end)

    output:
    file "*fits"
    //file "*hdr"
    //file "*vdif"

    beforeScript "module use $params.module_dir; module load vcstools/origbeam"

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
-d ${params.basedir}/${params.obsid}/combined -R ${point.split("_")[0]} -D ${point.split("_")[1]} \
-r 10000 -m ${params.basedir}/${params.obsid}/${params.obsid}_metafits_ppds.fits \
-p -u -z $utc
    """
}

process splice {
    publishDir "${params.basedir}/${params.obsid}/pointings/${unspliced[0].baseName.split("_")[2]}_${unspliced[0].baseName.split("_")[3]}", mode: 'copy', enabled: params.publish_fits
    publishDir "${params.scratch_basedir}/${params.obsid}/dpp_pointings/${unspliced[0].baseName.split("_")[2]}_${unspliced[0].baseName.split("_")[3]}", mode: 'copy', enabled: params.publish_fits_scratch
    label 'cpu'
    time '1h'

    input:
    val chan
    each file(unspliced)

    output:
    file "${params.obsid}*fits"
    val "${unspliced[0].baseName.split("_")[2]}_${unspliced[0].baseName.split("_")[3]}"

    beforeScript "module use $params.module_dir; module load vcstools/$params.vcstools_version; module load mwa_search/$params.mwa_search_version"

    """
    splice_wrapper.py -o ${params.obsid} -c ${chan.join(" ")}
    """
}


workflow pre_beamform {
    main:
        get_beg_end()
        get_channels()
        ensure_metafits()
        gps_to_utc( get_beg_end.out.map{ it.split(",") }.flatten().collect() )
        make_directories()
        combined_data_check(get_beg_end.out.map{ it.split(",") }.flatten().collect())
    emit:
        get_beg_end.out.map{ it.split(",") }.flatten().collect()
        get_channels.out.splitCsv()
        gps_to_utc.out
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
                make_beam.out.flatten().map { it -> [it.baseName.split("ch")[0], it ] }.groupTuple().map { it -> it[1] } )
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
                make_beam_ipfb.out[0].flatten().map { it -> [it.baseName.split("ch")[0], it ] }.groupTuple().map { it -> it[1] } )
    emit:
        make_beam_ipfb.out.flatten().map{ it -> [it.baseName.split("ch")[0], it ] }.groupTuple().map{ it -> it[1] }
        splice.out[0].flatten().map{ it -> [it.baseName.split("ch")[0], it ] }.groupTuple().map{ it -> it[1] }
        splice.out[1]
}
