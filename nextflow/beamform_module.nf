nextflow.preview.dsl = 2


params.obsid = null
params.pointings = null
params.calid = null

params.begin = null
params.end = null
params.all = false

params.summed = false
params.vcstools_version = 'master'

params.basedir = '/group/mwaops/vcs'
params.stratch_basedir = '/astro/mwaops/vcs'
params.didir = "${params.basedir}/${params.obsid}/cal/${params.calid}/rts"
params.channels = null
params.publish_fits = false
params.publish_fits_scratch = false


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
    //when:
    //params.all == true

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


process make_beam {
    label 'gpu'
    time '10h'
    errorStrategy 'retry'
    maxRetries 3

    input:
    each channel_pair
    val utc
    each point
    tuple val(begin), val(end)

    output:
    file "*/*fits"

    //TODO add other beamform options and flags -F
    """
    module use /group/mwa/software/modulefiles
    module load vcstools/$params.vcstools_version
    make_beam -o $params.obsid -b $begin -e $end -a 128 -n 128 \
-f ${channel_pair[0]} -J ${params.didir}/DI_JonesMatrices_node${channel_pair[1]}.dat \
-d ${params.basedir}/${params.obsid}/combined -P ${point.join(",")} \
-r 10000 -m ${params.basedir}/${params.obsid}/${params.obsid}_metafits_ppds.fits \
${bf_out} -z $utc
    """
}


process make_beam_ipfb {
    label 'gpu'
    time '10h'
    errorStrategy 'retry'
    maxRetries 3

    input:
    each channel_pair
    val utc
    each point
    tuple val(begin), val(end)

    output:
    file "*fits"

    //TODO add other beamform options and flags -F
    """
    module use /group/mwa/software/modulefiles
    module load vcstools/origbeam
    make_beam -o $params.obsid -b $begin -e $end -a 128 -n 128 \
-f ${channel_pair[0]} -J ${params.didir}/DI_JonesMatrices_node${channel_pair[1]}.dat \
-d ${params.basedir}/${params.obsid}/combined -R ${point.split("_")[0]} -D ${point.split("_")[1]}\
-r 10000 -m ${params.basedir}/${params.obsid}/${params.obsid}_metafits_ppds.fits \
${bf_out} -u -z $utc
    """
}

process splice {
    publishDir "${params.basedir}/${params.obsid}/pointings/${unspliced[0].baseName.split("_")[2]}_${unspliced[0].baseName.split("_")[3]}", mode: 'move', enabled: params.publish_fits
    publishDir "${params.scratch_basedir}/${params.obsid}/dpp_pointings/${unspliced[0].baseName.split("_")[2]}_${unspliced[0].baseName.split("_")[3]}", mode: 'move', enabled: params.publish_fits_scratch
    label 'cpu'
    time '1h'

    input:
    val chan
    each file(unspliced)

    output:
    file "${params.obsid}*fits"
    val "${unspliced[0].baseName.split("_")[2]}_${unspliced[0].baseName.split("_")[3]}"

    """
    module use /group/mwa/software/modulefiles
    module load mwa_search
    module load vcstools
    splice_wrapper.py -o ${params.obsid} -c ${chan.join(" ")}
    """
}


workflow pre_beamform {
    main:
        get_beg_end()
        get_channels()
        ensure_metafits()
        gps_to_utc( get_beg_end.out.map{ it.split(",") }.flatten().collect() )
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
        make_beam( channels | flatten() | merge(range),\
                   utc,\
                   pointings,\
                   obs_beg_end )
        splice( channels,\
                make_beam.out | flatten() | map { it -> [it.baseName.split("ch")[0], it ] } | groupTuple() | map { it -> it[1] } )
    emit:
        splice.out[0]
        splice.out[1]
}

workflow beamform_ipfb {
    take: 
        obs_beg_end
        channels
        utc
        pointings
    main:
        make_beam_ipfb( channels | flatten() | merge(range),\
                        utc,\
                        pointings.flatten(),\
                        obs_beg_end )
        splice( channels,\
                make_beam_ipfb.out | flatten() | map { it -> [it.baseName.split("ch")[0], it ] } | groupTuple() | map { it -> it[1] } )
    emit:
        splice.out[0]
        splice.out[1]
}
