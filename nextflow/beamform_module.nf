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
params.didir = "${params.basedir}/${params.obsid}/cal/${params.calid}/rts"
params.channels = null


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


process splice {
    publishDir "${params.basedir}/${params.obsid}/pointings/${unspliced[0].baseName.split("_")[2]}_${unspliced[0].baseName.split("_")[3]}", mode: 'move'

    label 'cpu'
    time '1h'

    input:
    val chan
    each file(unspliced)

    output:
    file "${params.obsid}*fits"
    """
    module use /group/mwa/software/modulefiles
    module load mwa_search
    module load vcstools
    splice_wrapper.py -o ${params.obsid} -c ${chan.join(" ")}
    """
}


workflow beamform_wf {
    take: 
        obs_beg_end
        pointings
        channels
    main:
        ensure_metafits()
        gps_to_utc( obs_beg_end )
        make_beam( channels | flatten() | merge(range),\
                   gps_to_utc.out,\
                   Channel.from(params.pointings.split(",")).collect().flatten().collate( 15 ),\
                   obs_beg_end )
        splice( channels,\
                make_beam.out | flatten() | map { it -> [it.baseName.split("ch")[0], it ] } | groupTuple() | map { it -> it[1] } )
}
