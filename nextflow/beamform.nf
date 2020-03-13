nextflow.preview.dsl = 2
include beamform_wf from './beamform_module'

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


pointings = Channel
    .from(params.pointings.split(","))
    .collect()
    .flatten()
    .collate( 15 )
    //.view()


// Handling begin and end times
process get_all_beg_end {
    when:
    params.all == true

    """
    #!/usr/bin/env python3

    from mwa_metadb_utils import obs_max_min

    beg, end = obs_max_min(${params.obsid})
    print("{},{}".format(beg, end), end="")
    """
}


workflow {
    get_all_beg_end()
    if( params.all )
        beamform_wf( get_all_beg_end.out.map{ it.split(",") }.flatten().collect(), pointings )
    else
        beamform_wf( Channel.from( params.begin, params.end ).collect() , pointings )
}