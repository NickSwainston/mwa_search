nextflow.preview.dsl = 2
include { pre_beamform; beamform } from './beamform_module'
include pulsar_search from './pulsar_search_module'
include classifier    from './classifier_module'

params.obsid = null
params.calid = null
params.pointings = null
params.pointing_file = null

params.begin = null
params.end = null
params.all = false

params.summed = true
params.vcstools_version = 'master'
params.mwa_search_version = 'master'

params.basedir = '/group/mwaops/vcs'
params.didir = "${params.basedir}/${params.obsid}/cal/${params.calid}/rts"
params.channels = null
params.out_dir = "${params.obsid}_candidates"

params.dm_min = 1
params.dm_max = 250

params.no_combined_check = false

if ( params.pointing_file ) {
    pointings = Channel
        .fromPath(params.pointing_file)
        .splitCsv()
        .collect()
        .flatten()
        .collate( params.max_pointings )
}
else if ( params.pointings ) {
    pointings = Channel
        .from(params.pointings.split(","))
        .collect()
        .flatten()
        .collate( params.max_pointings )
}
else {
    println "No pointings given. Either use --pointing_file or --pointings. Exiting"
    exit(1)
}

workflow {
    pre_beamform()
    beamform( pre_beamform.out[0],\
              pre_beamform.out[1],\
              pre_beamform.out[2],\
              pointings )
    pulsar_search( beamform.out[0],\
                   beamform.out[1] )
    classifier( pulsar_search.out[2].flatten().collate( 120 ) )
    publish:
        classifier.out to: params.out_dir
        pulsar_search.out to: params.out_dir, pattern: "*singlepulse*"
}
