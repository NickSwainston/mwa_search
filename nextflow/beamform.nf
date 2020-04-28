nextflow.preview.dsl = 2
include { pre_beamform; beamform; beamform_ipfb } from './beamform_module'

params.obsid = null
params.pointings = null
params.pointing_file = null
params.calid = null

params.begin = null
params.end = null
params.all = false

params.summed = false
params.ipfb = false
params.vcstools_version = 'master'
params.mwa_search_version = 'master'
params.channels = null

params.basedir = '/group/mwaops/vcs'
params.scratch_basedir = '/astro/mwaops/vcs'
params.didir = "${params.basedir}/${params.obsid}/cal/${params.calid}/rts"
params.publish_fits = true
params.publish_fits_scratch = false

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
    if ( params.ipfb ) {
        beamform_ipfb( pre_beamform.out[0],\
                       pre_beamform.out[1],\
                       pre_beamform.out[2],\
                       pointings )
    }
    else {
        beamform( pre_beamform.out[0],\
                  pre_beamform.out[1],\
                  pre_beamform.out[2],\
                  pointings )
    }
}