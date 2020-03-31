nextflow.preview.dsl = 2
include { pre_beamform; beamform } from './beamform_module'

params.obsid = null
params.pointings = null
params.calid = null

params.begin = null
params.end = null
params.all = false

params.summed = false
params.vcstools_version = 'master'
params.mwa_search_version = 'master'
params.channels = null

params.basedir = '/group/mwaops/vcs'
params.scratch_basedir = '/astro/mwaops/vcs'
params.didir = "${params.basedir}/${params.obsid}/cal/${params.calid}/rts"
params.publish_fits = true
params.publish_fits_scratch = false

params.no_combined_check = false

pointings = Channel
    .from(params.pointings.split(","))
    .collect()
    .flatten()
    .collate( 15 )
    //.view()


workflow {
    pre_beamform()
    beamform( pre_beamform.out[0],\
              pre_beamform.out[1],\
              pre_beamform.out[2],\
              pointings )
}