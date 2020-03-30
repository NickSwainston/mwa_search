nextflow.preview.dsl = 2
include beamform      from './beamform_module'
include pulsar_search from './pulsar_search_module'
include classifier    from './classifier_module'

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

params.dm_min = 1
params.dm_max = 250


pointings = Channel
    .from(params.pointings.split(","))
    .collect()
    .flatten()
    .collate( 15 )


workflow {
    pre_beamform()
    beamform( pre_beamform.out[0],\
              pre_beamform.out[1],\
              pre_beamform.out[2],\
              pointings )
    pulsar_search( beamform.out[0],\
                   beamform.out[1] )
    classifier( pulsar_search.out[2] )
}