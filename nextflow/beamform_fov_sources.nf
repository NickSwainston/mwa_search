nextflow.preview.dsl = 2
include { pre_beamform; beamform; beamform_ipfb; get_beg_end } from './beamform_module'
include { pulsar_search; single_pulse_search } from './pulsar_search_module'
include classifier from './classifier_module'

params.obsid = null
params.pointings = null
params.calid = null

params.begin = null
params.end = null
params.all = false

params.summed = false
params.vcstools_version = 'master'

params.basedir = '/group/mwaops/vcs'
params.scratch_basedir = '/astro/mwaops/vcs'
params.didir = "${params.basedir}/${params.obsid}/cal/${params.calid}/rts"
params.publish_fits = true

params.out_dir = "${params.basedir}/${params.obsid}/${params.obsid}_candidates"

params.no_combined_check = false



process find_pointings {
    input:
    tuple val(begin), val(end)

    output:
    file "${params.obsid}_fov_sources.csv"

    """
    pulsars_in_fov.py -o $params.obsid -b $begin -e $end
    """
}

workflow {
    get_beg_end()
    find_pointings( get_beg_end.out.map{ it.split(",") }.flatten().collect() )
    pre_beamform()
    beamform( pre_beamform.out[0],\
              pre_beamform.out[1],\
              pre_beamform.out[2],\
              //Grab the pointings for slow pulsars and single pulses
              find_pointings.out.splitCsv(skip: 1, limit: 1).mix(\
              find_pointings.out.splitCsv(skip: 5, limit: 1),\
              find_pointings.out.splitCsv(skip: 7, limit: 1)).collect().flatten().unique().collate( 15 ) )
    beamform_ipfb( pre_beamform.out[0],\
                   pre_beamform.out[1],\
                   pre_beamform.out[2],\
                   //Grab the pointings for slow pulsars and single pulses
                   find_pointings.out.splitCsv(skip: 3, limit: 1) )
    // Perform a search on all candidates (not known pulsars)
    // if pointing in fits file name is in pulsar search pointing list
    pulsar_search( find_pointings.out.splitCsv(skip: 5, limit: 1).flatten().merge(find_pointings.out.splitCsv(skip: 4, limit: 1).flatten()).\
                   concat(beamform.out[2]).groupTuple( size: 2, remainder: false ).map { it -> [ 'Blind_'+it[0], it[1][1] ] } )
    classifier( pulsar_search.out[2].flatten().collate( 120 ) )
    // Perform a single pulse search on all single pulse candidates
    single_pulse_search( find_pointings.out.splitCsv(skip: 7, limit: 1).flatten().merge(find_pointings.out.splitCsv(skip: 6, limit: 1).flatten()).\
                         concat(beamform.out[2]).groupTuple( size: 2, remainder: false ).map { it -> it[1] } )
    publish:
        classifier.out to: params.out_dir
        pulsar_search.out to: params.out_dir, pattern: "*singlepulse*"
        single_pulse_search.out to: params.out_dir
}
