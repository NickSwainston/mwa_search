#!/usr/bin/env nextflow

nextflow.preview.dsl = 2
include { pre_beamform; beamform; beamform_ipfb; get_beg_end } from './beamform_module'
include { pulsar_search; single_pulse_search } from './pulsar_search_module'
include classifier from './classifier_module'

params.obsid = null
params.calid = null

params.begin = null
params.end = null
params.all = false

params.vcstools_version = 'master'
params.mwa_search_version = 'master'

params.didir = "${params.basedir}/${params.obsid}/cal/${params.calid}/rts"
params.publish_fits = false
params.publish_fits_scratch = true

params.out_dir = "${params.basedir}/${params.obsid}/${params.obsid}_candidates"

params.no_combined_check = false

params.help = false
if ( params.help ) {
    help = """beamform_fov_sources.nf: A pipeline that will beamform on all pulsars in the FOV
             |                        and perform a search on all pulsar candidates.
             |Required argurments:
             |  --obsid     Observation ID you want to process [no default]
             |  --calid     Observation ID of calibrator you want to process [no default]
             |  --begin     First GPS time to process [no default]
             |  --end       Last GPS time to process [no default]
             |  --all       Use entire observation span. Use instead of -b & -e. [default: false]
             |  --publish_fits_scratch
             |              Publish to the scratch fits directory (/astro on Galaxy). Include
             |              this option.
             |
             |Optional arguments:
             |  --publish_fits
             |              Publish to the fits directory (/group on Galaxy). Use this instead
             |              of --publish_fits_scratch
             |  --vcstools_version
             |              The vcstools module version to use [default: master]
             |  --mwa_search_version
             |              The mwa_search module bersion to use [default: master]
             |  --no_combined_check
             |              Don't check if all the combined files are available [default: false]
             |  --out_dir   Where the search candidates will be output
             |              [default: /group/mwaops/vcs/<obsid>/<obsid>_candidates]
             |  -w          The Nextflow work directory. Delete the directory once the processs
             |              is finished [default: ${workDir}]""".stripMargin()
    println(help)
    exit(0)
}



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
                   concat(beamform.out[3]).groupTuple( size: 2, remainder: false ).map { it -> [ 'Blind_'+it[0], it[1][1] ] },\
                   pre_beamform.out[1] )
    classifier( pulsar_search.out[2].flatten().collate( 600 ) )
    // Perform a single pulse search on all single pulse candidates
    single_pulse_search( find_pointings.out.splitCsv(skip: 7, limit: 1).flatten().merge(find_pointings.out.splitCsv(skip: 6, limit: 1).flatten()).\
                         concat(beamform.out[3]).groupTuple( size: 2, remainder: false ).map { it -> it[1] },\
                         pre_beamform.out[1] )
    publish:
        classifier.out to: params.out_dir
        pulsar_search.out to: params.out_dir, pattern: "*singlepulse*"
        single_pulse_search.out to: params.out_dir
}
