#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

params.obsid = null
params.calid = null

params.begin = null
params.end = null
params.all = false

params.search_radius = 0.02
params.fwhm_deg = null

params.vcstools_version = 'master'
params.mwa_search_version = 'master'

params.didir = "${params.scratch_basedir}/${params.obsid}/cal/${params.calid}/rts"
params.publish_fits = false
params.publish_fits_scratch = false
params.publish_all_classifer_cands = false

params.out_dir = "${params.search_dir}/${params.obsid}_candidates"

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
             |
             |Optional arguments:
             |  --search_radius
             |              The radius to search (create beams within) in degrees to account for ionosphere.
             |              [default: 0.02 degrees]
             |  --publish_fits
             |              Publish to the fits directory (/group on Galaxy). Use this instead
             |              of --publish_fits_scratch
             |  --publish_fits_scratch
             |              Publish to the scratch fits directory (/astro on Galaxy). Include
             |              this option.
             |  --vcstools_version
             |              The vcstools module version to use [default: master]
             |  --mwa_search_version
             |              The mwa_search module bersion to use [default: master]
             |  --no_combined_check
             |              Don't check if all the combined files are available [default: false]
             |  --out_dir   Where the search candidates will be output
             |              [default: ${params.out_dir}]
             |  -w          The Nextflow work directory. Delete the directory once the processs
             |              is finished [default: ${workDir}]""".stripMargin()
    println(help)
    exit(0)
}


process fwhm_calc {
    input:
    val channels

    output:
    file "${params.obsid}_fwhm.txt"

    """
    #!/usr/bin/env python3

    from vcstools.metadb_utils import get_obs_array_phase
    from mwa_search.obs_tools import calc_ta_fwhm
    import csv

    if "${params.fwhm_deg}" == "null":
        oap = get_obs_array_phase(${params.obsid})
        centrefreq = 1.28 * float(${channels[0]} + ${channels[-1]}) / 2.
        fwhm = calc_ta_fwhm(centrefreq, array_phase=oap)
    else:
        fwhm = ${params.fwhm_deg}

    with open("${params.obsid}_fwhm.txt", "w") as outfile:
        spamwriter = csv.writer(outfile, delimiter=',')
        spamwriter.writerow([fwhm])
    """
}

process find_pointings {
    input:
    tuple val(begin), val(end)
    val fwhm

    output:
    file "${params.obsid}_fov_sources.csv"

    """
    pulsars_in_fov.py -o $params.obsid -b $begin -e $end --fwhm $fwhm --search_radius ${params.search_radius}
    """
}

include { pre_beamform; beamform; beamform_ipfb } from './beamform_module'
include { pulsar_search; single_pulse_search } from './pulsar_search_module'
include { classifier } from './classifier_module'

workflow {
    pre_beamform()
    fwhm_calc( pre_beamform.out[1] )
    find_pointings( pre_beamform.out[0],
                    fwhm_calc.out.splitCsv().flatten() )
    beamform( pre_beamform.out[0],\
              pre_beamform.out[1],\
              pre_beamform.out[2],\
              //Grab the pointings for slow pulsars and single pulses
              find_pointings.out.splitCsv(skip: 1, limit: 1).concat(\
              find_pointings.out.splitCsv(skip: 5, limit: 1),\
              find_pointings.out.splitCsv(skip: 7, limit: 1)).collect().flatten().unique().filter{ it != " " }.collate( params.max_pointings ) )
    beamform_ipfb( pre_beamform.out[0],\
                   pre_beamform.out[1],\
                   pre_beamform.out[2],\
                   //Grab the pointings for slow pulsars and single pulses
                   find_pointings.out.splitCsv(skip: 3, limit: 1) )

    // Perform a search on all candidates (not known pulsars)
    // if pointing in fits file name is in pulsar search pointing list
    pulsar_search( find_pointings.out.splitCsv(skip: 5, limit: 1).flatten().merge(find_pointings.out.splitCsv(skip: 4, limit: 1).flatten()).\
                   concat(beamform.out[3]).groupTuple( size: 2, remainder: false ).map { it -> [ "Blind_${params.obsid}_${it[0]}".toString(), it[1][1] ] } )
    classifier( pulsar_search.out[1].flatten().collate( 600 ) )
    // Perform a single pulse search on all single pulse candidates
    single_pulse_search( find_pointings.out.splitCsv(skip: 7, limit: 1).flatten().merge(find_pointings.out.splitCsv(skip: 6, limit: 1).flatten()).\
                         concat(beamform.out[3]).groupTuple( size: 2, remainder: false ).map { it -> it[1] } )
}