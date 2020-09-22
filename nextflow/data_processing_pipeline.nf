#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

params.obsid = null
params.calid = null

params.begin = null
params.end = null
params.all = false

params.search_radius = 0.02

params.vcstools_version = 'master'
params.mwa_search_version = 'master'

params.didir = "${params.scratch_basedir}/${params.obsid}/cal/${params.calid}/rts"
params.publish_fits = false
params.publish_fits_scratch = false

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

// Work out some estimated job times
if ( "$HOSTNAME".startsWith("farnarkle") ) {
    // In seconds
    search_dd_fft_acc_dur = obs_length * 5.0
    prepfold_dur = obs_length * 2.0
    presto_python_load = "module use ${params.presto_module_dir}; module load presto/${params.presto_module}; module load python/2.7.14; module load matplotlib/2.2.2-python-2.7.14"
}
else {
    search_dd_fft_acc_dur = 14400
    prepfold_dur = 7200
    presto_python_load = ""
}

process fwhm_calc {
    input:
    val channels

    output:
    file "${params.obsid}_fwhm.txt"

    """
    #!/usr/bin/env python3

    from mwa_metadb_utils import get_obs_array_phase
    from mwa_search.obs_tools import calc_ta_fwhm
    import csv

    oap = get_obs_array_phase(${params.obsid})
    centrefreq = 1.28 * float(${channels[0]} + ${channels[-1]}) / 2.
    fwhm = calc_ta_fwhm(centrefreq, array_phase=oap)

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

process make_yamls {
    input:
    tuple val(begin), val(end)
    tuple val(pointing), val(pulsar)

    output:
    file "*.yaml"

    """
    make_pulsar_yaml.py -o $params.obsid -O $params.calid --obs_beg $begin --obs_end $end --pointing ${pointing.join(" ")} --psr ${pulsar.join(" ")}\
    --mwa_search $params.mwa_search_version --vcstools $params.vcstools_version --label make_pulsar_yaml -d ./
    """
}


process pulsar_prepfold_cmd_make {
    input:
    file yaml_file

    output:
    file "*[sh,edited_ephemeris.eph]"
    file "*prep_cmd_make.yaml"
    // ephemeris files are formatted in the same way as the bash files

    """
    prepfold_cmd_make.py --yaml $yaml_file --label prep_cmd_make
    """
}


process pulsar_prepfold_run {
    label 'cpu'
    time "${prepfold_dur}s"
    errorStrategy 'retry'
    maxRetries 1

    input:
    tuple file(prepfold_cmd_file_and_eph), file(fits)

    output:
    file "*pfd*"

    if ( "$HOSTNAME".startsWith("farnarkle") ) {
        beforeScript "module use ${params.presto_module_dir}; module load presto/${params.presto_module}"
    }
    else if ( "$HOSTNAME".startsWith("x86") || "$HOSTNAME".startsWith("garrawarla") || "$HOSTNAME".startsWith("galaxy") ) {
        container = "file:///${params.containerDir}/presto/presto.sif"
    }
    else {
        container = "nickswainston/presto:realfft_docker"
    }
    """
    bash *sh
    """
}

process best_detection {
    label 'cpu'
    time '10m'
    errorStrategy 'retry'
    maxRetries 1

    input:
    file yaml_and_pfd

    output:
    file "*yaml"

    """
    hsdhjsj
    """
}


include { pre_beamform; beamform; beamform_ipfb } from './beamform_module'
include { pulsar_search; single_pulse_search } from './pulsar_search_module'
include { classifier } from './classifier_module'

workflow initial_fold {
    take:
        yaml_files
        fits_files
    main:
        // Create a bash file of the prepfold commands required
        pulsar_prepfold_cmd_make( yaml_files.flatten() )
        // Run the bash file
        pulsar_prepfold_run( // Work out pointings from the file names
                             pulsar_prepfold_cmd_make.out[0].\
                             map{ it -> [it.flatten().findAll{ it != null }[-1].baseName.split("_J")[0].split("prepfold_cmd_${params.obsid}_")[1], it ] }.groupTuple().\
                             // Group fits files by bash files with same pointings
                             concat( fits_files ).groupTuple( size: 2, remainder: false ).view().map{ it -> it[1] } )
        //if ( (params.search_radius - fwhm / 2) > (fwhm * 0.6) ){
            // If more than one loop of beams per source,
        //}
        // Run through the classfier
        classifier( pulsar_prepfold_run.out.flatten().collate( 120 ) )
        // Find the best detection for each pulsar
        best_detection( // Pair the classifier output witht their yaml file
                        classifier.out[0].flatten().map{ it -> [ it.baseName.split("_b")[0], it ]}.groupTuple().concat(
                        yaml_files.flatten().map{ it -> [ it.baseName.split("make_pulsar_yaml")[0], it ]}.groupTuple()).\
                        // Group by pulsar
                        map{ it -> [ it[0].split("_")[-1], it[1] ]}.groupTuple( size: 2, remainder: false  ).\
                        map{ it -> it[1][0] + it[1][1] } )
    emit:
        //classifier.out[0] //classifier files
        pulsar_prepfold_cmd_make.out[1] //yaml files
}


process decide_detections {
    input:
    file pfds
    file yamls

    output
    file *pfd*
    file *post_detection.yaml

    """
    post_fold_filter.py --yamls $yamls --pfds $pfds --label post_fold_filter
    """

    output:
    file "*post_fold_filter.yaml"
    file "*pfd*"

}


workflow post_fold{
    take:
        yaml_files
        fits_files
    main:
        pulsar_prepfold_cmd_make(yaml_files)
        pulsar_prepfold_run(prepfold_cmd_make.out, fits_files) //TODO: group the right prepfold commands with fits files
        decide_detections(pulsar_prepfold_run.out, pulsar_prepfold_cmd_make.out[1]) //figures out which bestprofs are detections and updates yaml file

    emit:
    //detections
    decide_detections.out[0]
    //yaml files
    decide_detections.out[1]
}


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
              find_pointings.out.splitCsv(skip: 7, limit: 1)).collect().flatten().unique().filter{ it != " " }.collate( params.max_pointings ).view() )
    beamform_ipfb( pre_beamform.out[0],\
                   pre_beamform.out[1],\
                   pre_beamform.out[2],\
                   //Grab the pointings for slow pulsars and single pulses
                   find_pointings.out.splitCsv(skip: 3, limit: 1) )

    // Make a yaml_file with all necessary info for each pointing
    make_yamls( pre_beamform.out[0],\
                find_pointings.out.splitCsv(skip: 1, limit: 1).concat( find_pointings.out.splitCsv(skip: 3, limit: 1) ).collect().map{ it -> [it] }.concat(\
                find_pointings.out.splitCsv(skip: 0, limit: 1).concat( find_pointings.out.splitCsv(skip: 2, limit: 1) ).collect().map{ it -> [it] }).collect() )

    // Perform processing pipeline on all known pulsars
    initial_fold( // yaml files
                  make_yamls.out,\
                  // fits files
                  beamform.out[3].concat(beamform_ipfb.out[3]) )

    //post_fold()
    post_fold( //yaml files
               initial_fold.out,\
               // fits files
               beamform.out[3].concat(beamform_ipfb.out[3]) )

    // Perform a search on all candidates (not known pulsars)
    // if pointing in fits file name is in pulsar search pointing list
    pulsar_search( find_pointings.out.splitCsv(skip: 5, limit: 1).flatten().merge(find_pointings.out.splitCsv(skip: 4, limit: 1).flatten()).\
                   concat(beamform.out[3]).groupTuple( size: 2, remainder: false ).map { it -> [ "Blind_${params.obsid}_${it[0]}".toString(), it[1][1] ] } )
    classifier( pulsar_search.out[1].flatten().collate( 600 ) )
    // Perform a single pulse search on all single pulse candidates
    single_pulse_search( find_pointings.out.splitCsv(skip: 7, limit: 1).flatten().merge(find_pointings.out.splitCsv(skip: 6, limit: 1).flatten()).\
                         concat(beamform.out[3]).groupTuple( size: 2, remainder: false ).map { it -> it[1] } )
}