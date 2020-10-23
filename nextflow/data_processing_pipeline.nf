#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

// The following allows * to perform a cartesian product on lists
class CartesianCategory {
    static Iterable multiply(Iterable a, Iterable b) {
        assert [a,b].every { it != null }
        def (m,n) = [a.size(),b.size()]
        (0..<(m*n)).inject([]) { prod, i -> prod << [a[i.intdiv(n)], b[i%n]].flatten() }
    }
}
Iterable.metaClass.mixin CartesianCategory

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

process make_yamls {
    input:
    tuple val(begin), val(end)
    tuple val(pointing), val(pulsar)

    output:
    file "*make_pulsar_yaml.yaml"

    """
    make_pulsar_yaml.py -o $params.obsid -O $params.calid --obs_beg $begin --obs_end $end --pointing ${pointing.join(" ")} --psr ${pulsar.join(" ")}\
    --mwa_search $params.mwa_search_version --vcstools $params.vcstools_version --label make_pulsar_yaml
    """
}


process pulsar_prepfold_cmd_make {
    input:
    each file(yaml_file)
    val label

    output:
    file "*[sh,edited_ephemeris.eph]"
    file "*${label}.yaml"
    // ephemeris files are formatted in the same way as the bash files

    """
    prepfold_cmd_make.py --yaml $yaml_file --label prepfold_cmd_make_${label}
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
    input:
    file yaml_and_pfd

    output:
    file "*yaml"
    file "*pfd*" includeInputs true

    """
    find_best_pointing.py --pfds *pfd* --yamls *yaml
    """
}


process decide_detections {
    input:
    file pfds_yamls

    output:
    file "*pfd*" includeInputs true
    file "*post_fold_filter.yaml"

    """
    post_fold_filter.py --yamls *yaml --pfds *pfd* --label post_fold_filter
    """
}

process polarimetry_call { //calls polarimetry which returns the apporpriate cmd file(s)
    input:
    tuple file (yaml) file(fits)
    val label

    output:
    tuple "*cmds.txt"
    tuple "*${label}.yaml", "*fits" includeInputs true

    """
    pulsar_polarimetry.py --yaml $yaml --fits $fits --label $label
    """
}

process fits_to_ar_and_back {
    label 'cpu'
    time "14400s"
    errorStrategy 'retry'
    maxRetries 1

    input:
    file fits
    file dspsr_fold_cmd
    file ar_to_fits_cmd

    output:
    file "*.newfits" //output fits file has .newfits extension

    container = "file:///${params.containerDir}/dspsr/dspsr.sif"

    """
    bash $dspsr_fold_cmd
    bash $ar_to_fits_cmd
    """
}

process baseline_removal {
    label 'cpu'
    time "600s"
    errorStrategy 'retry'
    maxRetries 1

    input:
    file newfits
    file debase_cmd

    beforeScript "module use ${params.presto_module_dir}; module load gcc/8.3.0; module load psrsalsa/master"

    output:
    file "*.debase.gg"

    """
    bash $debase_cmd
    """
}

process rm_synthesis {
    label 'cpu'
    time "3600s"
    errorStrategy 'retry'
    maxRetries 1

    input:
    file debased_fits
    file rmsynth_cmd

    beforeScript "module use ${params.presto_module_dir}; module load gcc/8.3.0; module load psrsalsa/master"

    output:
    file "*_plot.ps"
    file "*_map.ps"
    file "*.RMtable"

    """
    bash $rmsynth_cmd
    """
}

process defaraday_rotate {
    label 'cpu'
    time "1800s"
    errorStrategy 'retry'
    maxRetries 1

    input:
    file debased_fits
    file defarad_cmd

    beforeScript "module use ${params.presto_module_dir}; module load gcc/8.3.0; module load psrsalsa/master"

    output:
    file "*_profile.ps"
    file "*_polarimetry_profile.ps"
    file "*.paswing"

    """
    bash $defarad_cmd
    """
}

process rvm_fit {
    label 'cpu'
    time "3600s"
    errorStrategy 'retry'
    maxRetries 1

    input:
    file paswing
    file rvm_fit_cmd

    beforeScript "module use ${params.presto_module_dir}; module load gcc/8.3.0; module load psrsalsa/master"

    output:
    file "*_chigrid.ps"
    file "*.out"

    """
    bash $rvm_fit_cmd
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
        pulsar_prepfold_cmd_make( yaml_files.flatten(),\
                                  Channel.from("initial_fold") )
        // Run the bash file
        pulsar_prepfold_run( // Work out pointings from the file names
                             pulsar_prepfold_cmd_make.out[0].\
                             map{ it -> [it.flatten().findAll{ it != null }[-1].baseName.split("_J")[0].split("prepfold_cmd_${params.obsid}_")[1], it ] }.groupTuple().\
                             // Group fits files by bash files with same pointings
                             concat( fits_files ).groupTuple( size: 2, remainder: false ).\
                             // Then split them into a line per pulsar and format it to account for the optional .eph file
                             map{ it -> it[1][0] * it[1][1] }.flatMap().map{ it -> [it.init(), it.last()]} )
        //if ( (params.search_radius - fwhm / 2) > (fwhm * 0.6) ){
            // If more than one loop of beams per source,
        //}
        // Run through the classfier
        classifier( pulsar_prepfold_run.out.flatten().collate( 120 ) )
        // Find the best detection for each pulsar
        best_detection( // Pair the classifier output with their yaml file
                        classifier.out[0].flatten().map{ it -> [ it.baseName.split("_b")[0], it ]}.groupTuple().concat(
                        pulsar_prepfold_cmd_make.out[1].flatten().map{ it -> [ it.baseName.split("prepfold_cmd_make")[0], it ]}.groupTuple()).\
                        // Group by pulsar
                        map{ it -> [ it[0].split("_")[-1], it[1] ]}.groupTuple( size: 2, remainder: false  ).\
                        map{ it -> it[1][0] + it[1][1] } )
    emit:
        //classifier.out[0] //classifier files
        best_detection.out[0] //yaml files
        best_detection.out[1] //pfd files
}

workflow post_fold{
    take:
        yaml_files
        init_pfd_files
        fits_files
    main:
        pulsar_prepfold_cmd_make( yaml_files,\
                                  Channel.from("post_fold") )
        pulsar_prepfold_run(// Work out pointings from the file names
                             pulsar_prepfold_cmd_make.out[0].flatten().\
                             map{ it -> [it.flatten().findAll{ it != null }[-1].baseName.split("_J")[0].split("prepfold_cmd_${params.obsid}_")[1], it ] }.groupTuple().\
                             // Group fits files by bash files with same pointings
                             concat( fits_files ).groupTuple( size: 2, remainder: false ).\
                             // Then split them into a line per pulsar and format it to account for the optional .eph file
                             map{ it -> it[1][0] * it[1][1] }.flatMap().map{ it -> [it.init(), it.last()]} )
         //figures out which bestprofs are detections and updates yaml file
        decide_detections(// Group the one yaml file and the ~4 groups of pfd files for each pulsar
                          pulsar_prepfold_cmd_make.out[1].flatten().map{ it -> [ it.baseName.split("_prepfold")[0], it ] }.concat(
                          init_pfd_files.concat(pulsar_prepfold_run.out).flatten().map{ it -> [ it.baseName.split("_b")[0], it ] } ).groupTuple().
                          map{ it -> it[1]} )

    emit:
    //detections (pfds)
    decide_detections.out[0]
    //yaml files
    decide_detections.out[1]
}

workflow polarimetry{
    take:
        yaml_fits_tuple
    main:
        //polarimetry 1 - convert fits to archive, back to fits, baseline removal, RM synthesis
        polarimetry_call(yaml_fits_tuple, "polarimetry_one")
        fits_to_ar_and_back(polarimetry_call.out[1][1],  polarimetry_call.out[0].filter(~/$"dspsr_fold_cmds.sh"/), polarimetry_call.out[0].filter( ~/$"to_fits_cmds.sh"/))
        baseline_removal(fits_to_ar_and_back.out[0], polarimetry_call.out[0].filter(~/$"debase_cmds.sh"/))
        rm_synthesis(baseline_removal.out[0], polarimetry_call.out.filter[0].filter(~/$"initial_rm_synthesis_cmds.sh"/))
        //polarimetry 2 - Final RM synthesis
        polarimetry_call(polarimetry_call.out[1], "polarimetry_two")
        rm_synthesis(baseline_removal.out[0], polarimetry_call.out[0].filter(~/$"final_rm_synthesis_cmds.sh"/))
        //polarimetry 3 - Defaraday rotation
        polarimetry_call(polarimetry_call.out[1], "polarimetry_three")
        defaraday_rotate(baseline_removal.out[0], polarimetry_call.out[0].filter(~/$"defarad_cmds.sh"/))
        // polarimetry 4 - Initial RVM fitting
        polarimetry_call(polarimetry_call.out[1], "polarimetry_four")
        rvm_fit(defaraday_rotate.out[2], polarimetry_call.out[0].filter(~/$"initial_rvmfit_cmds.sh"/))
        // polarimetry 5 - Final RVM fitting
        polarimetry_call(polarimetry_call.out[1], "polarimetry_five")
        rvm_fit(defaraday_rotate.out[2], polarimetry_call.out[0].filter(~/$"final_rvmfit_cmds.sh"/))
        // polarimetry 6 - Reading the final RVM fit
        polarimetry_call(polarimetry_call.out[1], "polarimetry_six")

    emit:
        polarimetry_call.out[1][0]  //yaml file
        fits_to_ar_and_back.out[0]  // converted fits file
        baseline_removal.out[0]     // baseline femoved fits file
        rm_synthesis.out[0]         // RMsynth plot
        rm_synthesis.out[1]         // RMsynth map plot
        rm_synthesis.out[2]         // RMtable
        defaraday_rotate.out[1]     // polarimetry profile
        defaraday_rotate.out[2]     // paswing
        rvm_fit.out[0]              // Chi grid
        rvm_fit.out[1]              // Fit information

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
              find_pointings.out.splitCsv(skip: 7, limit: 1)).collect().flatten().unique().filter{ it != " " }.collate( params.max_pointings ) )
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
               initial_fold.out[0],\
               //initial pfd files
               initial_fold.out[1],\
               // fits files
               beamform.out[3].concat(beamform_ipfb.out[3]) )

    // Polarimetry
    polarimetry(//yaml + fits file tuple
                post_fold.out[1].map{ it -> [it.flatten().findAll{ it != null }[-1].split("_J").split("post_fold_filter_${params.obsid}_")[1], it]}.groupTuple().\
                                concat( beamform.out[3].concat(beamform_ipfb.out[3]) ).groupTuple( size: 2, remainder: false ).\
                                map{ it -> it[1][0] * it[1][1] }.flatMap().map{ it -> [it.init(), it.last()]} )


    // Perform a search on all candidates (not known pulsars)
    // if pointing in fits file name is in pulsar search pointing list
    pulsar_search( find_pointings.out.splitCsv(skip: 5, limit: 1).flatten().merge(find_pointings.out.splitCsv(skip: 4, limit: 1).flatten()).\
                   concat(beamform.out[3]).groupTuple( size: 2, remainder: false ).map { it -> [ "Blind_${params.obsid}_${it[0]}".toString(), it[1][1] ] } )
    classifier( pulsar_search.out[1].flatten().collate( 600 ) )
    // Perform a single pulse search on all single pulse candidates
    single_pulse_search( find_pointings.out.splitCsv(skip: 7, limit: 1).flatten().merge(find_pointings.out.splitCsv(skip: 6, limit: 1).flatten()).\
                         concat(beamform.out[3]).groupTuple( size: 2, remainder: false ).map { it -> it[1] } )
}