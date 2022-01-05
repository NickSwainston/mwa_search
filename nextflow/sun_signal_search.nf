#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.obsid = null
params.pointings = null
params.fitsdir = "/group/mwavcs/vcs/${params.obsid}/pointings"
params.fits_files
params.out_dir = "${params.search_dir}/${params.obsid}_candidates"
params.dm_min = 0
params.dm_max = 0.02
params.dm_min_step = 0.02

params.scratch = false
params.fits_file_dir = false

params.cand = "Blind"
params.sp = false
params.publish_all_prepfold = true

//Defaults for the accelsearch command
params.nharm = 16 // number of harmonics to search
params.min_period = 0.001 // min period to search for in sec (ANTF min = 0.0013)
params.max_period = 30 // max period to search for in sec  (ANTF max = 23.5)
params.zmax = 0

//Some math for the accelsearch command
//convert to freq
min_freq = 1 / params.max_period
max_freq = 1 / params.min_period
//adjust the freq to include the harmonics
min_f_harm = min_freq
max_f_harm = max_freq * params.nharm

// Work out some estimated job times
if ( "$HOSTNAME".startsWith("farnarkle") ) {
    presto_python_load = "module use ${params.presto_module_dir}; module load presto/${params.presto_module}; module load python/2.7.14; module load matplotlib/2.2.2-python-2.7.14"
}
else {
    presto_python_load = ""
}

pointing = Channel.from( params.pointings )
if ( params.fits_file_dir ) {
    fits_files = Channel.fromPath( "${params.fits_file_dir}/${params.obsid}_*.fits", checkIfExists: true )
    nfiles = new File("${params.fits_file_dir}").listFiles().findAll { it.name ==~ /.*1*fits/ }.size()
}
else if ( params.fits_files ) {
    fits_files = Channel.fromPath( "${params.fits_files}", checkIfExists: true )
    nfiles = new File("${params.fits_files}").listFiles().findAll { it.name ==~ /.*1*fits/ }.size()
    name_fits_files = Channel.from( fits_files.map{ it -> [ params.cand + '_' + it.getBaseName().split("/")[-1].split("_ch")[0], it ] } )
}
else {
    if ( params.scratch ) {
        fits_files = Channel.fromPath( "${params.scratch_basedir}/${params.obsid}/dpp_pointings/${params.pointings}/${params.obsid}_*.fits", checkIfExists: true )
        nfiles = new File("${params.scratch_basedir}/${params.obsid}/dpp_pointings/${params.pointings}").listFiles().findAll { it.name ==~ /.*1*fits/ }.size()
    }
    else {
        fits_files = Channel.fromPath( "${params.basedir}/${params.obsid}/pointings/${params.pointings}/${params.obsid}_*.fits", checkIfExists: true )
        nfiles = new File("${params.basedir}/${params.obsid}/pointings/${params.pointings}").listFiles().findAll { it.name ==~ /.*1*fits/ }.size()
    }
}

//name_fits_files = Channel.from( fits_files.map{ it -> [ params.cand + '_' + it.getBaseName().split("/")[-1].split("_ch")[0], it ] } )

params.help = false
if ( params.help ) {
    help = """pulsar_search.nf: A pipeline perform a pulsar search on the input fits files.
             |                  The fits files must be in the format
             |                  <obsid>_<pointing>_ch<min_chan>-<max_chan>_00??.fits
             |Required argurments:
             |  --obsid     Observation ID you want to process [no default]
             |  --pointings The pointing to search in with the RA and Dec seperated
             |              by _ in the format HH:MM:SS_+DD:MM:SS
             |
             |Optional arguments:
             |  --cand      Candidate name to do a targeted search [default: Blind]
             |  --sp        Perform a single pulse search [default false]
             |  --fits_file_dir
             |              Directory containing the fits files. Use this if the fits files
             |              are not in the default directory :
             |              ${params.basedir}/<obsid>/pointings/${params.pointings}
             |  --scratch   Change the default directory to:
             |              ${params.scratch_basedir}/<obsid>/pointings/${params.pointings}
             |  --dm_min    Minimum DM to search over [default: 1]
             |  --dm_max    Maximum DM to search over [default: 250]
             |  --dm_min_step
             |              Minimum DM step size (Delta DM) [default: 0.1]
             |  --out_dir   Output directory for the candidates files
             |              [default: ${params.search_dir}/<obsid>_candidates]
             |  --mwa_search_version
             |              The mwa_search module bersion to use [default: master]
             |  -w          The Nextflow work directory. Delete the directory once the processs
             |              is finished [default: ${workDir}]""".stripMargin()
    println(help)
    exit(0)
}

process search_dd_fft_acc {
    label 'cpu'
    time '5h'
    //Will ignore errors for now because I have no idea why it dies sometimes
    errorStrategy { task.attempt > 1 ? 'ignore' : 'retry' }
    maxForks 800
    publishDir "${params.scratch_basedir}/${params.obsid}/pointings/${name.split("${params.obsid}_")[1]}", mode: 'copy'

    input:
    tuple val(name), val(dm_values), file(fits_files)

    output:
    tuple val(name), file("*ACCEL_${params.zmax}"), file("*.dat"), file("*.inf"), file("*SpS"), file('*.cand')

    if ( "$HOSTNAME".startsWith("farnarkle") ) {
        clusterOptions { "--export=NONE "}
        scratch '$JOBFS'
        beforeScript "module use ${params.module_dir}; module load presto/min_path"
    }
    else if ( "$HOSTNAME".startsWith("x86") ) {
        container = "file:///${params.containerDir}/presto/presto.sif"
    }
    else if ( "$HOSTNAME".startsWith("galaxy") ) {
        container = "file:///${params.containerDir}/presto/presto.sif"
    }
    else if ( "$HOSTNAME".startsWith("garrawarla") ) {
        clusterOptions { "--export=NONE " }
        scratch '/nvmetmp'
        container = "file:///${params.containerDir}/presto/presto.sif"
    }
    else {
        container = "nickswainston/presto:realfft_docker"
    }

    //numout=${(int)(obs_length*10000/Float.valueOf(dm_values[5]))}
    // -numout \${numout} 
    """
    echo "lowdm highdm dmstep ndms timeres downsamp"
    echo ${dm_values}
    printf "\\n#Dedispersing the time series at \$(date +"%Y-%m-%d_%H:%m:%S") --------------------------------------------\\n"
    prepdata -ncpus $task.cpus -dm ${dm_values[0]} -o ${name}_DM0.00 ${params.obsid}_*.fits
    printf "\\n#Performing the FFTs at \$(date +"%Y-%m-%d_%H:%m:%S") -----------------------------------------------------\\n"
    realfft *dat
    printf "\\n#Performing the periodic search at \$(date +"%Y-%m-%d_%H:%m:%S") ------------------------------------------\\n"
    for i in \$(ls *.dat); do
        accelsearch -ncpus $task.cpus -zmax ${params.zmax} -flo $min_f_harm -fhi $max_f_harm -numharm $params.nharm \${i%.dat}.fft
    done
    ${presto_python_load}
    single_pulse_search.py -p -m 0.5 -b *.dat
    cat *.singlepulse > ${name}_DM${dm_values[0]}-${dm_values[1]}.SpS
    printf "\\n#Finished at \$(date +"%Y-%m-%d_%H:%m:%S") ----------------------------------------------------------------\\n"
    """
}

process accelsift {
    label 'cpu'
    time '25m'
    errorStrategy 'retry'
    maxRetries 1

    input:
    tuple val(name), file(accel_inf_single_pulse)

    output:
    tuple val(name), file("cands_*greped.txt")

    if ( "$HOSTNAME".startsWith("farnarkle") || "$HOSTNAME".startsWith("x86") ||\
              "$HOSTNAME".startsWith("garrawarla") || "$HOSTNAME".startsWith("galaxy") ) {
        container = "file:///${params.containerDir}/presto/presto.sif"
    }
    else {
        container = "nickswainston/presto:realfft_docker"
    }
    """
    ACCEL_sift.py --file_name ${name} --min_num_DMs 1 --low_DM_cutoff 0
    if [ -f cands_${name}.txt ]; then
        grep ${name} cands_${name}.txt > cands_${name}_greped.txt
    else
        #No candidates so make an empty file
        touch cands_${name}_greped.txt
    fi
    """
}

process single_pulse_searcher {
    label 'cpu_large_mem'
    time '2h'
    stageInMode = 'copy'
    publishDir params.out_dir, mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(name), file(sps), file(fits)

    output:
    file "*pdf" optional true

    if ( "$HOSTNAME".startsWith("farnarkle") || "$HOSTNAME".startsWith("x86") ||\
         "$HOSTNAME".startsWith("garrawarla") || "$HOSTNAME".startsWith("galaxy") ) {
        container = "file:///${params.containerDir}/sps/sps.sif"
    }
    else {
        container = "nickswainston/sps"
    }
    """
    #-SNR_min 4 -SNR_peak_min 4.5 -DM_cand 1.5 -N_min 3
    single_pulse_searcher.py -fits ${fits} -no_store -N_min 1 -plot_name ${name}_sps.pdf *.SpS
    """
}

include { ddplan; prepfold } from './pulsar_search_module'
include { classifier }   from './classifier_module'

workflow {
    ddplan( fits_files.map{ it -> [ params.cand + '_' + it.getBaseName().split("/")[-1].split("_ch")[0], it ] } )
    search_dd_fft_acc( // combine the fits files and ddplan with the matching name key (candidateName_obsid_pointing)
                        ddplan.out.splitCsv().map{ it -> [ it[0], [ it[1], it[2], it[3], it[4], it[5], it[6], it[7] ] ] }.\
                        concat(fits_files.map{ it -> [ params.cand + '_' + it.getBaseName().split("/")[-1].split("_ch")[0], it ] }).groupTuple( size: 2 ).\
                        // Find for each ddplan match that with the fits files and the name key then change the format to [val(name), val(dm_values), file(fits_files)]
                        map{ it -> [it[1].init(), [[it[0], it[1].last()]]].combinations() }.flatMap().\
                        map{ it -> [it[1][0], it[0], it[1][1]]} )
    // Get all the inf, ACCEL and single pulse files and sort them into groups with the same name key
    accelsift( search_dd_fft_acc.out.map{ it -> [it[0], [it[1]].flatten().findAll { it != null } + \
                                                        [it[3]].flatten().findAll { it != null }] }.\
                groupTuple( size: 1, remainder: true ).map{ it -> [it[0], it[1].flatten()]} )
    single_pulse_searcher( search_dd_fft_acc.out.map{ it -> [it[0], [it[4]].flatten().findAll { it != null }] }.\
                            groupTuple( size: 1, remainder: true ).map{ it -> [it[0], it[1].flatten()]}.\
                            // Add fits files
                            concat(fits_files.map{ it -> [ params.cand + '_' + it.getBaseName().split("/")[-1].split("_ch")[0], it ] }).groupTuple( size: 2 ).map{ it -> [it[0], it[1][0], it[1][1]]} )
    prepfold( fits_files.map{ it -> [ params.cand + '_' + it.getBaseName().split("/")[-1].split("_ch")[0], it ] }.cross(
                // Group all the accelsift lines together
                accelsift.out.map{ it -> it[1] }.splitCsv().flatten().map{ it -> [it.split()[0].split("_ACCEL")[0], it ] }.cross(
                // Group all the .cand and .inf files by their base names
                search_dd_fft_acc.out.map{ it -> [it[3]].flatten().findAll { it != null } }.
                flatten().map{ it -> [it.baseName.split(".inf")[0], it ] }.concat(
                search_dd_fft_acc.out.map{ it -> [it[5]].flatten().findAll { it != null } }.
                flatten().map{ it -> [it.baseName.split("_ACCEL")[0], it ] }).groupTuple( size: 2 )
                // match the cand and inf file with each accelsift line and reoraganise
                ).map{ it -> [it[0][0].split("_DM")[0], [it[0][1], it[1][1][0], it[1][1][1]]] }
                // Match with fits files and eogranise to val(cand_line), file(cand_file), file(cand_inf), file(fits_files)
                ).map{ it -> [it[1][1][0], it[1][1][2], it[1][1][1], it[0][1]] } )
}