#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.obsid = null
params.fits_file = null
params.out_dir = "${params.search_dir}/${params.obsid}_candidates"

params.dm_min = 1
params.dm_max = 250
params.dm_min_step = 0.02
params.dm_max_step = 0.5
params.max_dms_per_job = 5000

params.cand = "Blind"
params.sp = false

//Defaults for the accelsearch command
params.nharm = 16 // number of harmonics to search
params.min_period = 0.001 // min period to search for in sec (ANTF min = 0.0013)
params.max_period = 30 // max period to search for in sec  (ANTF max = 23.5)
params.zmax = 0 // don't do an accelsearch by default

params.help = false
if ( params.help ) {
    help = """pulsar_search.nf: A pipeline perform a pulsar search on a single input fits file.
             |                  The fits files must be in the format
             |                  <obsid>_<pointing>_ch<min_chan>-<max_chan>_00??.fits
             |Required argurments:
             |  --obsid     Observation ID you want to process   [no default]
             |  --fits_file The fits file to search              [no default]
             |  --dur       Duration of the fits file in seconds [no default]
             |
             | Dedispersion arguments (optional):
             |  --dm_min    Minimum DM to search over [default: ${params.dm_min}]
             |  --dm_max    Maximum DM to search over [default: ${params.dm_max}]
             |  --dm_min_step
             |              Minimum DM step size (Delta DM) [default: ${params.dm_min_step }]
             |  --dm_max_step
             |              Maximum DM step size (Delta DM) [default: ${params.dm_max_step }]
             |  --max_dms_per_job
             |              Maximum number of DM steps a single job will procces.
             |              Lowering this will reduce memory usage and increase parellelisation.
             |              [default: ${params.max_dms_per_job}]
             |
             |Optional arguments:
             |  --cand      Candidate name to do a targeted search [default: Blind]
             |  --sp        Perform a single pulse search [default false]
             |  --zmax      The acceleration range to search over. If you would like to perform
             |              an acceleration search I recomend you use 200 and set
             |              --max_dms_per_job 32
             |              [default: 0 (no acceleration search)]
             |  --out_dir   Output directory for the candidates files
             |              [default: ${params.search_dir}/<obsid>_candidates]
             |  --mwa_search_version
             |              The mwa_search module bersion to use [default: master]
             |  -w          The Nextflow work directory. Delete the directory once the processs
             |              is finished [default: ${workDir}]""".stripMargin()
    println(help)
    exit(0)
}

// Option parsing
if ( params.obsid == null ) {
    println("No obsid, please use --obsid. Exiting.")
    exit(0)
}
if ( params.fits_file ) {
    fits_file = Channel.fromPath( "${params.fits_file}", checkIfExists: true )
    //nfiles = new File("${params.fits_file}").listFiles().findAll { it.name ==~ /.*fits/ }.size()
    fits_file.view( it -> "Running search on ${it}" )
}
else {
    println("No fits file given, please use --fits_file. Exiting.")
    exit(0)
}
if ( params.dur ) {
    params.begin = 1
    params.end = params.dur
}
else {
    println("No duration given, please use --dur. Exiting.")
    exit(0)
}

// If doing an acceleration search, lower the number of DMs per job so the jobs don't time out
if ( params.zmax == 0 ) {
    total_dm_jobs = 6
}
else {
    total_dm_jobs = 24
    params.max_dms_per_job = 128
}

include {pulsar_search; single_pulse_search} from './pulsar_search_module'
include { classifier }   from './classifier_module'

workflow {
    if ( params.sp ) {
        //single_pulse_search( fits_file.toSortedList().map{ it -> [ params.cand + '_' + it[0].getBaseName().split("/")[-1].split("_ch")[0], it ] }.view() )
        single_pulse_search( fits_file.map{ it -> [ params.cand + '_' + it.getBaseName().split("/")[-1].split("_ch")[0], it ] } )
    }
    else {
        pulsar_search( fits_file.toSortedList().map{ it -> [ params.cand + '_' + it[0].getBaseName().split("/")[-1].split("_ch")[0], it ] } )
        classifier( pulsar_search.out[1].flatten().collate( 120 ) )
    }
}
