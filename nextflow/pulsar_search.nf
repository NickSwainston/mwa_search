#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

params.obsid = null
params.pointings = null
params.fitsdir = "/group/mwavcs/vcs/${params.obsid}/pointings"
params.out_dir = "${params.search_dir}/${params.obsid}_candidates"
params.dm_min = 1
params.dm_max = 250
params.dm_min_step = 0.02

params.scratch = false
params.fits_file_dir = false

params.cand = "Blind"
params.sp = false

//Defaults for the accelsearch command
params.nharm = 16 // number of harmonics to search
params.min_period = 0.001 // min period to search for in sec (ANTF min = 0.0013)
params.max_period = 30 // max period to search for in sec  (ANTF max = 23.5)

pointing = Channel.from( params.pointings )
if ( params.fits_file_dir ) {
    fits_files = Channel.fromPath( "${params.fits_file_dir}/${params.obsid}_*.fits", checkIfExists: true )
        nfiles = new File("${params.fits_file_dir}").listFiles().findAll { it.name ==~ /.*1*fits/ }.size()
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

// Work out length of obs, may over estimate up to 200 seconds
params.end = obs_length = nfiles * 200
params.begin = 1

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

include {pulsar_search; single_pulse_search} from './pulsar_search_module'
include { classifier }   from './classifier_module'

workflow {
    if ( params.sp ) {
        single_pulse_search( fits_files.toSortedList().map{ it -> [ params.cand + '_' + it[0].getBaseName().split("/")[-1].split("_ch")[0], it ] } )
    }
    else {
        pulsar_search( fits_files.toSortedList().map{ it -> [ params.cand + '_' + it[0].getBaseName().split("/")[-1].split("_ch")[0], it ] } )
        classifier( pulsar_search.out[1].flatten().collate( 120 ) )
    }
}