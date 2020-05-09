#!/usr/bin/env nextflow

nextflow.preview.dsl = 2
include { pre_beamform; beamform; beamform_ipfb } from './beamform_module'

params.obsid = null
params.calid = null
params.pointings = null
params.pointing_file = null

params.begin = null
params.end = null
params.all = false

params.summed = false
params.ipfb = false
params.vcstools_version = 'master'
params.mwa_search_version = 'master'

params.didir = "${params.basedir}/${params.obsid}/cal/${params.calid}/rts"
params.publish_fits = true
params.publish_fits_scratch = false

params.no_combined_check = false

params.help = false
if ( params.help ) {
    help = """beamform.nf: A pipeline that will beamform and splice on all input pointings.
             |Required argurments:
             |  --obsid     Observation ID you want to process [no default]
             |  --calid     Observation ID of calibrator you want to process [no default]
             |  --pointings A space sepertated list of pointings with the RA and Dec seperated
             |              by _ in the format HH:MM:SS_+DD:MM:SS, e.g.
             |              "19:23:48.53_-20:31:52.95 19:23:40.00_-20:31:50.00" [default: None]
             |  --pointing_file
             |              A file containing pointings with the RA and Dec seperated by _
             |              in the format HH:MM:SS_+DD:MM:SS on each line, e.g.
             |              "19:23:48.53_-20:31:52.95\\n19:23:40.00_-20:31:50.00" [default: None]
             |  --begin     First GPS time to process [no default]
             |  --end       Last GPS time to process [no default]
             |  --all       Use entire observation span. Use instead of -b & -e. [default: false]
             |  --publish_fits
             |              Publish to the fits directory (/group on Galaxy). Include this
             |              option.
             |
             |Optional arguments:
             |  --summed    Add this flag if you the beamformer output as summed polarisations
             |              (only Stokes I). This reduces the data size by a factor of 4.
             |              [default: False]
             |  --ipfb      Perform an the inverse PFB to produce high time resolution beamformed
             |              vdif files [default: false]
             |  --publish_fits_scratch
             |              Publish to the scratch fits directory (/astro on Galaxy). Use this
             |              instead of --publish_fits_scratch
             |  --vcstools_version
             |              The vcstools module version to use [default: master]
             |  --mwa_search_version
             |              The mwa_search module bersion to use [default: master]
             |  --no_combined_check
             |              Don't check if all the combined files are available [default: false]
             |  -w          The Nextflow work directory. Delete the directory once the processs
             |              is finished [default: ${workDir}]""".stripMargin()
    println(help)
    exit(0)
}

if ( params.pointing_file ) {
    pointings = Channel
        .fromPath(params.pointing_file)
        .splitCsv()
        .collect()
        .flatten()
        .collate( params.max_pointings )
}
else if ( params.pointings ) {
    pointings = Channel
        .from(params.pointings.split(","))
        .collect()
        .flatten()
        .collate( params.max_pointings )
}
else {
    println "No pointings given. Either use --pointing_file or --pointings. Exiting"
    exit(1)
}

workflow {
    pre_beamform()
    if ( params.ipfb ) {
        beamform_ipfb( pre_beamform.out[0],\
                       pre_beamform.out[1],\
                       pre_beamform.out[2],\
                       pointings )
    }
    else {
        beamform( pre_beamform.out[0],\
                  pre_beamform.out[1],\
                  pre_beamform.out[2],\
                  pointings )
    }
}