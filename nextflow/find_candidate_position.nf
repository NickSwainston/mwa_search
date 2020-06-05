#!/usr/bin/env nextflow

nextflow.preview.dsl = 2
include { pre_beamform; beamform } from './beamform_module'

params.obsid = null
params.calid = null
params.pointings = null
params.pointing_file = null

params.begin = null
params.end = null
params.all = false

params.summed = true
params.vcstools_version = 'master'
params.mwa_search_version = 'master'

params.basedir = '/group/mwavcs/vcs'
params.didir = "${params.basedir}/${params.obsid}/cal/${params.calid}/rts"
params.channels = null
params.out_dir = "${params.basedir}/${params.obsid}/pointings"

params.bins = 128
params.period = 0.90004
params.dm = 23.123
params.subint = 60
params.nchan = 48

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

process pdmp {
    label 'cpu'
    time '4h'

    input:
    file fits
    val pointings

    output:
    file "*pdmp*"

    beforeScript "module use ${params.presto_module_dir}; module load dspsr/master"

    """
    dspsr -b ${params.bins} -c ${params.period} -D ${params.dm} -L ${params.subint} -e subint -cont -U 4000 ${fits}
    psradd *.subint -o ${params.obsid}_${pointings}.ar
    pam --setnchn ${params.nchan} -m ${params.obsid}_${pointings}.ar
    pdmp -g ${params.obsid}_${pointings}_pdmp.ps/cps ${obsid}_${pointings}.ar
    """

}


workflow {
    pre_beamform()
    beamform( pre_beamform.out[0],\
              pre_beamform.out[1],\
              pre_beamform.out[2],\
              pointings )
    pdmp( beamform.out[1], beamform.out[2] )
    publish:
        pdmp.out to: params.out_dir
}