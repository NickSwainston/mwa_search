#!/usr/bin/env nextflow

include pdmp_wf from './dspsr_module'


params.out_dir = "${params.vcsdir}/${params.obsid}/pointings"


if ( params.pointings ) {
    pointings = Channel.from(params.pointings)
}

fits = Channel.fromPath( "${params.vcsdir}/${params.obsid}/pointings/${params.pointings}/${params.obsid}*fits" ).collect()

workflow {
    pdmp_wf( fits,
             pointings )
    publish:
        pdmp_wf.out to: params.out_dir
}