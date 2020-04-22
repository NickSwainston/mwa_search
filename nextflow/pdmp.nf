nextflow.preview.dsl = 2
include pdmp_wf from './dspsr_module'

params.obsid = null
params.pointings = null

params.out_dir = "${params.basedir}/${params.obsid}/pointings"

params.bins = 128
params.period = 0.90004
params.dm = 23.123
params.subint = 60
params.nchan = 48

if ( params.pointings ) {
    pointings = Channel.from(params.pointings)
}

fits = Channel.fromPath( "${params.basedir}/${params.obsid}/pointings/${params.pointings}/${params.obsid}*fits" ).collect()

workflow {
    pdmp_wf( fits,
             pointings )
    publish:
        pdmp_wf.out to: params.out_dir
}