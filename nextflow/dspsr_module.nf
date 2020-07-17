nextflow.preview.dsl = 2

params.obsid = null
params.pointings = null

params.out_dir = "${params.basedir}/${params.obsid}/pointings"

params.bins = 128
params.period = 0.90004
params.dm = 23.123
params.subint = 60
params.nchan = 48

process pdmp {
    label 'cpu'
    time '6h'

    input:
    file fits
    val pointings

    output:
    file "*pdmp*"

    beforeScript "module use ${params.presto_module_dir}; module load dspsr/master"

    """
    dspsr  -t $task.cpus -b ${params.bins} -c ${params.period} -D ${params.dm} -L ${params.subint} -e subint -cont -U 4000 ${fits}
    psradd *.subint -o ${params.obsid}_${pointings}.ar
    pam --setnchn ${params.nchan} -m ${params.obsid}_${pointings}.ar
    pdmp -g ${params.obsid}_${pointings}_pdmp.ps/cps ${params.obsid}_${pointings}.ar
    """

}

workflow pdmp_wf {
    take:
        fits
        pointings
    main:
        pdmp( fits,\
              pointings )
    emit:
        pdmp.out
}