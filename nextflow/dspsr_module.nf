nextflow.enable.dsl = 2


process pdmp {
    label 'cpu'
    time '6h'

    input:
    path fits
    val pointings

    output:
    path "*pdmp*"

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