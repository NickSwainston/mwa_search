nextflow.enable.dsl = 2

params.obsid = null
params.pointings = null

params.out_dir = "${params.scratch_basedir}/${params.obsid}/pointings"

params.bins = 128
params.period = 0.90004
params.dm = 23.123
params.subint = 60
params.nchan = 48


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
             |Optional arguments:
             |  --eph       The ephermis file to fold on the obs. If you don't have one use --period and --dm
             |  --period    The topo period of the pulsar in seconds. 
             |  --dm        The dispersion measure of the pulsar.
             |  --out_dir   Output directory for the candidates files
             |              [default: ${params.search_dir}/<obsid>_candidates]
             |  -w          The Nextflow work directory. Delete the directory once the processs
             |              is finished [default: ${workDir}]""".stripMargin()
    println(help)
    exit(0)
}

if ( params.fits_file ) {
    fits_file = Channel.fromPath( "${params.fits_file}", checkIfExists: true )
    //nfiles = new File("${params.fits_file}").listFiles().findAll { it.name ==~ /.*fits/ }.size()
    fits_file.view( it -> "Running pipeline on ${it}" )
}

process fits_to_archive_to_fits {
    label 'cpu'
    time '6h'

    input:
    file fits

    output:
    file "*pdmp*"

    container = "file:///${params.containerDir}/dspsr/dspsr.sif"

    """
    dspsr -t $task.cpus -b ${params.bins} -c ${params.period} -D ${params.dm} -L 5000 -A -cont -U 4000 ${fits}
    pam -a PSRFITS -e fits *ar
    """

}

workflow polarimetry {
    take:
        fits
    main:
        fits_to_archive_to_fits( fits )
    emit:
        fits_to_archive_to_fits.out
}

workflow {
    polarimetry(fits_file)
}