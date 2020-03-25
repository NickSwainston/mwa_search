nextflow.preview.dsl = 2
include pulsar_search from './pulsar_search_module'

params.obsid = 1253471952
params.pointing = '02:55:56.29_-53:04:21.27'
params.fitsdir = "/group/mwaops/vcs/${params.obsid}/pointings"
params.dm_min = 1
params.dm_max = 250

//Defaults for the accelsearch command
params.nharm = 16 // number of harmonics to search
params.min_period = 0.001 // min period to search for in sec (ANTF min = 0.0013)
params.max_period = 30 // max period to search for in sec  (ANTF max = 23.5)

pointing = Channel.from( params.pointing )
fits_files = Channel.fromPath( "${params.fitsdir}/${params.pointing}/${params.obsid}_*.fits" ).collect()



workflow {
    pulsar_search( fits_files,\
                   pointing )
    publish:
        pulsar_search.out to: "/group/mwaops/vcs/${params.obsid}/pointings/${pointing}" //Change maybe
}

