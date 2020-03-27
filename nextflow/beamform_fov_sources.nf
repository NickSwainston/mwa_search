nextflow.preview.dsl = 2
include { pre_beamform; beamform; beamform_ipfb; get_beg_end } from './beamform_module'

params.obsid = null
params.pointings = null
params.calid = null

params.begin = null
params.end = null
params.all = false

params.summed = false
params.vcstools_version = 'master'

params.basedir = '/group/mwaops/vcs'
params.scratch_basedir = '/astro/mwaops/vcs'
params.didir = "${params.basedir}/${params.obsid}/cal/${params.calid}/rts"
params.channels = null
params.publish_fits = true

params.no_combined_check = false



process find_pointings {
    input:
    tuple val(begin), val(end)

    output:
    file "${params.obsid}_fov_sources.csv"

    """
    pulsars_in_fov.py -o $params.obsid -b $begin -e $end
    """
}

workflow {
    get_beg_end()
    find_pointings( get_beg_end.out.map{ it.split(",") }.flatten().collect() )
    pre_beamform()
    beamform( pre_beamform.out[0],\
              pre_beamform.out[1],\
              pre_beamform.out[2],\
              //Grab the pointigngs for slow pulsars and single pulses
              find_pointings.out.splitCsv(skip: 1, limit: 1).mix(\
              find_pointings.out.splitCsv(skip: 5, limit: 1)).collect().unique().collate( 15 ) )
    beamform_ipfb( pre_beamform.out[0],\
                   pre_beamform.out[1],\
                   pre_beamform.out[2],\
                   //Grab the pointigngs for slow pulsars and single pulses
                   find_pointings.out.splitCsv(skip: 3, limit: 1) )
}