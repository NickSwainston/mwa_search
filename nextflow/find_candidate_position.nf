#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

params.obsid = null
params.calid = null
params.pointings = null
params.pointing_file = null
params.begin = 0
params.end = 0
params.all = false

params.pointing_grid = null
params.fwhm_deg = null
params.fraction = 0.8
params.loops = 1

params.summed = true
params.channels = null
params.vcstools_version = 'master'
params.mwa_search_version = 'master'

params.bins = 128
params.period = 0.90004
params.dm = 23.123
params.subint = 60
params.nchan = 48

params.no_pdmp = false
params.fwhm_ra = "None"
params.fwhm_dec = "None"



params.help = false
if ( params.help ) {
    help = """mwa_search_pipeline.nf: A pipeline that will beamform and perform a pulsar search
             |                        in the entire FOV.
             |Required argurments:
             |  --obsid     Observation ID you want to process [no default]
             |  --calid     Observation ID of calibrator you want to process [no default]
             |  --begin     First GPS time to process [no default]
             |  --end       Last GPS time to process [no default]
             |  --all       Use entire observation span. Use instead of -b & -e. [default: false]
             |  --summed    Add this flag if you the beamformer output as summed polarisations
             |              (only Stokes I). This reduces the data size by a factor of 4.
             |              [default: False]
             |
             |Required pointing arguments:
             |  --pointings A space sepertated list of pointings with the RA and Dec seperated
             |              by _ in the format HH:MM:SS_+DD:MM:SS, e.g.
             |              "19:23:48.53_-20:31:52.95 19:23:40.00_-20:31:50.00" [default: None]
             |  --pointing_file
             |              A file containing pointings with the RA and Dec seperated by _
             |              in the format HH:MM:SS_+DD:MM:SS on each line, e.g.
             |              "19:23:48.53_-20:31:52.95\\n19:23:40.00_-20:31:50.00" [default: None]
             |
             |Pointing grid arguments:
             | --pointing_grid
             |              Pointing which grid.py will make a loop of pointings around eg.
             |              "19:23:48.53_-20:31:52.95" [default: None]
             | --fwhm_deg   The FWHM of the observation in degrees (used by grid.py) [default: 0.021]
             | --fraction   The fraction of the FWHM to space the grid by [default: 0.8]
             | --loops      The number of loops of beamd to surround the centre pointing [default: 1]
             |
             |Presto and dspsr options:
             | --bins       Number of bins to use [default: 128]
             | --period     Period in seconds to fold on [default: 0.90004]
             | --dm         The dispersion measure to use [default: 23.123]
             | --subint     The number of subints to use in pmdp [default: 60]
             | --nchan      The number of subchans to use in pmdp [default: 48]
             |
             |Optional arguments:
             |  --out_dir   Output directory for the candidates files
             |              [default: ${params.search_dir}/<obsid>_candidates]
             |  --ipfb      Perform an the inverse PFB to produce high time resolution beamformed
             |              vdif files [default: false]
             |  --publish_fits
             |              Publish to the fits directory (/group on Galaxy).
             |  --publish_fits_scratch
             |              Publish to the scratch fits directory (/astro on Galaxy).
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

include { pre_beamform; beamform } from './beamform_module'
include { fwhm_calc } from './data_processing_pipeline'

params.didir = "${params.scratch_basedir}/${params.obsid}/cal/${params.calid}/rts"
params.out_dir = "${params.search_dir}/${params.obsid}_candidate_follow_up"

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
else if ( params.pointing_grid ) {
    pointing_grid = Channel.from(params.pointing_grid)
}
else {
    println "No pointings given. Either use --pointing_file, --pointings or --pointing_grid. Exiting"
    exit(1)
}

if ( params.no_pdmp ) {
    input_sn_option = " -b "
}
else {
    input_sn_option = " -p "
}

process grid {
    input:
    val pointings
    val fwhm

    output:
    file "*txt"

    """
    grid.py -o $params.obsid -d $fwhm -f $params.fraction -p $pointings -l $params.loops
    """

}

process prepfold {
    label 'cpu'
    time '3h'
    publishDir params.out_dir, mode: 'copy'

    input:
    tuple val(pointing), file(fits_files)

    output:
    file "*bestprof"
    file "*png"

    if ( "$HOSTNAME".startsWith("farnarkle") ) {
        beforeScript "module use ${params.presto_module_dir}; module load presto/${params.presto_module}"
    }
    else if ( "$HOSTNAME".startsWith("x86") || "$HOSTNAME".startsWith("garrawarla") || "$HOSTNAME".startsWith("galaxy") ) {
        container = "file:///${params.containerDir}/presto/presto.sif"
    }
    else {
        container = "nickswainston/presto:realfft_docker"
    }

    //no mask command currently
    """
    prepfold -ncpus $task.cpus -o ${params.obsid}_${pointing} -n ${params.bins} -noxwin -noclip -p ${params.period} -dm ${params.dm} -nsub 256 -npart 120 \
-dmstep 1 -pstep 1 -pdstep 2 -npfact 1 -ndmfact 1 -runavg ${params.obsid}*.fits
    """
}

process pdmp {
    label 'cpu'
    time '6h'
    publishDir params.out_dir, mode: 'copy'

    when:
    params.no_pdmp == false

    input:
    file bestprof
    file fits_files
    val pointings

    output:
    file "*ps"
    file "*posn"

    if ( "$HOSTNAME".startsWith("farnarkle") ) {
        beforeScript "module use ${params.presto_module_dir}; module load dspsr/master"
    }
    else if ( "$HOSTNAME".startsWith("x86") || "$HOSTNAME".startsWith("garrawarla") || "$HOSTNAME".startsWith("galaxy") ) {
        container = "file:///${params.containerDir}/dspsr/dspsr.sif"
    }
    else {
        container = "nickswainston/dspsr_docker"
    }

    //may need to add some channel names
    """
    DM=\$(grep DM *.bestprof | tr -s ' ' | cut -d ' ' -f 5)
    echo "DM: \$DM"
    period=\$(grep P_topo *.bestprof | tr -s ' ' | cut -d ' ' -f 5)
    period="\$(echo "scale=10;\${period}/1000"  |bc)"
    echo "period: \$period"
    samples="\$(grep "Data Folded" *.bestprof | tr -s ' ' | cut -d ' ' -f 5)"
    #One subint per 30 seconds
    subint=\$(python -c "print('{:d}'.format(int(\$samples/300000)))")
    if [ \$subint -lt 30 ]; then subint=30; fi
    echo "subint: \$subint"
    dspsr -t $task.cpus -b ${params.bins} -c \${period} -D \${DM} -L \${subint} -e subint -cont -U 4000 ${params.obsid}*.fits
    psradd *.subint -o ${params.obsid}_${pointings}.ar
    pam --setnchn ${params.nchan} -m ${params.obsid}_${pointings}.ar
    pdmp -g ${params.obsid}_${pointings}_pdmp.ps/cps ${params.obsid}_${pointings}.ar
    mv pdmp.posn ${params.obsid}_${pointings}_pdmp.posn
    """
}

process bestgridpos {
    publishDir params.out_dir, mode: 'copy'

    input:
    file posn_or_bestprof
    val fwhm

    output:
    file "*txt"
    file "*png"

    """
    if [[ ${params.fwhm_ra} == None || ${params.fwhm_dec} == None ]]; then
        bestgridpos.py -o ${params.obsid} ${input_sn_option} ./ -w -fr ${fwhm} -fd ${fwhm}
    else
        bestgridpos.py -o ${params.obsid} ${input_sn_option} ./ -w -fr ${params.fwhm_ra} -fd ${params.fwhm_dec}
    fi
    """
}

workflow find_pos {
    take:
        pointings
        pre_beamform_1
        pre_beamform_2
        pre_beamform_3
        fwhm
    main:
        beamform( pre_beamform_1,\
                  pre_beamform_2,\
                  pre_beamform_3,\
                  pointings )
        prepfold( beamform.out[3] )
        if ( params.no_pdmp ) {
            bestgridpos( prepfold.out[0].collect(),\
                         fwhm )
        }
        else {
            pdmp( prepfold.out[0],
                  beamform.out[1],
                  beamform.out[2] )
            bestgridpos( pdmp.out[1].collect(),\
                         fwhm )
        }
    emit:
        bestgridpos.out[0].splitCsv().collect().flatten().collate( params.max_pointings )
}

workflow {
    pre_beamform()
    fwhm_calc( pre_beamform.out[1] )
    if ( params.pointing_grid ) {
        grid( pointing_grid,\
              fwhm_calc.out.splitCsv().flatten() )
        find_pos( grid.out.splitCsv().collect().flatten().collate( params.max_pointings ),\
                  pre_beamform.out[0],\
                  pre_beamform.out[1],\
                  pre_beamform.out[2],\
                  fwhm_calc.out.splitCsv().flatten() )
    }
    else if ( params.pointing_file || params.pointings ) {
        find_pos( pointings,\
                  pre_beamform.out[0],\
                  pre_beamform.out[1],\
                  pre_beamform.out[2],\
                  fwhm_calc.out.splitCsv().flatten() )
    }
    beamform( pre_beamform.out[0],\
                pre_beamform.out[1],\
                pre_beamform.out[2],\
                find_pos.out.view() )
    prepfold( beamform.out[3] )
    pdmp( prepfold.out[0],
          beamform.out[1],
          beamform.out[2] )
}
