#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

params.obsid = null
params.calid = null

params.begin = null
params.end = null
params.all = false

params.vcstools_version = 'master'
params.mwa_search_version = 'master'

params.didir = "${params.scratch_basedir}/${params.obsid}/cal/${params.calid}/rts"
params.publish_fits = false
params.publish_fits_scratch = true

params.out_dir = "${params.search_dir}/${params.obsid}_candidates"

params.no_combined_check = false

params.help = false
if ( params.help ) {
    help = """beamform_fov_sources.nf: A pipeline that will beamform on all pulsars in the FOV
             |                        and perform a search on all pulsar candidates.
             |Required argurments:
             |  --obsid     Observation ID you want to process [no default]
             |  --calid     Observation ID of calibrator you want to process [no default]
             |  --begin     First GPS time to process [no default]
             |  --end       Last GPS time to process [no default]
             |  --all       Use entire observation span. Use instead of -b & -e. [default: false]
             |  --publish_fits_scratch
             |              Publish to the scratch fits directory (/astro on Galaxy). Include
             |              this option.
             |
             |Optional arguments:
             |  --publish_fits
             |              Publish to the fits directory (/group on Galaxy). Use this instead
             |              of --publish_fits_scratch
             |  --vcstools_version
             |              The vcstools module version to use [default: master]
             |  --mwa_search_version
             |              The mwa_search module bersion to use [default: master]
             |  --no_combined_check
             |              Don't check if all the combined files are available [default: false]
             |  --out_dir   Where the search candidates will be output
             |              [default: /group/mwavcs/vcs/<obsid>/<obsid>_candidates]
             |  -w          The Nextflow work directory. Delete the directory once the processs
             |              is finished [default: ${workDir}]""".stripMargin()
    println(help)
    exit(0)
}

// Work out some estimated job times
if ( "$HOSTNAME".startsWith("farnarkle") ) {
    // In seconds
    search_dd_fft_acc_dur = obs_length * 5.0
    prepfold_dur = obs_length * 2.0
    presto_python_load = "module use ${params.presto_module_dir}; module load presto/${params.presto_module}; module load python/2.7.14; module load matplotlib/2.2.2-python-2.7.14"
}
else {
    search_dd_fft_acc_dur = 14400
    prepfold_dur = 7200
    presto_python_load = ""
}

process find_pointings {
    input:
    tuple val(begin), val(end)

    output:
    file "${params.obsid}_fov_sources.csv"

    """
    pulsars_in_fov.py -o $params.obsid -b $begin -e $end
    """
}

process known_pulsar_prepfold {
    label 'cpu'
    time "${prepfold_dur}s"
    errorStrategy 'retry'
    maxRetries 1

    input:
    tuple file(fits_files), val(pulsar)

    output:
    file "*pfd*"

    if ( "$HOSTNAME".startsWith("farnarkle") ) {
        beforeScript "module use ${params.presto_module_dir}; module load presto/${params.presto_module}"
    }
    else if ( "$HOSTNAME".startsWith("x86") || "$HOSTNAME".startsWith("garrawarla") || "$HOSTNAME".startsWith("galaxy") ) {
        container = "file:///${config.containerDir}/presto/presto.sif"
    }
    else {
        container = "nickswainston/presto:realfft_docker"
    }
    //no mask command currently
    """
    echo "pulsar: ${pulsar}"
    # Set up the prepfold options to match the ML candidate profiler
    period=\$(psrcat -c p0 -x ${pulsar} | cut -d ' ' -f 1)
    psrcat -e ${pulsar} > ${pulsar}.eph
    if (( \$(echo "\$period > 0.01" | bc -l) )); then
        nbins=100
        ntimechunk=120
        dmstep=1
        period_search_n=1
    else
        # bin size is smaller than time resolution so reduce nbins
        nbins=50
        ntimechunk=40
        dmstep=3
        period_search_n=2
    fi

    prepfold -ncpus $task.cpus -o initial_fold -n \$nbins -noxwin -noclip -par ${pulsar}.eph -nsub 256 \
-npart \$ntimechunk -dmstep \$dmstep -pstep 1 -pdstep 2 -npfact \$period_search_n -ndmfact 1 -runavg *.fits
    """
}

include { pre_beamform; beamform; beamform_ipfb; get_beg_end } from './beamform_module'
include { pulsar_search; single_pulse_search } from './pulsar_search_module'
include { classifier } from './classifier_module'

workflow initial_fold {
    take:
        pointings
        pulsar_names
        fits_files
    main:
        known_pulsar_prepfold( // Combine pointings with pulsar names
                               pointings.merge( pulsar_names ).\
                               // Then group them with available fits files
                               concat( fits_files ).groupTuple( size: 2, remainder: false ).map{ it -> [it[1][1], it[1][0]] } )
        // Run through the classfier
        classifier( known_pulsar_prepfold.out.flatten().collate( 120 ) )
    emit:
        classifier.out[0]
}


workflow {
    pre_beamform()
    find_pointings( pre_beamform.out[0] )
    beamform( pre_beamform.out[0],\
              pre_beamform.out[1],\
              pre_beamform.out[2],\
              //Grab the pointings for slow pulsars and single pulses
              find_pointings.out.splitCsv(skip: 1, limit: 1).mix(\
              find_pointings.out.splitCsv(skip: 5, limit: 1),\
              find_pointings.out.splitCsv(skip: 7, limit: 1)).collect().flatten().unique().collate( params.max_pointings ) )
    beamform_ipfb( pre_beamform.out[0],\
                   pre_beamform.out[1],\
                   pre_beamform.out[2],\
                   //Grab the pointings for slow pulsars and single pulses
                   find_pointings.out.splitCsv(skip: 3, limit: 1) )

    // Perform processing pipeline on all known pulsars
    initial_fold( // pointings
                  find_pointings.out.splitCsv(skip: 1, limit: 1).mix( find_pointings.out.splitCsv(skip: 3, limit: 1) ).collect().flatten(),\
                  // pulsar names
                  find_pointings.out.splitCsv(skip: 0, limit: 1).mix( find_pointings.out.splitCsv(skip: 2, limit: 1) ).collect().flatten(),\
                  // fits files
                  beamform.out[3].mix(beamform_ipfb.out[3]) )

    // Perform a search on all candidates (not known pulsars)
    // if pointing in fits file name is in pulsar search pointing list
    pulsar_search( find_pointings.out.splitCsv(skip: 5, limit: 1).flatten().merge(find_pointings.out.splitCsv(skip: 4, limit: 1).flatten()).\
                   concat(beamform.out[3]).groupTuple( size: 2, remainder: false ).map { it -> [ "Blind_${params.obsid}_${it[0]}".toString(), it[1][1] ] } )
    classifier( pulsar_search.out[1].flatten().collate( 600 ) )
    // Perform a single pulse search on all single pulse candidates
    single_pulse_search( find_pointings.out.splitCsv(skip: 7, limit: 1).flatten().merge(find_pointings.out.splitCsv(skip: 6, limit: 1).flatten()).\
                         concat(beamform.out[3]).groupTuple( size: 2, remainder: false ).map { it -> it[1] } )
}