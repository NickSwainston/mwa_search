nextflow.preview.dsl = 2

params.vcstools_version = 'master'
params.scratch_basedir = '/astro/mwavcs/vcs'
params.version_compare = "V2.3_FEE2016"
params.version_compare_base_dir = "${params.scratch_basedir}/beamformer_tests"
params.version_compare_dir = "${params.version_compare_base_dir}/${params.version_compare}"
params.publish_version = false


// All tests are 300 seconds long
obs_length = 300
max_job_pointings = 5


//Beamforming ipfb duration calc
mb_ipfb_dur = ( obs_length * (params.bm_read + 3 * (params.bm_cal + params.bm_beam) + params.bm_write) + 200 ) * 1.2

//Beamforming duration calc
mb_dur = ( obs_length * (params.bm_read + params.bm_cal + max_job_pointings * (params.bm_beam +params.bm_write)) + 200 ) * 1.2

//Required temp SSD mem required for gpu jobs
temp_mem = (int) (0.0048 * obs_length * max_job_pointings + 1)
temp_mem_single = (int) (0.0096 * obs_length + 2)

bf_out = " -p "



params.help = false
if ( params.help ) {
    help = """test_beamformer.nf: A pipeline that will test the beamformer by comparing the binaires of the output psrfits 
             |                    and vdif files with previous runs and the signal-to-noise of pulsar detections
             |
             |Optional arguments:
             |  --vcstools_version
             |              The vcstools module version to use [default: master]
             |  --version_compare
             |              The beamformer version to compare the current beamformer with.
             |              The previous runs are contained in ${params.version_compare_base_dir}
             |              [default: ${params.version_compare}]
             |  --publish_version
             |              Publish the current beamformer run to ${params.version_compare_base_dir} so it can be compared in future runs.
             |              Make sure the --version_compare is different than any other previous runs.
             |              [default: false]
             |  -w          The Nextflow work directory. Delete the directory once the processs
             |              is finished [default: ${workDir}]""".stripMargin()
    println(help)
    exit(0)
}


process make_beam {
    label 'gpu'
    //time '2h'
    time "${mb_dur*task.attempt}s"
    errorStrategy 'retry'
    maxRetries 1
    if ( "$HOSTNAME".startsWith("garrawarla") ) {
        maxForks 70
    }
    else {
        maxForks 120
    }

    if ( "$HOSTNAME".startsWith("farnarkle") ) {
        clusterOptions = "--gres=gpu:1  --tmp=${temp_mem}GB"
        scratch '$JOBFS'
        beforeScript "module use ${params.module_dir}; module load vcstools/${params.vcstools_version}"
    }
    else if ( "$HOSTNAME".startsWith("x86") ) {
        clusterOptions = "--gres=gpu:1"
        scratch '/ssd'
        //container = "file:///${params.containerDir}/vcstools/vcstools_${params.vcstools_version}.sif"
        beforeScript "module use ${params.module_dir}; module load vcstools/${params.vcstools_version}"
    }
    else if ( "$HOSTNAME".startsWith("garrawarla") ) {
        clusterOptions = "--gres=gpu:1  --tmp=${temp_mem}GB"
        scratch '/nvmetmp'
        //container = "file:///${params.containerDir}/vcstools/vcstools_${params.vcstools_version}.sif"
        beforeScript "module use ${params.module_dir}; module load vcstools/${params.vcstools_version}"
    }
    else if ( "$HOSTNAME".startsWith("galaxy") ) {
        beforeScript "module use ${params.module_dir}; module load vcstools/${params.vcstools_version}"
    }
    else {
        container = "cirapulsarsandtransients/vcstools:${params.vcstools_version}"
    }

    input:
    val obsid
    val calid
    val label
    val channel_pair
    val utc
    val point
    tuple val(begin), val(end)

    output:
    file "*fits"

    """
    make_beam -o $obsid -b $begin -e $end -a 128 -n 128 \
-f ${channel_pair[0]} -J ${params.scratch_basedir}/${obsid}/cal/${calid}/rts/DI_JonesMatrices_node${channel_pair[1]}.dat \
-d ${params.scratch_basedir}/${obsid}_beamformer_test_data/combined -P ${point.join(",")} \
-r 10000 -m ${params.scratch_basedir}/${obsid}/${obsid}_metafits_ppds.fits \
${bf_out} -t 6000 -z $utc
    
    # Label all outputs
    for i in \$(ls */*); do mv \${i} ${label}_\${i##*/}; done
    """
}


process make_beam_ipfb {
    publishDir "${params.basedir}/${obsid}/pointings/${point}", mode: 'copy', enabled: params.publish_version, pattern: "*hdr"
    publishDir "${params.basedir}/${obsid}/pointings/${point}", mode: 'copy', enabled: params.publish_version, pattern: "*vdif"

    label 'gpu'
    //time '2h'
    time "${mb_ipfb_dur*task.attempt}s"
    errorStrategy 'retry'
    maxRetries 1
    if ( "$HOSTNAME".startsWith("garrawarla") ) {
        maxForks 70
    }
    else {
        maxForks 120
    }
    
    if ( "$HOSTNAME".startsWith("farnarkle") ) {
        clusterOptions = "--gres=gpu:1  --tmp=${temp_mem_single}GB"
        scratch '$JOBFS'
        beforeScript "module use ${params.module_dir}; module load vcstools/${params.vcstools_version}"
    }
    else if ( "$HOSTNAME".startsWith("x86") ) {
        clusterOptions = "--gres=gpu:1"
        scratch '/ssd'
        container = "file:///${params.containerDir}/vcstools/vcstools_${params.vcstools_version}.sif"
    }
    else if ( "$HOSTNAME".startsWith("garrawarla") ) {
        clusterOptions = "--gres=gpu:1  --tmp=${temp_mem_single}GB"
        scratch '/nvmetmp'
        //container = "file:///${params.containerDir}/vcstools/vcstools_${params.vcstools_version}.sif"
        beforeScript "module use ${params.module_dir}; module load vcstools/${params.vcstools_version}"
    }
    else if ( "$HOSTNAME".startsWith("galaxy") ) {
        beforeScript "module use ${params.module_dir}; module load vcstools/${params.vcstools_version}"
    }
    else {
        container = "cirapulsarsandtransients/vcstools:${params.vcstools_version}"
    }

    input:
    val obsid
    val calid
    val label
    val channel_pair
    val utc
    each point
    tuple val(begin), val(end)

    output:
    file "*fits"
    file "*vdif"

    """
    make_beam -o $obsid -b $begin -e $end -a 128 -n 128 \
-f ${channel_pair[0]} -J ${params.scratch_basedir}/${obsid}/cal/${calid}/rts/DI_JonesMatrices_node${channel_pair[1]}.dat \
-d ${params.scratch_basedir}/${obsid}_beamformer_test_data/combined -P ${point} \
-r 10000 -m ${params.scratch_basedir}/${obsid}/${obsid}_metafits_ppds.fits \
-p -v -t 6000 -z $utc
    # Label all outputs
    for i in \$(ls */*); do mv \${i} ${label}_\${i##*/}; done
    for i in \$(ls *vdif); do mv \${i} ${label}_\${i##*/}; done
    """
}

process compare_repeats {
    echo true
    validExitStatus 0,2

    input:
    file fits_1
    file vdif_1
    file fits_2
    file vdif_2
    file fits_3
    file vdif_3

    """
    # Remove headers
    tail -n +2 ${fits_1} > no_header_1.fits
    tail -n +2 ${fits_2} > no_header_2.fits
    tail -n +2 ${fits_3} > no_header_3.fits
    tail -n +2 ${vdif_1} > no_header_1.vdif
    tail -n +2 ${vdif_2} > no_header_2.vdif
    tail -n +2 ${vdif_3} > no_header_3.vdif

    # Fits check
    fits_diff_success=True
    diff no_header_1.fits no_header_2.fits
    if [ \$? != 0 ]; then fits_diff_success=False; fi
    diff no_header_3.fits no_header_2.fits
    if [ \$? != 0 ]; then fits_diff_success=False; fi
    diff no_header_1.fits no_header_3.fits
    if [ \$? != 0 ]; then fits_diff_success=False; fi
    if [ \$fits_diff_success == True ]; then echo "PASS: fits files are repeatable"; fi
  
    # vdif check
    vdif_diff_success=True
    diff no_header_1.vdif no_header_2.vdif
    if [ \$? != 0 ]; then vdif_diff_success=False; fi
    diff no_header_3.vdif no_header_2.vdif
    if [ \$? != 0 ]; then vdif_diff_success=False; fi
    diff no_header_1.vdif no_header_3.vdif
    if [ \$? != 0 ]; then vdif_diff_success=False; fi
    if [ \$vdif_diff_success == True ]; then echo "PASS: vdif files are repeatable"; fi
    """

}

workflow repeat_1 {
    main:
        make_beam_ipfb( Channel.from("1150234552"), Channel.from("1150234232"),
                        Channel.from("repeat_1"), Channel.from([["148", "009"]]),
                        Channel.from("2016-06-17T21:46:25"),
                        Channel.from("00:34:21.83_-05:34:36.72"),
                        Channel.of([1150235202, 1150235501]) )
    emit:
        make_beam_ipfb.out[0]
        make_beam_ipfb.out[1]

}

workflow repeat_2 {
    main:
        make_beam_ipfb( Channel.from("1150234552"), Channel.from("1150234232"),
                        Channel.from("repeat_2"), Channel.from([["148", "009"]]),
                        Channel.from("2016-06-17T21:46:25"),
                        Channel.from("00:34:21.83_-05:34:36.72"),
                        Channel.of([1150235202, 1150235501]) )
    emit:
        make_beam_ipfb.out[0]
        make_beam_ipfb.out[1]
}

workflow repeat_3 {
    main:
        make_beam_ipfb( Channel.from("1150234552"), Channel.from("1150234232"),
                        Channel.from("repeat_3"), Channel.from([["148", "009"]]),
                        Channel.from("2016-06-17T21:46:25"),
                        Channel.from("00:34:21.83_-05:34:36.72"),
                        Channel.of([1150235202, 1150235501]) )
    emit:
        make_beam_ipfb.out[0]
        make_beam_ipfb.out[1]
}

workflow repeatability_test {
    repeat_1()
    repeat_2()
    repeat_3()
    compare_repeats( repeat_1.out[0],
                     repeat_1.out[1],
                     repeat_2.out[0],
                     repeat_2.out[1],
                     repeat_3.out[0],
                     repeat_3.out[1] )
}

workflow prepfold_compare {

}

workflow tests_1150234552 {
    make_beam( // obsid
               Channel.from("1150234552"),
               // calid
               Channel.from("1150234232"),
               // label
               Channel.from("tests_1150234552"),
               // coarse channels id and number
               Channel.of(["148", "009"]),
               // utc
               Channel.from("2016-06-17T21:46:25"),
               // pointings
               Channel.of(["00:34:21.83_-05:34:36.72", "00:34:08.87_-07:21:53.40"]),
               // begin, end
               Channel.of([1150235202, 1150235501]) )
    make_beam_ipfb( // obsid
                    Channel.from("1150234552"),
                    // calid
                    Channel.from("1150234232"),
                    // label
                    Channel.from("tests_1150234552_vdif"),
                    // coarse channels id and number
                    Channel.of(["148", "009"]),
                    // utc
                    Channel.from("2016-06-17T21:46:25"),
                    // pointings
                    Channel.of(["00:34:21.83_-05:34:36.72", "00:34:08.87_-07:21:53.40"]),
                    // begin, end
                    Channel.of([1150235202, 1150235501]) )
}

workflow {
    // Test that multiple runs of the beamformer creates the same result
    repeatability_test()
    // J0034-0721 and J0034-0721 tests
    tests_1150234552()
}