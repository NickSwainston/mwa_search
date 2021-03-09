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

// Option input checks
def folder = new File(params.version_compare_dir)
if ( folder.exists() && params.publish_version ) {
    println("This version has already been published in ${params.version_compare_dir}. Please label your version differently if you'd like to publish it")
    exit(0)
}

process print_version {
    echo true
    """
    echo "Comparing the beamformer in branch ${params.vcstools_version} with the previous run labeled ${params.version_compare}."
    """
}

process make_beam {
    publishDir "${params.version_compare_dir}", mode: 'copy', enabled: params.publish_version

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
-f ${channel_pair[0]} -J ${params.version_compare_base_dir}/${obsid}/cal/${calid}/rts/DI_JonesMatrices_node${channel_pair[1]}.dat \
-d ${params.version_compare_base_dir}/${obsid}/combined -P ${point.join(",")} \
-r 10000 -m ${params.version_compare_base_dir}/${obsid}/${obsid}_metafits_ppds.fits \
${bf_out} -t 6000 -F ${params.version_compare_base_dir}/${obsid}/cal/${calid}/rts/flagged_tiles.txt -z $utc
    
    # Label all outputs
    for i in \$(ls */*); do mv \${i} ${label}_\${i##*/}; done
    """
}


process make_beam_ipfb {
    publishDir "${params.version_compare_dir}", mode: 'copy', enabled: params.publish_version

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
    file "*hdr"

    """
    make_beam -o $obsid -b $begin -e $end -a 128 -n 128 \
-f ${channel_pair[0]} -J ${params.version_compare_base_dir}/${obsid}/cal/${calid}/rts/DI_JonesMatrices_node${channel_pair[1]}.dat \
-d ${params.version_compare_base_dir}/${obsid}/combined -P ${point} \
-r 10000 -m ${params.version_compare_base_dir}/${obsid}/${obsid}_metafits_ppds.fits \
-p -v -t 6000 -F ${params.version_compare_base_dir}/${obsid}/cal/${calid}/rts/flagged_tiles.txt -z $utc
    # Label all outputs
    for i in \$(ls */*); do mv \${i} ${label}_\${i##*/}; done
    vdif_file=\$(ls *vdif)
    mv \${vdif_file} ${label}_\${vdif_file}
    mv \${vdif_file%.vdif}.hdr ${label}_\${vdif_file%.vdif}.hdr
    # Update header DATAFILE name
    sed -i "s/\${vdif_file}/${label}_\${vdif_file}/" ${label}_\${vdif_file%.vdif}.hdr
    """
}

process compare_repeats {
    echo true
    errorStrategy { task.exitStatus=2 ? 'ignore' : 'terminate' }

    input:
    file fits_1
    file vdif_1
    file fits_2
    file vdif_2
    file fits_3
    file vdif_3

    """
    # Remove headers from fits files
    tail -n +2 ${fits_1} > no_header_1.fits
    tail -n +2 ${fits_2} > no_header_2.fits
    tail -n +2 ${fits_3} > no_header_3.fits

    # Fits check
    fits_diff_success=True
    diff no_header_1.fits no_header_2.fits
    if [ \$? != 0 ]; then fits_diff_success=False; fi
    diff no_header_3.fits no_header_2.fits
    if [ \$? != 0 ]; then fits_diff_success=False; fi
    diff no_header_1.fits no_header_3.fits
    if [ \$? != 0 ]; then fits_diff_success=False; fi
    if [ \$fits_diff_success == True ]; then echo "PASS: fits files are repeatable 3 times"; fi
  
    # vdif check
    vdif_diff_success=True
    diff ${vdif_1} ${vdif_2}
    if [ \$? != 0 ]; then vdif_diff_success=False; fi
    diff ${vdif_3} ${vdif_2}
    if [ \$? != 0 ]; then vdif_diff_success=False; fi
    diff ${vdif_1} ${vdif_3}
    if [ \$? != 0 ]; then vdif_diff_success=False; fi
    if [ \$vdif_diff_success == True ]; then echo "PASS: vdif files are repeatable 3 times"; fi
    """

}

process prepfold_and_compare {
    echo true
    publishDir "${params.version_compare_dir}", mode: 'copy', enabled: params.publish_version

    label 'cpu'
    time "1h"

    input:
    tuple val(label), val(pulsar), file(fits)

    output:
    file "*pfd*"

    if ( "$HOSTNAME".startsWith("farnarkle") || "$HOSTNAME".startsWith("galaxy") ) {
        beforeScript "module use ${params.presto_module_dir}; module load presto/${params.presto_module}"
    }
    else if ( "$HOSTNAME".startsWith("x86") || "$HOSTNAME".startsWith("garrawarla") ) {
        container = "file:///${params.containerDir}/presto/presto.sif"
    }
    else {
        container = "nickswainston/presto:realfft_docker"
    }

    script:
    if ( params.publish_version )
        """
        psrcat -e ${pulsar} | grep -v TCB > ${pulsar}.eph
        prepfold -n 100 -noxwin -noclip -o ${label}_${pulsar} -nsub 8 -timing ${pulsar}.eph \
    -npart 120 -npfact 1 -pstep 1 -pdstep 2 -ndmfact 1 -runavg -noxwin *.fits  &> prepfold.out
        """
    else
        """
        psrcat -e ${pulsar} | grep -v TCB > ${pulsar}.eph
        prepfold -n 100 -noxwin -noclip -o ${label}_${pulsar} -nsub 8 -timing ${pulsar}.eph \
    -npart 120 -npfact 1 -pstep 1 -pdstep 2 -ndmfact 1 -runavg -noxwin *.fits  &> prepfold.out

        # Get sigma from current detection
        sigma_line=\$(grep sigma *prof)
        sigma_line_2=\${sigma_line#*~}
        current_sn=\${sigma_line_2% sigma)}

        # Get sigma from comparison detection
        sigma_line=\$(grep sigma ${params.version_compare_dir}/${label}_${pulsar}*prof)
        sigma_line_2=\${sigma_line#*~}
        compare_sn=\${sigma_line_2% sigma)}

        # Compare the signal to noise ratios
        if [ \$current_sn == \$compare_sn ]; then
            echo "PASS: ${label} prepfold detections of ${pulsar} have the same SN"
        else
            echo "WARN: ${label} prepfold detections of ${pulsar} have different SN!"
            echo "      Current SN: \$current_sn"
            echo "      ${params.version_compare} SN: \$compare_sn"
        fi
        """
}

process pdmp_and_compare {
    echo true
    publishDir "${params.version_compare_dir}", mode: 'copy', enabled: params.publish_version

    label 'cpu'
    time "1h"
    errorStrategy 'retry'
    maxRetries 1
    //stageInMode 'copy'

    input:
    tuple val(label), val(pulsar), file(vdif), file(hdr)

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

    script:
    if ( params.publish_version )
        """
        psrcat -e ${pulsar} | grep -v TCB > ${pulsar}.eph
        dspsr -t ${task.cpus} -b 100 -E ${pulsar}.eph -L 30 -e subint -cont -U 4000 *.hdr &> dspsr.out
        psradd *.subint -o ${label}_${pulsar}.ar &> prsadd.out
        pam --setnchn 32 -m ${label}_${pulsar}.ar &> pam.out
        pdmp -g ${label}_${pulsar}_pdmp.ps/cps ${label}_${pulsar}.ar &> pdmp.out
        mv pdmp.posn ${label}_${pulsar}_pdmp.posn
        """
    else
        """
        psrcat -e ${pulsar} | grep -v TCB > ${pulsar}.eph
        dspsr -t ${task.cpus} -b 100 -E ${pulsar}.eph -L 30 -e subint -cont -U 4000 *.hdr &> dspsr.out
        psradd *.subint -o ${label}_${pulsar}.ar &> prsadd.out
        pam --setnchn 32 -m ${label}_${pulsar}.ar &> pam.out
        pdmp -g ${label}_${pulsar}_pdmp.ps/cps ${label}_${pulsar}.ar &> pdmp.out
        mv pdmp.posn ${label}_${pulsar}_pdmp.posn

        # Get sigma from current detection
        current_sn=\$(cat *.posn | cut -d ' ' -f 7)

        # Get sigma from comparison detection
        compare_sn=\$(cat ${params.version_compare_dir}/${label}_${pulsar}*.posn | cut -d ' ' -f 7)

        # Compare the signal to noise ratios
        if [ \$current_sn == \$compare_sn ]; then
            echo "PASS: ${label} pdmp detections of ${pulsar} have the same SN"
        else
            echo "WARN: ${label} pdmp detections of ${pulsar} have different SN!"
            echo "      Current SN: \$current_sn"
            echo "      ${params.version_compare} SN: \$compare_sn"
        fi
        """
}

process dspsr {
    publishDir "${params.version_compare_dir}", mode: 'copy', enabled: params.publish_version

    label 'cpu'
    time "1h"
    errorStrategy 'retry'
    maxRetries 1

    input:
    tuple val(label), val(pulsar), file(vdif), file(hdr)

    output:
    file "*ascii"

    if ( "$HOSTNAME".startsWith("farnarkle") ) {
        beforeScript "module use ${params.presto_module_dir}; module load dspsr/master"
    }
    else if ( "$HOSTNAME".startsWith("x86") || "$HOSTNAME".startsWith("garrawarla") || "$HOSTNAME".startsWith("galaxy") ) {
        container = "file:///${params.containerDir}/dspsr/dspsr.sif"
    }
    else {
        container = "nickswainston/dspsr_docker"
    }

    """
    psrcat -e ${pulsar} | grep -v TCB > ${pulsar}.eph
    dspsr -O ${label}_${pulsar} -t ${task.cpus} -b 100 -E ${pulsar}.eph -cont -U 4000 *.hdr
    pdv -tK ${label}_${pulsar}.ar > ${label}_${pulsar}.ascii
    """
}


process compare_ascii {
    echo true

    when:
    params.publish_version == false

    input:
    file ascii

    """
    ascii_file=${ascii}
    label_pulsar=\${ascii_file%*.ascii}
    period=\$(psrcat -e2 \${label_pulsar#*_} | grep P0 | cut -d ' ' -f15)

    # Get sigma from current detection
    prof_estimate.py --ascii ${ascii} --auto --period \$period &> prof_utils.out
    current_sn=\$(grep estimate prof_utils.out | cut -d ':' -f 6 | cut -d '+' -f 1 | xargs)

    # Get sigma from comparison detection
    prof_estimate.py --ascii ${params.version_compare_dir}/${ascii} --auto --period \$period &> prof_utils.out
    compare_sn=\$(grep estimate prof_utils.out | cut -d ':' -f 6 | cut -d '+' -f 1 | xargs)

    # Compare the signal to noise ratios
    if [ \$current_sn == \$compare_sn ]; then
        echo "PASS: \${label_pulsar%_*} dspsr detections of \${label_pulsar#*_} have the same SN"
    else
        echo "WARN: \${label_pulsar%_*} dspsr detections of \${label_pulsar#*_} have different SN!"
        echo "      Current SN: \$current_sn"
        echo "      ${params.version_compare} SN: \$compare_sn"
    fi
        """
}

process compare_fits {
    echo true

    input:
    tuple val(label), val(pulsar), file(fits)

    """
    # Compare fits files binaries
    tail -n +2 ${fits} > no_header.fits
    tail -n +2 ${params.version_compare_dir}/${fits} > no_header_compare.fits
    trap 'if [[ \$? == 2 ]]; then echo "WARN: ${label} ${pulsar} fits files differ!"; exit 0; else echo "PASS: ${label} ${pulsar} fits files do not differ"; fi' EXIT
    diff no_header.fits no_header_compare.fits &> fits_diff.txt
    """
}

process compare_vdif {
    echo true

    input:
    tuple val(label), val(pulsar), file(vdif), file(hdr)

    """
    # Compare fits files binaries
    trap 'if [[ \$? == 2 ]]; then echo "WARN: ${label} ${pulsar} vdif files differ!"; exit 0; else echo "PASS: ${label} ${pulsar} vdif files do not differ"; fi' EXIT
    diff ${vdif} ${params.version_compare_dir}/${vdif} &> vdif_diff.txt
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
    prepfold_and_compare( // [label, pulsar, fits]
                          Channel.of(["00:34:21.83_-05:34:36.72", "1150234552"],
                                     ["00:34:08.87_-07:21:53.40", "1150234552"],
                                     ["00:34:21.83_-05:34:36.72", "J0034-0534"],
                                     ["00:34:08.87_-07:21:53.40", "J0034-0721"]).concat(\
                          make_beam.out.flatten().map { it -> [it.baseName.split("_ch")[0].split("1150234552_")[-1], it ] }).\
                          groupTuple().map{ it -> [it[1][0], it[1][1], it[1][2]]})
    compare_fits( // [label, pulsar, fits]
                  Channel.of(["00:34:21.83_-05:34:36.72", "1150234552"],
                             ["00:34:08.87_-07:21:53.40", "1150234552"],
                             ["00:34:21.83_-05:34:36.72", "J0034-0534"],
                             ["00:34:08.87_-07:21:53.40", "J0034-0721"]).concat(\
                  make_beam.out.flatten().map { it -> [it.baseName.split("_ch")[0].split("1150234552_")[-1], it ] }).\
                  groupTuple().map{ it -> [it[1][0], it[1][1], it[1][2]]})
    // We don't test J0034-0534 because it's not bright enough to fit a good SN with dspsr
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
                    Channel.of(["00:34:08.87_-07:21:53.40"]).flatten(),
                    // begin, end
                    Channel.of([1150235202, 1150235501]) )
    dspsr( // [label, pulsar, vdif, hdr]
           Channel.of(["00:34:08.87_-07:21:53.40", "1150234552"],
                      ["00:34:08.87_-07:21:53.40", "J0034-0721"]).concat(\
           make_beam_ipfb.out[1].flatten().concat(make_beam_ipfb.out[2].flatten()).\
           map { it -> [it.baseName.split("_ch")[0].split("1150234552_")[-1], it ] }).\
           groupTuple().map{ it -> [it[1][0], it[1][1], it[1][2], it[1][3]]} )
    compare_ascii( dspsr.out.flatten() )
    compare_vdif( // [label, pulsar, vdif, hdr]
                  Channel.of(["00:34:08.87_-07:21:53.40", "1150234552"],
                             ["00:34:08.87_-07:21:53.40", "J0034-0721"]).concat(\
                  make_beam_ipfb.out[1].flatten().concat(make_beam_ipfb.out[2].flatten()).\
                  map { it -> [it.baseName.split("_ch")[0].split("1150234552_")[-1], it ] }).\
                  groupTuple().map{ it -> [it[1][0], it[1][1], it[1][2], it[1][3]]} )
}

workflow tests_1260302416 {
    make_beam( // obsid
               Channel.from("1260302416"),
               // calid
               Channel.from("1260305032"),
               // label
               Channel.from("tests_1260302416"),
               // coarse channels id and number
               Channel.of(["119", "011"]),
               // utc
               Channel.from("2019-12-13T20:06:22"),
               // pointings
               Channel.of(["09:53:09.31_+07:55:35.75"]),
               // begin, end
               Channel.of([1260302800, 1260302899]) )
    prepfold_and_compare( // [label, pulsar, fits]
                          Channel.of(["09:53:09.31_+07:55:35.75", "1260302416"],
                                     ["09:53:09.31_+07:55:35.75", "J0953+0755"],).concat(\
                          make_beam.out.flatten().map { it -> [it.baseName.split("_ch")[0].split("1260302416_")[-1], it ] }).\
                          groupTuple().map{ it -> [it[1][0], it[1][1], it[1][2]]})
    make_beam_ipfb( // obsid
                    Channel.from("1260302416"),
                    // calid
                    Channel.from("1260305032"),
                    // label
                    Channel.from("tests_1260302416_vdif"),
                    // coarse channels id and number
                    Channel.of(["119", "011"]),
                    // utc
                    Channel.from("2019-12-13T20:06:22"),
                    // pointings
                    Channel.of(["09:53:09.31_+07:55:35.75"]).flatten(),
                    // begin, end
                    Channel.of([1260302800, 1260302899]) )
    dspsr( // [label, pulsar, fits]
           Channel.of(["09:53:09.31_+07:55:35.75", "1260302416"],
                       ["09:53:09.31_+07:55:35.75", "J0953+0755"]).concat(\
           make_beam_ipfb.out[1].flatten().concat(make_beam_ipfb.out[2].flatten()).\
           map { it -> [it.baseName.split("_ch")[0].split("1260302416_")[-1], it ] }).\
           groupTuple().map{ it -> [it[1][0], it[1][1], it[1][2], it[1][3]]} )
    compare_ascii( dspsr.out.flatten() )
}

workflow {
    print_version()
    // Test that multiple runs of the beamformer creates the same result
    repeatability_test()
    // J0034-0534 and J0034-0721 tests
    tests_1150234552()
    // J0953+0755 test
    tests_1260302416()
}