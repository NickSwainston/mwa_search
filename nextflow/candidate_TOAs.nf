#!/usr/bin/env nextflow

nextflow.preview.dsl = 2
include { pre_beamform; beamform } from './beamform_module'

params.obsid = null
params.calid = null
params.pointings = null
params.pointing_file = null

params.begin = null
params.end = null
params.all = false

params.summed = true
params.vcstools_version = 'master'
params.mwa_search_version = 'master'

params.didir = "${params.scratch_basedir}/${params.obsid}/cal/${params.calid}/rts"
params.channels = null
params.out_dir = "${params.search_dir}/${params.obsid}_toas"

params.bins = 128
params.period = 0.90004
params.dm = 23.123
params.nchan = 48
params.ncchan = 1
params.subint = ""
params.eph = ""

params.no_beamform = false
params.no_combined_check = false

std_profile = Channel.fromPath("/group/mwavcs/nswainston/pulsar_timing/1255444104_cand_0.90004_23.1227_archive_24chan_profile.pTP")

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
else {
    println "No pointings given. Either use --pointing_file or --pointings. Exiting"
    exit(1)
}

if ( ! ( params.chan_split || params.time_split) ) {
    println "Please use either --chan_split or --time_split"
    exit(1)
}

if ( params.subint == "" ) {
    subint_command = '\$(python -c "print(\'{:d}\'.format(int((8.0/\$sn)**2*\$samples/10000)))")'
}
else {
    subint_command = params.subint
}
if ( params.eph == "" ) {
    eph_command = "-c \${period} -D \${DM}"
}
else {
    eph_command = " -E ${params.eph}"
}

process prepfold_ch {
    label 'cpu'
    time '2h'

    input:
    file fits_files
    each chans

    output:
    file "*bestprof"
    file fits_files
    
    if ( "$HOSTNAME".startsWith("farnarkle") ) {
        beforeScript "module use ${params.presto_module_dir}; module load presto/${params.presto_module}"
    }
    else if ( "$HOSTNAME".startsWith("x86") || "$HOSTNAME".startsWith("garrawarla") || "$HOSTNAME".startsWith("galaxy") ) {
        container = "${config.containDir}/presto/presto.sif"
    }
    else {
        container = "nickswainston/presto:realfft_docker"
    }

    //no mask command currently
    """
    prepfold  -o pulsar_timing_check_ch${chans[0]} -n ${params.bins} -noxwin -noclip -p ${params.period} -dm ${params.dm} -nsub ${params.ncchan*8} -npart 120 \
-dmstep 1 -pstep 1 -pdstep 2 -npfact 1 -ndmfact 1 -runavg G*_${params.obsid}*ch${chans[0]}*.fits
    """
}

process dspsr_ch {
    label 'cpu'
    time '2h'

    input:
    each file(bestprof)
    file fits_files

    output:
    file "*pTDF"

    if ( "$HOSTNAME".startsWith("farnarkle") ) {
        beforeScript "module use ${params.presto_module_dir}; module load dspsr/master"
    }
    else if ( "$HOSTNAME".startsWith("x86") || "$HOSTNAME".startsWith("garrawarla") || "$HOSTNAME".startsWith("galaxy") ) {
        container = "${config.containDir}/presto/presto.sif"
    }
    else {
        container = "nickswainston/presto:realfft_docker"
    }

    //may need to add some channel names
    """
    chans=\$(ls *.bestprof | cut -d 'h' -f 3 | cut -d '_' -f 1)
    echo "chans: \$chans"
    DM=\$(grep DM *.bestprof | tr -s ' ' | cut -d ' ' -f 5)
    echo "DM: \$DM"
    period=\$(grep P_topo *.bestprof | tr -s ' ' | cut -d ' ' -f 5)
    period="\$(echo "scale=10;\${period}/1000"  |bc)"
    echo "period: \$period"
    dspsr -t $task.cpus -b ${params.bins} -c \${period} -D \${DM} -O ${params.obsid}_b${params.bins}_ch\${chans} -cont -U 4000 G*_${params.obsid}*ch\${chans}*.fits
    pam -pTF -e pTDF --name J0036-1033 *.ar
    """
}

process prepfold_time {
    label 'cpu'
    time '2h'

    input:
    file fits_files

    output:
    file "*bestprof"

    if ( "$HOSTNAME".startsWith("farnarkle") ) {
        beforeScript "module use ${params.presto_module_dir}; module load presto/${params.presto_module}"
    }
    else if ( "$HOSTNAME".startsWith("x86") || "$HOSTNAME".startsWith("garrawarla") || "$HOSTNAME".startsWith("galaxy") ) {
        container = "${config.containDir}/presto/presto.sif"
    }
    else {
        container = "nickswainston/presto:realfft_docker"
    }

    //no mask command currently
    """
    prepfold -o pulsar_timing_check -n ${params.bins} -noxwin -noclip -p ${params.period} -dm ${params.dm} -nsub 256 -npart 120 \
-dmstep 1 -pstep 1 -pdstep 2 -npfact 1 -ndmfact 1 -runavg ${params.obsid}*fits
    """
}


process dspsr_time {
    cpus = 5
    label 'cpu'
    time '6h'

    input:
    file bestprof
    file fits_files

    output:
    file "*pTDF"

    if ( "$HOSTNAME".startsWith("farnarkle") ) {
        beforeScript "module use ${params.presto_module_dir}; module load dspsr/master"
    }
    else if ( "$HOSTNAME".startsWith("x86") || "$HOSTNAME".startsWith("garrawarla") || "$HOSTNAME".startsWith("galaxy") ) {
        container = "${config.containDir}/presto/presto.sif"
    }
    else {
        container = "nickswainston/presto:realfft_docker"
    }

    //may need to add some channel names
    """
    sn="\$(grep sigma *.bestprof | tr -s ' ' | cut -d ' ' -f 5 | cut -d '~' -f 2)"
    samples="\$(grep "Data Folded" *.bestprof | tr -s ' ' | cut -d ' ' -f 5)"
    period=\$(grep P_topo *.bestprof | tr -s ' ' | cut -d ' ' -f 5)
    period="\$(echo "scale=10;\${period}/1000"  |bc)"
    DM=\$(grep DM *.bestprof | tr -s ' ' | cut -d ' ' -f 5)
    echo "DM: \$DM   Period: \$period   SN: \$sn"
    dspsr -t $task.cpus -b ${params.bins} ${eph_command} -L ${subint_command} -e subint -cont -U 4000 ${params.obsid}*fits
    pam -pTF -e pTDF --name J0036-1033 *.subint
    """
}

process get_toas {
    publishDir params.out_dir, pattern: "*ps", mode: 'copy'

    input:
    each file(archive)
    file std_profile

    output:
    file "*tim"
    file "*ps"

    beforeScript "module use ${params.presto_module_dir}; module load dspsr/master; module load tempo2"

    """
    archive_name=${archive}
    pat -s ${std_profile} ${archive} -f tempo2 > \${archive_name%pTDF}tim
    pav -CDFTp -g \${archive_name%pTDF}ps/cps ${archive}
    """
}

process combine_toas {
    publishDir params.out_dir, mode: 'copy'

    input:
    file toa_tims

    output:
    file "*all.tim"

    """
    cat *tim > temp.tim
    awk  '/FORMAT 1/&&c++>0 {next} 1' temp.tim > ${params.obsid}_all.tim
    """
}


workflow {
    pre_beamform()
    if ( params.no_beamform ) {
        fits_files = Channel.fromPath("${params.basedir}/${params.obsid}/pointings/${params.pointings}/${params.obsid}*fits").collect()
    }
    else {
        beamform( pre_beamform.out[0],\
                  pre_beamform.out[1],\
                  pre_beamform.out[2],\
                  pointings )
        beamform.out[1].collect().set{ fits_files }
    }
    if ( params.chan_split ) {
        prepfold_ch( beamform.out[0],\
                  pre_beamform.out[1].flatten().collate( params.ncchan ) )
        dspsr_ch( prepfold.out[0],\
                  beamform.out[0] )
        get_toas( dspsr_ch.out,\
                  std_profile )
    }
    else if ( params.time_split ) {
        prepfold_time( fits_files )
        dspsr_time( prepfold_time.out[0],\
                    fits_files.collect() )
        get_toas( dspsr_time.out,\
                    std_profile )
    }
    combine_toas( get_toas.out[0].collect() )
}