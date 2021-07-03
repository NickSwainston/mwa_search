#!/usr/bin/env nextflow

nextflow.preview.dsl = 2
params.no_combined_check = true
include { pre_beamform; beamform } from './beamform_module'

//params.out_dir = "${params.search_dir}/${params.obsid}_toas"
params.out_dir = "${params.search_dir}/psr2_timing/${params.obsid}_toas"

params.bins = 128
params.a_period = 1.30624302355253
params.a_dm = 20.703
params.period = ""
params.dm = ""
params.nchan = 48
params.ncchan = 1
params.subint = ""
params.eph = ""

params.fits_file = "None"
//params.std_profile = "/astro/mwavcs/nswainston/pulsar_timing/1255444104_cand_0.90004_23.1227_archive_24chan_profile.pTP"
params.std_profile = "/astro/mwavcs/pulsar_search/psr2_timing/1274143152_J0024-1932_profile.ar"

params.help = false
if ( params.help ) {
    help = """candidate_TOAs.nf: A pipeline that will generate pulsar TOAs to be used for timing in tempo2
             |Argurments:
             |  --obsid     Observation ID you want to process [no default]
             |  --fits_file The fits file to process. If this is not supplied, the pipeline will beamform
             |              for you if you have the required beamforming arguments (see beamform.nf -h)
             |  --std_profile
             |              The standard profile to convolve your pulses with to get the TOAs
             |              [default: ${params.std_profile}]
             |
             | --time_split Split the observation in time of size --subint to get several TOAs
             |  --subint    Size in seconds to split the observation into if using --time_split.
             |              If not supplied will make a guesstimate of a reasonable subint with the presto SN.
             |
             | --chan_split Split the observation in frequency of size --nchan coarse channels
             |              This method isn't well test and may break
             |  --ncchan    Number of coarse channels to split the obs into if using --chan_split.
             |              [default: 1]
             |
             |  --eph       The ephermis file to fold on the obs. If you don't have one use --period and --dm
             |  --a_period  The approximate period of the pulsar. Prepfold will be used find a more accurate one.
             |  --a_dm      The approximate dispersion measure of the pulsar . Prepfold will be used find a more accurate one.
             |  --period    The topo period of the pulsar in seconds. 
             |  --dm        The dispersion measure of the pulsar.
             |
             |  --out_dir   Where the TOAs will be output
             |              [default: ${params.out_dir}]
             |  -w          The Nextflow work directory. Delete the directory once the processs
             |              is finished [default: ${workDir}]""".stripMargin()
    println(help)
    exit(0)
}

std_profile = Channel.fromPath(params.std_profile)

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
    if ( params.period == "" && params.dm == "" ) {
        // approx period and dm from bestprof
        eph_command = "-c \${period} -D \${DM}"
    }
    else {
        // input period and dm
        eph_command = "-c ${params.period} -D ${params.dm}"
    }
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
        container = "file:///${params.containerDir}/presto/presto.sif"
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
        container = "file:///${params.containerDir}/presto/presto.sif"
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
        container = "file:///${params.containerDir}/presto/presto.sif"
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
    publishDir params.out_dir, mode: 'copy'
    label 'cpu_any'
    cpus = 10
    time '6h'

    input:
    file bestprof
    file fits_files

    output:
    file "*pTDF"
    file "*subint"

    if ( "$HOSTNAME".startsWith("farnarkle") ) {
        beforeScript "module use ${params.presto_module_dir}; module load dspsr/master"
    }
    else if ( "$HOSTNAME".startsWith("x86") || "$HOSTNAME".startsWith("garrawarla") || "$HOSTNAME".startsWith("galaxy") ) {
        container = "file:///${params.containerDir}/dspsr/dspsr.sif"
    }
    //may need to add some channel names
    """
    sn="\$(grep sigma *.bestprof | tr -s ' ' | cut -d ' ' -f 5 | cut -d '~' -f 2)"
    samples="\$(grep "Data Folded" *.bestprof | tr -s ' ' | cut -d ' ' -f 5)"
    period=\$(grep P_topo *.bestprof | tr -s ' ' | cut -d ' ' -f 5)
    period="\$(echo "scale=10;\${period}/1000"  |bc)"
    DM=\$(grep DM *.bestprof | tr -s ' ' | cut -d ' ' -f 5)
    echo "DM: \$DM   Period: \$period   SN: \$sn"
    dspsr -t $task.cpus -b ${params.bins} ${eph_command} -L ${subint_command} -e subint -cont -U 600 ${params.obsid}*fits
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

    if ( "$HOSTNAME".startsWith("farnarkle") ) {
        beforeScript "module use ${params.presto_module_dir}; module load dspsr/master; module load tempo2"
    }
    else if ( "$HOSTNAME".startsWith("x86") || "$HOSTNAME".startsWith("garrawarla") || "$HOSTNAME".startsWith("galaxy") ) {
        container = "file:///${params.containerDir}/dspsr/dspsr.sif"
    }

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
    file subints

    output:
    file "*all.tim"
    file "*.ar"

    if ( "$HOSTNAME".startsWith("farnarkle") ) {
        beforeScript "module use ${params.presto_module_dir}; module load dspsr/master; module load tempo2"
    }
    else if ( "$HOSTNAME".startsWith("x86") || "$HOSTNAME".startsWith("garrawarla") || "$HOSTNAME".startsWith("galaxy") ) {
        container = "file:///${params.containerDir}/dspsr/dspsr.sif"
    }

    """
    cat *tim > temp.tim
    awk  '/FORMAT 1/&&c++>0 {next} 1' temp.tim > ${params.obsid}_all.tim
    psradd -f ${params.obsid}.ar *.subint
    """
}


workflow {
    // Get fits files
    pre_beamform()
    if ( params.fits_file == "None" ) {
        // Send off beamforming
        if ( params.pointings ) {
            pointings = Channel
                .from(params.pointings.split(","))
                .collect()
                .flatten()
                .collate( params.max_pointings )
        }
        else {
            println "No pointings given. Either use --pointings. Exiting"
            exit(1)
        }
        beamform( pre_beamform.out[0],\
                  pre_beamform.out[1],\
                  pre_beamform.out[2],\
                  pointings )
        fits_files = beamform.out[1].collect()
    }
    else {
        fits_files = Channel.fromPath(params.fits_file).collect()
    }

    // Either split the obs in time of frequency
    if ( params.chan_split ) {
        prepfold_ch( beamform.out[0],\
                     pre_beamform.out[1].flatten().collate( params.ncchan ) )
        dspsr_ch( prepfold.out[0],\
                  beamform.out[0] )
        get_toas( dspsr_ch.out,\
                  std_profile )
    }
    else if ( params.time_split ) {
        if ( params.subint == "" || ( params.eph == "" && params.period == "" && params.dm == "" ) ) {
            prepfold_time( fits_files )
            bestprof_file = prepfold_time.out[0]
        }
        else {
            // dummy bestprof
            bestprof_file = Channel.fromPath("/astro/mwavcs/nswainston/J0036-1033_detections/1275085816_00:36:11.58_-10:33:56.44_900.04ms_Cand.pfd.bestprof")
        }
        dspsr_time( bestprof_file,\
                    fits_files )
        get_toas( dspsr_time.out[0].flatten(),
                  std_profile )
    }
    combine_toas( get_toas.out[0].collect(),
                  dspsr_time.out[1].collect() )
}