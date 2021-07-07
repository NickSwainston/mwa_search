#!/usr/bin/env nextflow

nextflow.preview.dsl = 2
params.no_combined_check = true
include { pre_beamform; beamform } from './beamform_module'

//params.out_dir = "${params.search_dir}/${params.obsid}_toas"
//params.out_dir = "${params.search_dir}/psr2_timing/${params.obsid}_toas"
params.out_dir = "${params.search_dir}/psr2_timing"

params.bins = 128
params.period = ""
params.dm = ""
params.nchan = 48
params.ncchan = 1
params.subint = ""
params.eph = ""
params.dspsr_options = ""

params.chan_split = false
params.time_split = false

params.fits_file = "None"
params.fits_file_dir = "None"
//params.std_profile = "/astro/mwavcs/nswainston/pulsar_timing/1255444104_cand_0.90004_23.1227_archive_24chan_profile.pTP"
params.std_profile = "/astro/mwavcs/pulsar_search/psr2_timing/1274143152_J0024-1932_profile.ar"
params.label = "psr2"

params.help = false
if ( params.help ) {
    help = """candidate_TOAs.nf: A pipeline that will generate pulsar TOAs to be used for timing in tempo2
             |Argurments:
             |  --fits_file The fits file to process. If this is not supplied, the pipeline will beamform
             |              for you if you have the required beamforming arguments (see beamform.nf -h)
             |  --fits_file_dir
             |              A base directory of all fits files to process. Will search subdirectory to
             |              find fits files.
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
             |  --period    The topo period of the pulsar in seconds. 
             |  --dm        The dispersion measure of the pulsar.
             |
             |  --out_dir   Where the TOAs will be output
             |              [default: ${params.out_dir}/<obsid>]
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
    publishDir "${params.out_dir}/${fits_files.baseName.split("_")[0]}", mode: 'copy'
    label 'cpu'
    time '2h'

    input:
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
    dspsr -t $task.cpus -b ${params.bins} -c \${period} -D \${DM} -O ${params.obsid}_b${params.bins}_ch\${chans} -cont -U 4000 ${params.dspsr_options} G*_${params.obsid}*ch\${chans}*.fits
    pam -pTF -e pTDF --name J0036-1033 *.ar
    """
}


process dspsr_time {
    publishDir "${params.out_dir}/${fits_files.baseName.split("_")[0]}", mode: 'copy'
    label 'cpu_any'
    cpus = 1
    time '12h'

    input:
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
    dspsr -t $task.cpus -b ${params.bins} ${eph_command} -L ${subint_command} -e subint -cont -U 4000 ${params.dspsr_options} *fits
    pam -pTF -e pTDF --name ${params.label} *.subint

    # Update file names
    for i in \$(ls); do
        mv \$i ${fits_files.baseName}_\$i
    done
    """
}

process get_toas {
    publishDir "${params.out_dir}/${archive.baseName.split("_")[0]}", pattern: "*ps", mode: 'copy'

    input:
    //each file(archive)
    //file std_profile
    tuple file(archive), file(std_profile)

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

process combine_obs_toas {
    publishDir "${params.out_dir}/${toa_tims_and_subints[0].baseName.split("_")[0]}", mode: 'copy'

    input:
    //file toa_tims
    //file subints
    file toa_tims_and_subints

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
    awk  '/FORMAT 1/&&c++>0 {next} 1' temp.tim > ${toa_tims_and_subints[0].baseName.split("_")[0]}_all.tim
    psradd -f ${toa_tims_and_subints[0].baseName.split("_")[0]}.ar *.subint
    """
}

process combine_all_toas {
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    file toa_tims

    output:
    file "*all.tim"

    if ( "$HOSTNAME".startsWith("farnarkle") ) {
        beforeScript "module use ${params.presto_module_dir}; module load dspsr/master; module load tempo2"
    }
    else if ( "$HOSTNAME".startsWith("x86") || "$HOSTNAME".startsWith("garrawarla") || "$HOSTNAME".startsWith("galaxy") ) {
        container = "file:///${params.containerDir}/dspsr/dspsr.sif"
    }

    """
    cat *tim > temp.tim
    awk  '/FORMAT 1/&&c++>0 {next} 1' temp.tim > ${params.label}_all.tim
    """
}


workflow {
    // Get fits files
    if ( params.fits_file == "None" && params.fits_file_dir == "None" ) {
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
        pre_beamform()
        beamform( pre_beamform.out[0],\
                  pre_beamform.out[1],\
                  pre_beamform.out[2],\
                  pointings )
        fits_files = beamform.out[1].flatten()
    }
    else if ( params.fits_file_dir != "None" ) {
        fits_files = Channel.fromPath(params.fits_file_dir + "/**fits").flatten()
    }
    else {
        fits_files = Channel.fromPath(params.fits_file).flatten()
    }

    // Either split the obs in time of frequency
    if ( params.chan_split ) {
        pre_beamform()
        dspsr_ch( fits_files )
        dspsr_out_pTDF   = dspsr_ch.out[0]
        dspsr_out_subint = dspsr_ch.out[1]
    }
    else if ( params.time_split ) {
        dspsr_time( fits_files )
        dspsr_out_pTDF   = dspsr_time.out[0]
        dspsr_out_subint = dspsr_time.out[1]
    }
    //get_toas( dspsr_out_pTDF.flatten().view(),
    //          std_profile )
    get_toas( dspsr_out_pTDF.flatten().combine(std_profile) )
    combine_obs_toas( get_toas.out[0].flatten().concat(dspsr_out_subint.flatten()).map{ it -> [ it.baseName.split("_ch")[0], it ] }.\
                      groupTuple().map{ it -> it[1] } )
    combine_all_toas( combine_obs_toas.out[0].view().collect() )
}