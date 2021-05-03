nextflow.preview.dsl = 2

params.out_dir = "${params.search_dir}/${params.obsid}_candidates"

params.vcstools_version = 'master'
params.mwa_search_version = 'master'
params.publish_all_prepfold = false

params.begin = 0
params.end = 0
params.all = false

params.dm_min = 1
params.dm_max = 250
params.dm_min_step = 0.02
params.dm_max_step = 0.5
params.max_dms_per_job = 5000

//Defaults for the accelsearch command
params.nharm = 16 // number of harmonics to search
params.min_period = 0.001 // min period to search for in sec (ANTF min = 0.0013)
params.max_period = 30 // max period to search for in sec  (ANTF max = 23.5)
params.zmax = 0

//Some math for the accelsearch command
//convert to freq
min_freq = 1 / params.max_period
max_freq = 1 / params.min_period
//adjust the freq to include the harmonics
min_f_harm = min_freq
max_f_harm = max_freq * params.nharm

// Work out total obs time
if ( params.all ) {
    // an estimation since there's no easy way to make this work
    obs_length = 4805
}
else {
    obs_length = params.end - params.begin + 1
}

// If doing an acceleration search, lower the number of DMs per job so the jobs don't time out
if ( params.zmax == 0 ) {
    total_dm_jobs = 6
}
else {
    total_dm_jobs = 24
}

// Work out some estimated job times
if ( "$HOSTNAME".startsWith("farnarkle") ) {
    // In seconds
    search_dd_fft_acc_dur = obs_length * 5.0
    prepfold_dur = obs_length * 16.0
    presto_python_load = "module use ${params.presto_module_dir}; module load presto/${params.presto_module}; module load python/2.7.14; module load matplotlib/2.2.2-python-2.7.14"
}
else if ( "$HOSTNAME".startsWith("garrawarla") ) {
    // In seconds
    search_dd_fft_acc_dur = obs_length * 5.0
    prepfold_dur = obs_length * 16.0
    presto_python_load = ""
}
else {
    search_dd_fft_acc_dur = 14400
    prepfold_dur = 7200
    presto_python_load = ""
}

process ddplan {
    label 'ddplan'

    input:
    //tuple val(name), val(fits_files) //fits_files is actauly files but I assume this will save me link
    each val(name), val(fits_files)
    tuple val(begin), val(end)

    output:
    file 'DDplan.txt'

    """
    #!/usr/bin/env python3

    from vcstools.catalogue_utils import grab_source_alog
    from mwa_search.dispersion_tools import dd_plan
    import csv

    #obsid_pointing = "${fits_files[0]}".split("/")[-1].split("_ch")[0]
    #print(obsid_pointing)

    if '$name'.startswith('Blind'):
        output = dd_plan(150., 30.72, 3072, 0.1, $params.dm_min, $params.dm_max,
                         min_DM_step=$params.dm_min_step, max_DM_step=$params.dm_max_step,
                         max_dms_per_job=$params.max_dms_per_job)
    else:
        if '$name'.startswith('dm_'):
            dm = float('$name'.split('dm_')[-1].split('_')[0])
        elif '$name'.startswith('FRB'):
            dm = grab_source_alog(source_type='FRB',
                 pulsar_list=['$name'], include_dm=True)[0][-1]
        else:
            # Try RRAT first
            rrat_temp = grab_source_alog(source_type='RRATs',
                        pulsar_list=['$name'.split("_")[0]], include_dm=True)
            if len(rrat_temp) == 0:
                #No RRAT so must be pulsar
                dm = grab_source_alog(source_type='Pulsar',
                     pulsar_list=['$name'.split("_")[0]], include_dm=True)[0][-1]
            else:
                dm = rrat_temp[0][-1]
        dm_min = float(dm) - 2.0
        if dm_min < 1.0:
            dm_min = 1.0
        dm_max = float(dm) + 2.0
        output = dd_plan(150., 30.72, 3072, 0.1, dm_min, dm_max,
                         min_DM_step=$params.dm_min_step, max_DM_step=$params.dm_max_step,
                         max_dms_per_job=$params.max_dms_per_job)
    with open("DDplan.txt", "w") as outfile:
        spamwriter = csv.writer(outfile, delimiter=',')
        for o in output:
            spamwriter.writerow(['${name}'] + o)

    # Upload beam to the MWA pulsar database
    from vcstools.client import upload_beam
    upload_beam('${name}', ${begin}, ${end}, mwa_search_command='${workflow.commandLine}')
    """
}


process search_dd_fft_acc {
    label 'cpu'
    if ( params.zmax == 0 ) {
        time { search_dd_fft_acc_dur * (0.006*Float.valueOf(dm_values[3]) + 1) < 86400 ? \
                   "${search_dd_fft_acc_dur * (0.006*Float.valueOf(dm_values[3]) + 1)}s" :
                   "86400s"}
    }
    else {
        time { 4 * search_dd_fft_acc_dur * (0.006*Float.valueOf(dm_values[3]) + 1) < 86400 ? \
                   "${4 * search_dd_fft_acc_dur * (0.006*Float.valueOf(dm_values[3]) + 1)}s" :
                   "86400s"}
    }
    //Will ignore errors for now because I have no idea why it dies sometimes
    errorStrategy { task.attempt > 2 ? 'ignore' : 'retry' }
    maxRetries 2
    if ( "$HOSTNAME".startsWith("garrawarla") ) {
        maxForks 400
    }
    else {
        maxForks 800
    }

    input:
    tuple val(name), val(dm_values), file(fits_files)

    output:
    tuple val(name), file("*ACCEL_${params.zmax}"), file("*.inf"), file("*.subSpS"), file('*.cand')
    //file "*ACCEL_0" optional true
    //Will have to change the ACCEL_0 if I do an accelsearch

    if ( "$HOSTNAME".startsWith("farnarkle") ) {
        clusterOptions { "--export=NONE --tmp=${ (int) ( 0.08 * obs_length * Float.valueOf(dm_values[3]) / Float.valueOf(dm_values[5]) ) }MB" }
        scratch '$JOBFS'
        beforeScript "module use ${params.module_dir}; module load presto/min_path"
    }
    else if ( "$HOSTNAME".startsWith("x86") ) {
        //scratch '/ssd'
        container = "file:///${params.containerDir}/presto/presto.sif"
    }
    else if ( "$HOSTNAME".startsWith("galaxy") ) {
        container = "file:///${params.containerDir}/presto/presto.sif"
    }
    else if ( "$HOSTNAME".startsWith("garrawarla") ) {
        clusterOptions { "--export=NONE --tmp=${ (int) ( 0.12 * obs_length * Float.valueOf(dm_values[3]) / Float.valueOf(dm_values[5]) ) }MB" }
        scratch '/nvmetmp'
        container = "file:///${params.containerDir}/presto/presto.sif"
    }
    else {
        container = "nickswainston/presto:realfft_docker"
    }


    """
    echo "lowdm highdm dmstep ndms timeres downsamp"
    echo ${dm_values}
    numout=${(int)(obs_length*10000/Float.valueOf(dm_values[5]))}
    if (( \$numout % 2 != 0 )) ; then
        numout=\$(expr \$numout + 1)
    fi
    printf "\\n#Dedispersing the time series at \$(date +"%Y-%m-%d_%H:%m:%S") --------------------------------------------\\n"
    prepsubband -ncpus $task.cpus -lodm ${dm_values[0]} -dmstep ${dm_values[2]} -numdms ${dm_values[3]} -zerodm -nsub ${dm_values[6]} \
-downsamp ${dm_values[5]} -numout \${numout} -o ${name} ${params.obsid}_*.fits
    printf "\\n#Performing the FFTs at \$(date +"%Y-%m-%d_%H:%m:%S") -----------------------------------------------------\\n"
    realfft *dat
    printf "\\n#Performing the periodic search at \$(date +"%Y-%m-%d_%H:%m:%S") ------------------------------------------\\n"
    for i in \$(ls *.dat); do
        accelsearch -ncpus $task.cpus -zmax ${params.zmax} -flo $min_f_harm -fhi $max_f_harm -numharm $params.nharm \${i%.dat}.fft
    done
    ${presto_python_load}
    single_pulse_search.py -p -m 0.5 -b *.dat
    cat *.singlepulse > ${name}_DM${dm_values[0]}-${dm_values[1]}.subSpS
    printf "\\n#Finished at \$(date +"%Y-%m-%d_%H:%m:%S") ----------------------------------------------------------------\\n"
    """
}


process accelsift {
    label 'cpu'
    time '25m'
    errorStrategy 'retry'
    maxRetries 1

    input:
    tuple val(name), file(accel_inf_single_pulse)

    output:
    tuple val(name), file("cands_*greped.txt")

    if ( "$HOSTNAME".startsWith("farnarkle") || "$HOSTNAME".startsWith("x86") ||\
              "$HOSTNAME".startsWith("garrawarla") || "$HOSTNAME".startsWith("galaxy") ) {
        container = "file:///${params.containerDir}/presto/presto.sif"
    }
    else {
        container = "nickswainston/presto:realfft_docker"
    }
    """
    ACCEL_sift.py --file_name ${name}
    if [ -f cands_${name}.txt ]; then
        grep ${name} cands_${name}.txt > cands_${name}_greped.txt
    else
        #No candidates so make an empty file
        touch cands_${name}_greped.txt
    fi
    #cat *.subSpS > ${name}.SpS
    """
}


process single_pulse_searcher {
    label 'cpu_large_mem'
    time '2h'
    publishDir params.out_dir, mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(name), file(sps), file(fits)

    output:
    file "*pdf" optional true
    file "*.SpS"

    if ( "$HOSTNAME".startsWith("farnarkle") || "$HOSTNAME".startsWith("x86") ||\
         "$HOSTNAME".startsWith("garrawarla") || "$HOSTNAME".startsWith("galaxy") ) {
        container = "file:///${params.containerDir}/sps/sps.sif"
    }
    else {
        container = "nickswainston/sps"
    }
    """
    cat *.subSpS > ${name}.SpS
    #-SNR_min 4 -SNR_peak_min 4.5 -DM_cand 1.5 -N_min 3
    single_pulse_searcher.py -fits ${fits} -no_store -plot_name ${name}_sps.pdf ${name}.SpS
    """
}


process prepfold {
    publishDir params.out_dir, mode: 'copy', enabled: params.publish_all_prepfold
    label 'cpu'
    time "${prepfold_dur}s"
    errorStrategy 'retry'
    maxRetries 1

    input:
    tuple val(cand_line), file(cand_file), file(cand_inf), file(fits_files)

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
    //no mask command currently
    //${cand_line.split()[0].substring(0, cand_line.split()[0].lastIndexOf(":")) + '.cand'}
    """
    echo "${cand_line.split()}"
    # Set up the prepfold options to match the ML candidate profiler
    temp_period=${Float.valueOf(cand_line.split()[7])/1000}
    period=\$(printf "%.8f" \$temp_period)
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

    # Work out how many dmfacts to use to search +/- 2 DM
    ddm=`echo "scale=10;0.000241*138.87^2*\${dmstep} / (1/\$period *\$nbins)" | bc`
    ndmfact=`echo "1 + 1/(\$ddm*\$nbins)" | bc`
    echo "ndmfact: \$ndmfact   ddm: \$ddm"

    #-p \$period
    prepfold -ncpus $task.cpus -o ${cand_line.split()[0]} -n \$nbins -dm ${cand_line.split()[1]} -noxwin -noclip -nsub 256 \
-accelfile ${cand_line.split()[0].substring(0, cand_line.split()[0].lastIndexOf(":"))}.cand -accelcand ${cand_line.split()[0].split(":")[-1]} \
-npart \$ntimechunk -dmstep \$dmstep -pstep 1 -pdstep 2 -npfact \$period_search_n -ndmfact \$ndmfact -runavg *.fits
    """
}


process search_dd {
    label 'cpu'
    time '4h'
    //Will ignore errors for now because I have no idea why it dies sometimes
    errorStrategy { task.attempt > 1 ? 'ignore' : 'retry' }
    if ( "$HOSTNAME".startsWith("garrawarla") ) {
        maxForks 300
    }
    else {
        maxForks 800
    }

    input:
    tuple val(name), val(dm_values), file(fits_files)

    output:
    tuple val(name), file("*.inf"), file("*.subSpS")
    //Will have to change the ACCEL_0 if I do an accelsearch

    if ( "$HOSTNAME".startsWith("farnarkle") ) {
        clusterOptions { "--export=NONE --tmp=${ (int) ( 0.08 * obs_length * Float.valueOf(dm_values[3]) / Float.valueOf(dm_values[5]) ) }MB" }
        scratch '$JOBFS'
        beforeScript "module use ${params.module_dir}; module load presto/min_path"
    }
    else if ( "$HOSTNAME".startsWith("x86") ) {
        scratch '/ssd'
        container = "file:///${params.containerDir}/presto/presto.sif"
    }
    else if ( "$HOSTNAME".startsWith("galaxy") ) {
        container = "file:///${params.containerDir}/presto/presto.sif"
    }
    else if ( "$HOSTNAME".startsWith("garrawarla") ) {
        clusterOptions { "--export=NONE --tmp=${ (int) ( 0.08 * obs_length * Float.valueOf(dm_values[3]) / Float.valueOf(dm_values[5]) ) }MB" }
        scratch '/nvmetmp'
        container = "file:///${params.containerDir}/presto/presto.sif"
    }
    else {
        container = "nickswainston/presto:realfft_docker"
    }

    """
    echo "lowdm highdm dmstep ndms timeres downsamp"
    echo ${dm_values}
    printf "\\n#Dedispersing the time series at \$(date +"%Y-%m-%d_%H:%m:%S") --------------------------------------------\\n"
    prepsubband -ncpus $task.cpus -lodm ${dm_values[0]} -dmstep ${dm_values[2]} -numdms ${dm_values[3]} -zerodm -nsub ${dm_values[6]} \
-downsamp ${dm_values[5]} -numout ${(int)(obs_length*10000/Float.valueOf(dm_values[5]))} -o ${name.replaceAll("\\*","")} ${params.obsid}_*.fits
    single_pulse_search.py -p -m 0.5 -b *.dat
    cat *.singlepulse > ${name}_DM${dm_values[0]}-${dm_values[1]}.subSpS
    """
}


workflow pulsar_search {
    take:
        name_fits_files // [val(candidateName_obsid_pointing), file(fits_files)]
    main:
        ddplan( name_fits_files )
        search_dd_fft_acc( // combine the fits files and ddplan with the matching name key (candidateName_obsid_pointing)
                           ddplan.out.splitCsv().map{ it -> [ it[0], [ it[1], it[2], it[3], it[4], it[5], it[6], it[7] ] ] }.\
                           concat(name_fits_files).groupTuple().\
                           // Find for each ddplan match that with the fits files and the name key then change the format to [val(name), val(dm_values), file(fits_files)]
                           map{ it -> [it[1].init(), [[it[0], it[1].last()]]].combinations() }.flatMap().\
                           map{ it -> [it[1][0], it[0], it[1][1]]} )
        // Get all the inf, ACCEL and single pulse files and sort them into groups with the same name key
        accelsift( search_dd_fft_acc.out.map{ it -> [it[0], [it[1]].flatten().findAll { it != null } + \
                                                            [it[2]].flatten().findAll { it != null }] }.\
                   groupTuple( size: total_dm_jobs, remainder: true ).map{ it -> [it[0], it[1].flatten()]} )
        single_pulse_searcher( search_dd_fft_acc.out.map{ it -> [it[0], [it[3]].flatten().findAll { it != null }] }.\
                               groupTuple( size: total_dm_jobs, remainder: true ).map{ it -> [it[0], it[1].flatten()]}.\
                               // Add fits files
                               concat(name_fits_files).groupTuple( size: 2 ).map{ it -> [it[0], it[1][0], it[1][1]]} )
        prepfold( name_fits_files.cross(
                  // Group all the accelsift lines together
                  accelsift.out.map{ it -> it[1] }.splitCsv().flatten().map{ it -> [it.split()[0].split("_ACCEL")[0], it ] }.cross(
                  // Group all the .cand and .inf files by their base names
                  search_dd_fft_acc.out.map{ it -> [it[2]].flatten().findAll { it != null } }.
                  flatten().map{ it -> [it.baseName.split(".inf")[0], it ] }.concat(
                  search_dd_fft_acc.out.map{ it -> [it[4]].flatten().findAll { it != null } }.
                  flatten().map{ it -> [it.baseName.split("_ACCEL")[0], it ] }).groupTuple( size: 2 )
                  // match the cand and inf file with each accelsift line and reoraganise
                  ).map{ it -> [it[0][0].split("_DM")[0], [it[0][1], it[1][1][0], it[1][1][1]]] }
                  // Match with fits files and eogranise to val(cand_line), file(cand_file), file(cand_inf), file(fits_files)
                  ).map{ it -> [it[1][1][0], it[1][1][2], it[1][1][1], it[0][1]] } )
    emit:
        accelsift.out
        prepfold.out
}

workflow single_pulse_search {
    take:
        name_fits_files
    main:
        ddplan( name_fits_files )
        search_dd( // combine the fits files and ddplan witht he matching name key (candidateName_obsid_pointing)
                   ddplan.out.splitCsv().map{ it -> [ it[0], [ it[1], it[2], it[3], it[4], it[5], it[6], it[7] ] ] }.concat(name_fits_files).groupTuple().\
                   // Find for each ddplan match that with the fits files and the name key then change the format to [val(name), val(dm_values), file(fits_files)]
                   map{ it -> [it[1].init(), [[it[0], it[1].last()]]].combinations() }.flatMap().map{ it -> [it[1][0], it[0], it[1][1]]} )
        single_pulse_searcher( search_dd.out.map{ it -> [it[0], [it[1]].flatten().findAll { it != null } + [it[2]].flatten().findAll { it != null }] }.\
                               groupTuple( size: total_dm_jobs, remainder: true).map{ it -> [it[0], it[1].flatten()] }.\
                               // Add fits files
                               concat(name_fits_files).groupTuple( size: 2 ).map{ it -> [it[0], it[1][0], it[1][1]]}  )
        // Get all the inf and single pulse files and sort them into groups with the same basename (obsid_pointing)
    emit:
        single_pulse_searcher.out[0]
        single_pulse_searcher.out[1]
}
