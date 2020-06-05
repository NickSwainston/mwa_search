nextflow.preview.dsl = 2

params.out_dir = "${params.search_dir}/${params.obsid}_candidates"

params.vcstools_version = 'master'
params.mwa_search_version = 'master'

params.begin = 0
params.end = 0
params.all = false

params.dm_min = 1
params.dm_max = 250
params.dm_min_step = 0.02

//Defaults for the accelsearch command
params.nharm = 16 // number of harmonics to search
params.min_period = 0.001 // min period to search for in sec (ANTF min = 0.0013)
params.max_period = 30 // max period to search for in sec  (ANTF max = 23.5)

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

// Work out some estimated job times
if ( "$HOSTNAME".startsWith("galaxy") ) {
    // In seconds
    search_dd_fft_acc_dur = 14400
    prepfold_dur = 7200
    presto_python_load = ""
}
else{
    search_dd_fft_acc_dur = obs_length * 5.0
    prepfold_dur = obs_length * 2.0
    presto_python_load = "module use ${params.presto_module_dir}; module load presto/${params.presto_module}; module load python/2.7.14; module load matplotlib/2.2.2-python-2.7.14"
}

process ddplan {
    label 'ddplan'

    input:
    tuple val(name), val(fits_files) //fits_files is actauly files but I assume this will save me link

    output:
    file 'DDplan.txt'
    
    """
    #!/usr/bin/env python3

    import find_pulsar_in_obs as fpio
    from lfDDplan import dd_plan
    import csv
    
    #obsid_pointing = "${fits_files[0]}".split("/")[-1].split("_ch")[0]
    #print(obsid_pointing)

    if '$name'.startswith('Blind'):
        output = dd_plan(150., 30.72, 3072, 0.1, $params.dm_min, $params.dm_max,
                         min_DM_step=$params.dm_min_step)
    else:
        if '$name'.startswith('FRB'):
            dm = fpio.grab_source_alog(source_type='FRB',
                 pulsar_list=['$name'.split("_")[0]], include_dm=True)[0][-1]
        else:
            # Try RRAT first
            rrat_temp = fpio.grab_source_alog(source_type='RRATs',
                        pulsar_list=['$name'.split("_")[0]], include_dm=True)
            if len(rrat_temp) == 0:
                #No RRAT so must be pulsar
                dm = fpio.grab_source_alog(source_type='Pulsar',
                     pulsar_list=['$name'.split("_")[0]], include_dm=True)[0][-1]
            else:
                dm = rrat_temp[0][-1]
        dm_min = float(dm) - 2.0
        if dm_min < 1.0:
            dm_min = 1.0
        dm_max = float(dm) + 2.0
        output = dd_plan(150., 30.72, 3072, 0.1, dm_min, dm_max,
                         min_DM_step=$params.dm_min_step)
    with open("DDplan.txt", "w") as outfile:
        spamwriter = csv.writer(outfile, delimiter=',')
        for o in output:
            spamwriter.writerow(['${name}'] + o)
    """ 
}


process search_dd_fft_acc {
    if ( "$HOSTNAME".startsWith("farnarkle") ) {
        scratch '$JOBFS'
        clusterOptions { "--export=NONE --tmp=${ (int) ( 0.08 * obs_length * Float.valueOf(dm_values[3]) / Float.valueOf(dm_values[5]) ) }MB" }
    }
    else {
        //container = "nickswainston/presto"
        container = "presto.sif"
        //stageInMode = 'copy'
    }
    label 'cpu'
    time { "${search_dd_fft_acc_dur * (0.006*Float.valueOf(dm_values[3]) + 1)}s" }
    //Will ignore errors for now because I have no idea why it dies sometimes
    errorStrategy { task.attempt > 1 ? 'ignore' : 'retry' }
    maxForks 800

    input:
    tuple val(name), val(dm_values), file(fits_files)

    output:
    tuple val(name), file("*ACCEL_0"), file("*.inf"), file("*.singlepulse")
    //file "*ACCEL_0" optional true
    //Will have to change the ACCEL_0 if I do an accelsearch

    if ( "$HOSTNAME".startsWith("farnarkle") ) {
        beforeScript "module use ${params.module_dir}; module load presto/min_path"
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
        accelsearch -ncpus $task.cpus -zmax 0 -flo $min_f_harm -fhi $max_f_harm -numharm $params.nharm \${i%.dat}.fft
    done
    ${presto_python_load}
    single_pulse_search.py -p *.dat
    printf "\\n#Finished at \$(date +"%Y-%m-%d_%H:%m:%S") ----------------------------------------------------------------\\n"
    """
}


process accelsift {
    if ( "$HOSTNAME".startsWith("galaxy") ) {
        //container = "nickswainston/presto"
        container = "presto.sif"
        //stageInMode = 'copy'
    }
    label 'cpu'
    time '25m'
    publishDir params.out_dir, pattern: "*_singlepulse.tar.gz", mode: 'copy'
    publishDir params.out_dir, pattern: "*_singlepulse.ps", mode: 'copy'

    input:
    tuple val(name), file(accel_inf_single_pulse)

    output:
    tuple val(name), file("cands_*greped.txt"), file("*_singlepulse.tar.gz"), file("*_singlepulse.ps")

    if ( "$HOSTNAME".startsWith("farnarkle") ) {
        beforeScript "module use ${params.presto_module_dir}; module load presto/${params.presto_module};"+\
                     "module use $params.module_dir; module load mwa_search/py2_scripts"
    }

    """
    ACCEL_sift.py --file_name ${name}
    if [ -f cands_${name}.txt ]; then
        grep ${name} cands_${name}.txt > cands_${name}_greped.txt
    else
        #No candidates so make an empty file
        touch cands_${name}_greped.txt
    fi
    single_pulse_search.py *.singlepulse
    tar -czvhf singlepulse.tar.gz *DM*.singlepulse
    mv singlepulse.tar.gz ${name}_singlepulse.tar.gz
    """
}


process prepfold {
    label 'cpu'
    time "${prepfold_dur}s"

    input:
    tuple file(fits_files), val(cand_line)

    output:
    file "*pfd*"

    beforeScript "module use ${params.presto_module_dir}; module load presto/${params.presto_module}"

    //no mask command currently
    """
    echo "${cand_line.split()}"
    # Set up the prepfold options to match the ML candidate profiler
    period=${Float.valueOf(cand_line.split()[7])/1000}
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

    prepfold -ncpus $task.cpus -o ${cand_line.split()[0]} -n \$nbins -noxwin -noclip -p \$period -dm ${cand_line.split()[1]} -nsub 256 \
-npart \$ntimechunk -dmstep \$dmstep -pstep 1 -pdstep 2 -npfact \$period_search_n -ndmfact 1 -runavg *.fits
    """
}


process search_dd {
    if ( "$HOSTNAME".startsWith("farnarkle") ) {
        scratch '$JOBFS'
        clusterOptions { "--tmp=${ (int) ( 0.04 * obs_length * Float.valueOf(dm_values[3]) / Float.valueOf(dm_values[5]) ) }MB" }
    }
    else {
        //container = "nickswainston/presto"
        container = "presto.sif"
        //stageInMode = 'copy'
    }
    label 'cpu'
    time '4h'
    //Will ignore errors for now because I have no idea why it dies sometimes
    errorStrategy { task.attempt > 1 ? 'ignore' : 'retry' }

    input:
    tuple val(name), val(dm_values), file(fits_files)

    output:
    tuple val(name), file("*.inf"), file("*.singlepulse")
    //Will have to change the ACCEL_0 if I do an accelsearch

    if ( "$HOSTNAME".startsWith("farnarkle") ) {
        beforeScript "module use ${params.presto_module_dir}; module load presto/${params.presto_module};"+\
                     "module load python/2.7.14; module load matplotlib/2.2.2-python-2.7.14; module load numpy/1.16.3-python-2.7.14"
    }

    """
    echo "lowdm highdm dmstep ndms timeres downsamp"
    echo ${dm_values}
    printf "\\n#Dedispersing the time series at \$(date +"%Y-%m-%d_%H:%m:%S") --------------------------------------------\\n"
    prepsubband -ncpus $task.cpus -lodm ${dm_values[0]} -dmstep ${dm_values[2]} -numdms ${dm_values[3]} -zerodm -nsub ${dm_values[6]} \
-downsamp ${dm_values[5]} -numout ${(int)(obs_length*10000/Float.valueOf(dm_values[5]))} -o ${name.replaceAll("\\*","")} ${params.obsid}_*.fits
    single_pulse_search.py -p *.dat
    """
}


process assemble_single_pulse {
    if ( "$HOSTNAME".startsWith("galaxy") ) {
        //container = "nickswainston/presto"
        container = "presto.sif"
        //stageInMode = 'copy'
    }
    label 'cpu_backup'
    time '10m'
    publishDir params.out_dir, mode: 'move'

    input:
    tuple val(name), file(inf_single_pulse)

    output:
    tuple val(name), file("*_singlepulse.tar.gz"), file("*_singlepulse.ps")

    if ( "$HOSTNAME".startsWith("farnarkle") ) {
        beforeScript "module use ${params.presto_module_dir}; module load presto/${params.presto_module};"+\
                     "module load python/2.7.14; module load matplotlib/2.2.2-python-2.7.14; module load numpy/1.16.3-python-2.7.14"
    }

    """
    single_pulse_search.py *.singlepulse
    tar -czvhf singlepulse.tar.gz *DM*.singlepulse
    mv singlepulse.tar.gz ${name.replaceAll("\\*","")}_singlepulse.tar.gz
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
                                                            [it[2]].flatten().findAll { it != null } + \
                                                            [it[3]].flatten().findAll { it != null }] }.\
                   groupTuple( size: 6, remainder: true ).map{ it -> [it[0], it[1].flatten()]} )
        // Make a pair of accelsift out lines and fits files that match
        prepfold( name_fits_files.cross(accelsift.out.map{ it -> it[1] }.splitCsv().flatten().map{ it -> [it.split()[0].split("_DM")[0], it ] }).\
                  map{ it -> [it[0][1], it[1][1]] } )
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
        // Get all the inf and single pulse files and sort them into groups with the same basename (obsid_pointing)
        assemble_single_pulse( search_dd.out.map{ it -> [it[0], [it[1]].flatten().findAll { it != null } + [it[2]].flatten().findAll { it != null }] }.\
                               groupTuple( size: 6, remainder: true).map{ it -> [it[0], it[1].flatten()] } )
    emit:
        assemble_single_pulse.out
}
