nextflow.preview.dsl = 2

params.obsid = 1253471952
params.fitsdir = "/group/mwaops/vcs/${params.obsid}/pointings"

params.vcstools_version = 'master'
params.mwa_search_version = 'master'

params.begin = 0
params.end = 0
params.all = false

params.dm_min = 1
params.dm_max = 250

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

//Work out total obs time
if ( params.all ) {
    // an estimation since there's no easy way to make this work
    obs_length = 4805
}
else {
    obs_length = params.end - params.begin + 1
}

if ( "$HOSTNAME".startsWith("galaxy") ) {
    accel_sift_loc = ""
}
else {
    accel_sift_loc = "python /home/nswainst/code/mwa_search/${params.mwa_search_version}"
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
    
    if '$name'.startswith('Blind'):
        output = dd_plan(150., 30.72, 3072, 0.1, $params.dm_min, $params.dm_max)
    else:
        if '$name'.startswith('FRB'):
            dm = fpio.grab_source_alog(source_type='FRB',
                 pulsar_list=['$name'], include_dm=True)[0][-1]
        else:
            # Try RRAT first
            rrat_temp = fpio.grab_source_alog(source_type='RRATs',
                        pulsar_list=['$name'], include_dm=True)
            if len(rrat_temp) == 0:
                #No RRAT so must be pulsar
                dm = fpio.grab_source_alog(source_type='Pulsar',
                     pulsar_list=['$name'], include_dm=True)[0][-1]
            else:
                dm = rrat_temp[0][-1]
        dm_min = float(dm) - 2.0
        if dm_min < 1.0:
            dm_min = 1.0
        dm_max = float(dm) + 2.0
        output = dd_plan(150., 30.72, 3072, 0.1, dm_min, dm_max)
    with open("DDplan.txt", "w") as outfile:
        spamwriter = csv.writer(outfile, delimiter=',')
        for o in output:
            spamwriter.writerow(['$name'] + o)
    """ 
}


process search_dd_fft_acc {
    if ( "$HOSTNAME".startsWith("farnarkle") ) {
        scratch '$JOBFS'
        clusterOptions "--tmp=100GB"
    }
    else {
        container = "nickswainston/presto"
        stageInMode = 'copy'
    }
    label 'cpu'
    time '4h'
    //Will ignore errors for now because I have no idea why it dies sometimes
    errorStrategy 'ignore'

    input:
    tuple val(dm_values), file(fits_files), val(chan)

    output:
    file "*ACCEL_0" optional true
    file "*.inf"
    file "*.singlepulse"
    //Will have to change the ACCEL_0 if I do an accelsearch

    beforeScript "module load singularity/${params.singularity_module}"

    """
    echo "lowdm highdm dmstep ndms timeres downsamp"
    echo ${dm_values}
    file_name=\$(ls *fits | head)
    file_name=\${file_name%%_ch*}
    nsub=\$(calc_nsub.py -f ${(Float.valueOf(chan[0]) + Float.valueOf(chan[-1]))/2*1.28} -dm ${dm_values[1]})
    printf "\\n#Dedispersing the time series at \$(date +"%Y-%m-%d_%H:%m:%S") --------------------------------------------\\n"
    prepsubband -ncpus $task.cpus -lodm ${dm_values[0]} -dmstep ${dm_values[2]} -numdms ${dm_values[3]} -zerodm -nsub \$nsub -numout ${obs_length*10000} -o \${file_name} ${params.obsid}_*.fits
    printf "\\n#Performing the FFTs at \$(date +"%Y-%m-%d_%H:%m:%S") -----------------------------------------------------\\n"
    for i in \$(ls *.dat); do
        realfft \${i}
    done
    printf "\\n#Performing the periodic search at \$(date +"%Y-%m-%d_%H:%m:%S") ------------------------------------------\\n"
    for i in \$(ls *.dat); do
        accelsearch -ncpus $task.cpus -zmax 0 -flo $min_f_harm -fhi $max_f_harm -numharm $params.nharm \${i%.dat}.fft
    done
    single_pulse_search.py -p *.dat
    printf "\\n#Finished at \$(date +"%Y-%m-%d_%H:%m:%S") ----------------------------------------------------------------\\n"
    """
}


process accelsift {
    if ( "$HOSTNAME".startsWith("galaxy") ) {
        container = "nickswainston/presto"
        stageInMode = 'copy'
    }
    label 'cpu'
    time '1h'

    input:
    file accel_inf_single_pulse

    output:
    file "cands_*greped.txt"
    file "${params.obsid}_*_singlepulse.tar.gz"
    file "${params.obsid}_*_singlepulse.ps"

    if ( "$HOSTNAME".startsWith("galaxy") ) {
        beforeScript "module load singularity/${params.singularity_module}"
    }
    else {
        beforeScript "module use ${params.presto_module_dir}; module load presto/${params.presto_module}; module load python/2.7.14; matplotlib/2.2.2-python-2.7.14"
    }

    """
    file_name=\$(ls *singlepulse | head)
    file_name=\${file_name%%_DM*}
    ${accel_sift_loc}ACCEL_sift.py --file_name \${file_name}
    grep \${file_name} cands_\${file_name}.txt > cands_\${file_name}_greped.txt
    single_pulse_search.py *.singlepulse
    tar -czvf singlepulse.tar.gz *DM*.singlepulse
    mv singlepulse.tar.gz \${file_name}_singlepulse.tar.gz
    """
}


process prepfold {
    label 'cpu'
    time '2h'

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
    if [ \$period -gt 0.01 ]; then
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

    prepfold  -o ${cand_line.split()[0]} \
-n \$nbins -noxwin -noclip -p \$period -dm ${cand_line.split()[1]} -nsub 256 -npart \$ntimechunk \
-dmstep \$dmstep -pstep 1 -pdstep 2 -npfact \$period_search_n -ndmfact 1 -runavg ${cand_line.split()[0].split("_DM")[0]}_*.fits
    """
}


process search_dd {
    if ( "$HOSTNAME".startsWith("farnarkle") ) {
        scratch '$JOBFS'
        clusterOptions "--tmp=100GB"
    }
    else {
        container = "nickswainston/presto"
        stageInMode = 'copy'
    }
    label 'cpu'
    time '4h'
    //Will ignore errors for now because I have no idea why it dies sometimes
    errorStrategy 'ignore'

    input:
    tuple val(dm_values), file(fits_files), val(chan)

    output:
    file "*.inf"
    file "*.singlepulse"
    //Will have to change the ACCEL_0 if I do an accelsearch

    beforeScript "module load singularity/${params.singularity_module}"

    """
    echo "lowdm highdm dmstep ndms timeres downsamp"
    echo ${dm_values}
    file_name=\$(ls *fits | head)
    file_name=\${file_name%%_ch*}
    nsub=\$(calc_nsub.py -f ${(Float.valueOf(chan[0]) + Float.valueOf(chan[-1]))/2*1.28} -dm ${dm_values[1]})
    printf "\\n#Dedispersing the time series at \$(date +"%Y-%m-%d_%H:%m:%S") --------------------------------------------\\n"
    prepsubband -ncpus $task.cpus -lodm ${dm_values[0]} -dmstep ${dm_values[2]} -numdms ${dm_values[3]} -zerodm -nsub \$nsub -numout ${obs_length*10000} -o \${file_name} ${params.obsid}_*.fits
    single_pulse_search.py -p *.dat
    """
}


process assemble_single_pulse {
    if ( "$HOSTNAME".startsWith("galaxy") ) {
        container = "nickswainston/presto"
        stageInMode = 'copy'
    }
    label 'cpu'
    time '1h'

    input:
    file accel_inf_single_pulse

    output:
    file "${params.obsid}_*_singlepulse.tar.gz"
    file "${params.obsid}_*_singlepulse.ps"

    if ( "$HOSTNAME".startsWith("galaxy") ) {
        beforeScript "module load singularity/${params.singularity_module}"
    }
    else {
        beforeScript "module use ${params.presto_module_dir}; module load presto/${params.presto_module}; module load python/2.7.14; matplotlib/2.2.2-python-2.7.14"
    }

    """
    file_name=\$(ls *singlepulse | head)
    file_name=\${file_name%%_DM*}
    single_pulse_search.py *.singlepulse
    tar -czvf singlepulse.tar.gz *DM*.singlepulse
    mv singlepulse.tar.gz \${file_name}_singlepulse.tar.gz
    """
}


workflow pulsar_search {
    take:
        name_fits_files
        channels
    main:
        ddplan( name_fits_files )
        search_dd_fft_acc( ddplan.out.splitCsv().map{ it -> [ it[0], [ it[1], it[2], it[3], it[4], it[5], it[6] ] ] }.concat(name_fits_files).groupTuple().\
                           map{ it -> [it[1].init(), [it[1].last()]].combinations() }.flatMap().\
                           combine(channels.map{ it -> [it]}) )
        // Get all the inf, ACCEL and single pulse files and sort them into groups with the same basename (obsid_pointing)
        accelsift( search_dd_fft_acc.out[0].concat(search_dd_fft_acc.out[1], search_dd_fft_acc.out[2]).\
                   flatten().map{ it -> [it.baseName.split("DM")[0], it ] }.groupTuple().map{ it -> it[1] } )
        // Make a pair of accelsift out lines and fits files that match
        prepfold( name_fits_files.map{ it -> it[1] }.flatten().map{ it -> [it.baseName.split("ch")[0], it ] }.groupTuple().cross(\
                  accelsift.out[0].splitText().flatten().map{ it -> [it.split()[0].split("DM")[0], it ] }).\
                  map{ it -> [ it[0][1], it[1][1] ]} )
    emit:
        accelsift.out[1]
        accelsift.out[2]
        prepfold.out
}

workflow single_pulse_search {
    take:
        name_fits_files
        channels
    main:
        ddplan( name_fits_files )
        search_dd( ddplan.out.splitCsv().map{ it -> [ it[0], [ it[1], it[2], it[3], it[4], it[5], it[6] ] ] }.concat(name_fits_files).groupTuple().map{ it -> it[1] }.\
                       combine(channels.map{ it -> [it]}) )
        // Get all the inf and single pulse files and sort them into groups with the same basename (obsid_pointing)
        assemble_single_pulse( search_dd.out[0].concat(search_dd.out[1]).\
                               flatten().map{ it -> [it.baseName.split("DM")[0], it ] }.groupTuple().map{ it -> it[1] } )
    emit:
        assemble_single_pulse.out[0]
        assemble_single_pulse.out[1]
}