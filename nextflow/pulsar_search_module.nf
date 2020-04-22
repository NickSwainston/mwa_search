nextflow.preview.dsl = 2

params.obsid = 1253471952
params.fitsdir = "/group/mwaops/vcs/${params.obsid}/pointings"

params.vcstools_version = 'master'
params.mwa_search_version = 'master'

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


process ddplan {
    label 'ddplan'

    output:
    file 'DDplan.txt'
    
    """
    #!/usr/bin/env python3

    from lfDDplan import dd_plan
    import csv
    
    output = dd_plan(150., 30.72, 3072, 0.1, $params.dm_min, $params.dm_max)
    with open("DDplan.txt", "w") as outfile:
        spamwriter = csv.writer(outfile, delimiter=',')
        for o in output:
            spamwriter.writerow(o)
    """ 
}


process search_dd_fft_acc {
    scratch '$JOBFS'
    clusterOptions "--tmp=100GB"
    label 'cpu'
    time '4h'
    //Will ignore errors for now because I have no idea why it dies sometimes
    errorStrategy 'ignore'
    
    input:
    tuple val(lowdm), val(highdm), val(dmstep), val(ndms), val(timeres), val(downsamp)
    each file(fits_files)

    output:
    file "*ACCEL_0" optional true
    file "*.inf"
    file "*.singlepulse"
    //Will have to change the ACCEL_0 if I do an accelsearch

    beforeScript "module use ${params.presto_module_dir}; module load presto/${params.presto_module}"

    """
    file_name=\$(ls *fits | head)
    file_name=\${file_name%%_ch*}
    printf "\\n#Dedispersing the time series at \$(date +"%Y-%m-%d_%H:%m:%S") --------------------------------------------\\n"
    prepsubband -ncpus $task.cpus -lodm $lowdm -dmstep $dmstep -numdms $ndms -zerodm -o \${file_name} ${params.obsid}_*.fits
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
    //container = '/group/mwaops/nswainston/.singularity/cache/oci-tmp/1bb06e133e447f4e70fb325fc6bd7f4dd80987fbfd5fe376dd547592e9b57846/accel_sift_latest.sif'
    //singularity.enabled = true
    label 'cpu'
    time '1h'

    input:
    file accel_inf_single_pulse

    output:
    file "cands_*greped.txt"
    file "${params.obsid}_*_singlepulse.tar.gz"
    file "${params.obsid}_*_singlepulse.ps"

    beforeScript "module use $params.presto_module_dir; module load presto/${params.presto_module}; module load python/2.7.14, matplotlib/2.2.2-python-2.7.14"

    """
    file_name=\$(ls *singlepulse | head)
    file_name=\${file_name%%_DM*}
    python /home/nswainst/code/mwa_search/$params.mwa_search_version/ACCEL_sift.py --file_name \${file_name}
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

workflow pulsar_search {
    take:
        fits_files
    main:
        ddplan()
        search_dd_fft_acc( ddplan.out.splitCsv(),\
                           fits_files )
        // Get all the inf, ACCEL and single pulse files and sort them into groups with the same basename (obsid_pointing)
        accelsift( search_dd_fft_acc.out[0].concat(search_dd_fft_acc.out[1], search_dd_fft_acc.out[2]).\
                   flatten().map{ it -> [it.baseName.split("DM")[0], it ] }.groupTuple().map{ it -> it[1] } )
        // Make a pair of accelsift out lines and fits files that match
        prepfold( fits_files.flatten().map{ it -> [it.baseName.split("ch")[0], it ] }.groupTuple().cross(\
                  accelsift.out[0].splitText().flatten().map{ it -> [it.split()[0].split("DM")[0], it ] }).\
                  map{ it -> [ it[0][1], it[1][1] ]} )
    emit:
        accelsift.out[1]
        accelsift.out[2]
        prepfold.out
}