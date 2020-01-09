params.obsid = 1253471952
params.pointing = '02:55:56.29_-53:04:21.27'
params.fitsdir = "/group/mwaops/vcs/${params.obsid}/pointings/${params.pointing}"
params.dm_min = 1
params.dm_max = 250


process ddplan {
    label 'ddplan'

    output:
    file 'DDplan.txt' into ddplan
    
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


ddplan
    .splitCsv()
    .into { ddlist_ch; ddlist_calc_ch }
/*
process calc_dm_range {
    input:
    tuple val(lowdm), val(highdm), val(dmstep), val(ndms), val(timeres), val(downsamp) from ddlist_calc_ch

    //output:
    //file 'DM_list.csv' into dm_csv
    """
    echo
    """
    
}*/



process dedisperse {
    //publishDir 'results', mode: 'copy', saveAs: { filename -> "out_${params.obsid}.out"     }

    tag "${lowdm}"
    label 'dedisperse'
    
    input:
    tuple val(lowdm), val(highdm), val(dmstep), val(ndms), val(timeres), val(downsamp) from ddlist_ch

    output: 
    file "*_DM*" into dm_files
    file "*DM*.inf" into dm_infs

    //TODO add mask command
    //TODO add suband calc -nsub {2:d}
    //TODO add -numout {5}
    """
    echo prepsubband -ncpus $task.cpus -lodm $lowdm -dmstep $dmstep -numdms $ndms -zerodm -o $params.obsid $params.fitsdir/1*.fits
    prepsubband -ncpus $task.cpus -lodm $lowdm -dmstep $dmstep -numdms $ndms -zerodm -o $params.obsid $params.fitsdir/${params.obsid}_0001*.fits
    """
}

dm_files
    .flatten()
    .map { it -> [it.baseName, it ]}
    .groupTuple(size: 2, sort: true)
    .map {it -> it[1]}
    .collate(200)
    .map { it -> it.flatten()}
    .set{ dm_files }

//Some math for the accelsearch command
nharm = 16 // number of harmonics to search
min_period = 0.001 // min period to search for in sec (ANTF min = 0.0013)
max_period = 30 // max period to search for in sec  (ANTF max = 23.5)
//convert to freq
min_freq = 1 / max_period
max_freq = 1 / min_period
//adjust the freq to include the harmonics
min_f_harm = min_freq
max_f_harm = max_freq * nharm

process fft_accelsearch {
    label 'fft'

    input:
    file dm_file from dm_files

    output:
    //file "*_DM*.fft" into fft_files
    file "*ACCEL*" into accel_files

    """
    for i in \$(ls *.dat); do
        realfft \${i}
        accelsearch -ncpus $task.cpus -zmax 0 -flo $min_f_harm -fhi $max_f_harm -numharm $nharm \${i%.dat}.fft
    done
    """
}


/*
process accelsearch {
    label 'accelsearch'
    
    input:
    tuple val(basename), file(fft_file) from fft_paired

    output:
    file "*ACCEL*" into accel_files
    
    //TODO add zmax option
    """
    echo ${fft_file[1]}
    accelsearch -ncpus $task.cpus -zmax 0 -flo $min_f_harm -fhi $max_f_harm -numharm $nharm ${basename}.fft
    """
}*/


process accelsift{
    container = '/group/mwaops/nswainston/.singularity/cache/oci-tmp/1bb06e133e447f4e70fb325fc6bd7f4dd80987fbfd5fe376dd547592e9b57846/accel_sift_latest.sif'
    //singularity.enabled = true

    input:
    file accel_file from accel_files.collect()

    """
    ACCEL_sift.py ./
    """
}

