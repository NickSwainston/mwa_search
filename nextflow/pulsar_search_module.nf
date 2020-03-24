nextflow.preview.dsl = 2

params.obsid = 1253471952
params.fitsdir = "/group/mwaops/vcs/${params.obsid}/pointings"
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
    //scratch true
    label 'cpu'
    
    input:
    tuple val(lowdm), val(highdm), val(dmstep), val(ndms), val(timeres), val(downsamp)
    file fits_files

    output:
    file "*ACCEL*"
    //file "*.singlepulse"

    beforeScript 'module use /group/mwa/software/modulefiles; module load mwa_search/master'

    """
    printf "\\n#Dedispersing the time series at \$(date +"%Y-%m-%d_%H:%m:%S") --------------------------------------------\\n"
    prepsubband -ncpus $task.cpus -lodm $lowdm -dmstep $dmstep -numdms $ndms -zerodm -o ${params.obsid} ${params.obsid}_*.fits
    printf "\\n#Performing the FFTs at \$(date +"%Y-%m-%d_%H:%m:%S") -----------------------------------------------------\\n"
    for i in \$(ls *.dat); do
        realfft \${i}
    done
    printf "\\n#Performing the periodic search at \$(date +"%Y-%m-%d_%H:%m:%S") ------------------------------------------\\n"
    for i in \$(ls *.dat); do
        accelsearch -ncpus $task.cpus -zmax 0 -flo $min_f_harm -fhi $max_f_harm -numharm $params.nharm \${i%.dat}.fft
    done
    #single_pulse_search.py -p *.dat
    printf "\\n#Finished at \$(date +"%Y-%m-%d_%H:%m:%S") ----------------------------------------------------------/------\\n"
    """
}


process accelsift {
    //container = '/group/mwaops/nswainston/.singularity/cache/oci-tmp/1bb06e133e447f4e70fb325fc6bd7f4dd80987fbfd5fe376dd547592e9b57846/accel_sift_latest.sif'
    //singularity.enabled = true

    input:
    //file single_pulse
    file accel

    output:
    file "cands*.txt" //This is probably wrong
    //file "${params.obsid}_${params.pointing}_singlepulse.tar.gz"
    //file "${params.obsid}_${params.pointing}_singlepulse.ps"

    //removed the parts that don't work which is most of this
    """
    #ACCEL_sift.py ./
    cp /group/mwaops/nswainston/pulsar_search/cand_files/test_cand.txt cands_00:36:10.34_-10:33:25.93_test.txt
    #single_pulse_search.py *.singlepulse
    #mv ${params.obsid}_singlepulse.ps ${params.obsid}_${params.pointing}_singlepulse.ps
    #tar -czvf singlepulse.tar.gz *DM*.singlepulse
    #mv singlepulse.tar.gz ${params.obsid}_${params.pointing}_singlepulse.tar.gz
    """
}


process prepfold {
    input:
    each cand_line
    //tuple val(accel_file_name), val(cand_num), val(cand_DM), val(period)
    val cand_file_name
    file fits_files

    output:
    file "*pfd*"

    //no mask command currently
    """
    echo ${cand_line.split()}
    # Set up the prepfold options to match the ML candidate profiler
    period=${Float.valueOf(cand_line.split()[3])/1000}
    if [ \$period -gt 10 ]; then
        nbins=100
        ntimechunk=120
        dmstep=1
        period_search_n=1
    else; then
        # bin size is smaller than time resolution so reduce nbins
        nbins=50
        ntimechunk=40
        dmstep=3
        period_search_n=2
    fi

    prepfold  -o ${cand_line.split()[0].replace(':', '_')}_${cand_file_name.split("_")[1]}_${cand_file_name.split("_")[1]} \
-n \$nbins -noxwin -noclip -p \$period -dm ${cand_line.split()[2]} -nsub 256 -npart \$ntimechunk \
-dmstep \$dmstep -pstep 1 -pdstep 2 -npfact \$period_search_n -ndmfact 1 -runavg ${params.obsid}_*.fits
    """
}

workflow pulsar_search {
    take:
        fits_files
        pointing
    main:
        ddplan()
        search_dd_fft_acc( ddplan.out.splitCsv(),\
                           fits_files )
        accelsift( search_dd_fft_acc.out[0].collect())//,\
                   //search_dd_fft_acc.out[1] )
        prepfold( accelsift.out[0].splitText(),\
                  accelsift.out[0],\
                  fits_files )
    emit:
        //accelsift.out[1]
        //accelsift.out[2]
        prepfold.out
}

