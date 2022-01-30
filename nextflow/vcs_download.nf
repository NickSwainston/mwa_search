#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.obsid = 'no_obsid'

params.begin = null
params.end = null
params.all = false

params.no_combined_check = true
params.ozstar_transfer = false

params.increment = 32
params.parallel_dl = 3
params.untar_jobs = 2
params.keep_tarball = false
params.keep_raw = false
params.max_jobs = 12

params.vcstools_version = 'master'
params.mwa_voltage_version = 'master'

if ( params.keep_tarball ) {
    keep_tarball_command = "-k"
}
else {
    keep_tarball_command = ""
}

params.help = false
if ( params.help ) {
    help = """vcs_download.nf: A pipeline that will download vcs data and untar or recombine it if required.
             |Required argurments:
             |  --obsid     Observation ID you want to process [no default]
             |  --begin     First GPS time to process [no default]
             |  --end       Last GPS time to process [no default]
             |  --all       Use entire observation span. Use instead of -b & -e. [default: false]
             |
             |Optional arguments:
             |  --increment Increment in seconds (how much we process at once). [default: 64]
             |  --max_jobs  Number of maximum jobs of each type to run at once (to limit IO). [default: 12]
             |  --parallel_dl
             |              Number of parallel downloads to envoke. [default: 3]
             |  --untar_jobs
             |              Number of parallel jobs when untaring downloaded tarballs. [default: 2]
             |  --keep_tarball
             |              Keep the tarballs after unpacking. [default: false]
             |  --keep_raw
             |              Keep the raw data after recombining. [default: false]
             |  --vcstools_version
             |              The vcstools module version to use [default: master]
             | --mwa_voltage_version
             |              The mwa-voltage module version to use [default: master]
             |  -w          The Nextflow work directory. Delete the directory once the processs
             |              is finished [default: ${workDir}]""".stripMargin()
    println(help)
    exit(0)
}

// Argument checking
if ( params.all && ( params.begin != null || params.begin != null ) ) {
    println("Please use either --all or --begin and --end. Not both!")
    exit(0)
}
else if ( params.all == false && ( params.begin == null || params.begin == null ) ) {
    println("Please use either --all or --begin and --end.")
    exit(0)
}
if ( params.obsid == 'no_obsid') {
    println("Please use --obsid.")
    exit(0)
}
if ( params.increment > params.max_cpus_per_node ) {
    println("Please us an --increment less than the max number of cpus per node which is ${params.max_cpus_per_node}.")
    exit(0)
}

process check_data_format {
    input:
    tuple val(begin), val(end)

    output:
    file "${params.obsid}_data_type_dir.txt"
    file "${params.obsid}_time_increments.txt"

    """
    #!/usr/bin/env python

    import os
    import csv

    import vcstools.metadb_utils as meta
    from vcstools.general_utils import mdir

    # Ensure the metafits files is there
    meta.ensure_metafits("${params.basedir}/${params.obsid}", "${params.obsid}",\
                         "${params.scratch_basedir}/${params.obsid}/${params.obsid}_metafits_ppds.fits")

    data_dir = '${params.scratch_basedir}/${params.obsid}'
    obsinfo = meta.getmeta(service='obs', params={'obs_id':'${params.obsid}'})
    comb_del_check = meta.combined_deleted_check(${params.obsid}, begin=${begin}, end=${end})
    data_format = obsinfo['dataquality']
    if data_format == 1 or (comb_del_check and data_format == 6):
        # either only the raw data is available (data_format == 1)
        # or there was combined files but they were deleted (comb_del_check and data_format == 6)
        target_dir = link = 'raw'
        data_type = 11
        dl_dir = os.path.join(data_dir, target_dir)
        dir_description = "Raw"
        # also make combined dir
        mdir(os.path.join(data_dir, 'combined'), "Combined")
    elif data_format == 6:
        target_dir = link = 'combined'
        data_type = 16
        dl_dir = os.path.join(data_dir, target_dir)
        dir_description = "Combined"
    else:
        logger.error("Unable to determine data format from archive. Exiting")
        sys.exit(1)
    mdir(dl_dir, dir_description)

    with open("${params.obsid}_data_type_dir.txt", "w") as outfile:
        spamwriter = csv.writer(outfile, delimiter=',')
        spamwriter.writerow([data_type])
        spamwriter.writerow([dl_dir])

    # Work out time increments
    with open("${params.obsid}_time_increments.txt", "w") as outfile:
        spamwriter = csv.writer(outfile, delimiter=',')
        for time_to_get in range(${begin}, ${end}, ${params.increment}):
            if time_to_get + ${params.increment} > ${end}:
                increment = ${end} - time_to_get + 1
            else:
                increment = ${params.increment}
            spamwriter.writerow([time_to_get, increment])
    """
}


process volt_download {
    label 'download'
    time { "${500*params.increment*task.attempt + 900}s" }
    errorStrategy { task.attempt > 3 ? 'finish' : 'retry' }
    maxRetries 3
    maxForks params.max_jobs

    input:
    val data_type
    val dl_dir
    each begin_time_increment

    output:
    val begin_time_increment

    beforeScript "module use /group/mwa/software/modulefiles; module load vcstools/${params.vcstools_version}; module load mwa-voltage/${params.mwa_voltage_version}"
    """
    voltdownload.py --obs=$params.obsid --type=$data_type --from=${begin_time_increment[0]} --duration=${begin_time_increment[1] - 1} --parallel=$params.parallel_dl --dir=$dl_dir
    checks.py -m download -o $params.obsid -w $dl_dir -b ${begin_time_increment[0]} -i ${begin_time_increment[1]} --data_type $data_type
    """
}

process untar {
    label 'cpu'
    time { "${200*params.increment*task.attempt + 900}s" }
    errorStrategy { task.attempt > 3 ? 'finish' : 'retry' }
    maxRetries 3
    maxForks params.max_jobs

    beforeScript "module use ${params.module_dir}; module load vcstools/${params.vcstools_version}"

    input:
    val data_type
    val begin_time_increment

    output:
    val begin_time_increment

    when:
    data_type == '16'

    """
    untar.sh -w ${params.scratch_basedir}/${params.obsid}/combined -o $params.obsid -j $params.untar_jobs \
-b ${begin_time_increment[0]} -e ${begin_time_increment[0] + begin_time_increment[1] - 1} $keep_tarball_command
    """
}

process recombine {
    label 'cpu'
    time { "${500*params.increment*task.attempt + 900}s" }
    errorStrategy { task.attempt > 3 ? 'finish' : 'retry' }
    maxRetries 3
    maxForks params.max_jobs
    clusterOptions {"--nodes=1 --ntasks-per-node=${begin_time_increment[1]}"}

    beforeScript "module use ${params.module_dir}; module load vcstools/${params.vcstools_version}; module load mwa-voltage/${params.mwa_voltage_version}; module load gcc/8.3.0; module load cfitsio; module load mpi4py"

    input:
    val data_type
    val begin_time_increment

    output:
    val begin_time_increment

    when:
    data_type == '11'

    script:
    """
    srun --export=all recombine.py -o ${params.obsid} -s ${begin_time_increment[0]} -w ${params.scratch_basedir}/${params.obsid} -e recombine
    checks.py -m recombine -o ${params.obsid} -w ${params.scratch_basedir}/${params.obsid}/combined/ -b ${begin_time_increment[0]} -i ${begin_time_increment[1]}
    if ! ${params.keep_raw}; then
        # Loop over each second and delete raw files
        for gps in \$(seq ${begin_time_increment[0]} ${begin_time_increment[0] + begin_time_increment[1] - 1}); do
            rm ${params.scratch_basedir}/${params.obsid}/raw/${params.obsid}_\${gps}_vcs*.dat
        done
    fi
    """
}

process ozstar_transfer {
    label 'download'
    time { "${500*params.increment*task.attempt + 900}s" }
    errorStrategy { task.attempt > 3 ? 'finish' : 'retry' }
    maxRetries 3
    maxForks 3

    input:
    val begin_time_increment

    """
    start=${begin_time_increment[0]}
    end=${begin_time_increment[0] + begin_time_increment[1] - 1}
    echo "obsid: ${params.obsid} start: \${start} end: \${end}"

    ls --format single-column /astro/mwavcs/vcs/${params.obsid}/combined/*{${begin_time_increment[0]}..${begin_time_increment[0] + begin_time_increment[1] - 1}}*dat | xargs -n1 basename > temp_file_list.txt
    rsync -vhu --files-from=temp_file_list.txt /astro/mwavcs/vcs/${params.obsid}/combined/ ozstar:/fred/oz125/vcs/1221399680/combined
    rm /astro/mwavcs/vcs/${params.obsid}/combined/*{${begin_time_increment[0]}..${begin_time_increment[0] + begin_time_increment[1] - 1}}*dat
    """
}

include { pre_beamform } from './beamform_module'

workflow {
    pre_beamform()
    check_data_format( pre_beamform.out[0] )
    volt_download( check_data_format.out[0].splitCsv().collect().map{ it -> it[0] },
                   check_data_format.out[0].splitCsv().collect().map{ it -> it[1] },
                   check_data_format.out[1].splitCsv().map{ it -> [Integer.valueOf(it[0]), Integer.valueOf(it[1])] } )
    untar( check_data_format.out[0].splitCsv().collect().map{ it -> it[0] },
           volt_download.out )
    recombine( check_data_format.out[0].splitCsv().collect().map{ it -> it[0] },
               volt_download.out )
    if ( params.ozstar_transfer ) {
        ozstar_transfer( untar.out.concat( recombine.out ) )
    }
}
