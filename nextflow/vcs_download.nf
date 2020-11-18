#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

params.obsid = 'no_obsid'

params.begin = null
params.end = null
params.all = false

params.no_combined_check = true

params.increment = 64
params.parallel_dl = 3
params.untar_jobs = 2
params.keep_tarball = false

params.vcstools_version = 'master'

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
             |  --parallel_dl
             |              Number of parallel downloads to envoke. [default: 3]
             |  --untar_jobs
             |              Number of parallel jobs when untaring downloaded tarballs. [default: 2]
             |  --keep_tarball
             |              Keep the tarballs after unpacking. [default: false]
             |  --vcstools_version
             |              The vcstools module version to use [default: master]
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

    import mwa_metadb_utils as meta
    from mdir import mdir

    data_dir = '${params.scratch_basedir}/${params.obsid}'
    obsinfo = meta.getmeta(service='obs', params={'obs_id':'${params.obsid}'})
    data_format = obsinfo['dataquality']
    if data_format == 1:
        target_dir = link = 'raw'
        data_type = 11
        dl_dir = os.path.join(data_dir, target_dir)
        dir_description = "Raw"
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

    input:
    val data_type
    val dl_dir
    each begin_time_increment

    output:
    val begin_time_increment
    
    beforeScript "module use /group/mwa/software/modulefiles; module load vcstools/${params.vcstools_version}; module load mwa-voltage/master"
    """
    voltdownload.py --obs=$params.obsid --type=$data_type --from=${begin_time_increment[0]} --duration=${begin_time_increment[1] - 1} --parallel=$params.parallel_dl --dir=$dl_dir
    checks.py -m download -o $params.obsid -w $dl_dir -b ${begin_time_increment[0]} -i ${begin_time_increment[1]} --data_type $data_type
    """
}

process untar {
    label 'cpu'
    time { "${50*params.increment*task.attempt + 900}s" }
    errorStrategy 'retry'
    beforeScript "module use ${params.module_dir}; module load vcstools/${params.vcstools_version}"

    input:
    val data_type
    val begin_time_increment

    when:
    data_type == '16'

    """
    untar.sh -w ${params.scratch_basedir}/${params.obsid}/combined -o $params.obsid -j $params.untar_jobs \
-b ${begin_time_increment[0]} -e ${begin_time_increment[0] + begin_time_increment[1] - 1} $keep_tarball_command
    """
}

process recombine {
    label 'gpu'
    time { "${500*params.increment + 900}s" }
    
    if ( "$HOSTNAME".startsWith("garrawarla") ) {
        if ( { params.max_cpus_per_node > begin_time_increment[1] } ) {
            clusterOptions {"--gres=gpu:1 --nodes=${( params.increment - (params.increment % begin_time_increment[1]) ) / begin_time_increment[1] + 1} "+\
                            "--ntasks-per-node=${begin_time_increment[1] / 2}"}
        }
        else {
            clusterOptions {"--gres=gpu:1 --nodes=${1} "+\
                            "--ntasks-per-node=${params.max_cpus_per_node}"}
        }
    }
    else {
        if ( { params.max_cpus_per_node > begin_time_increment[1] } ) {
            clusterOptions {"--nodes=${( params.increment - (params.increment % begin_time_increment[1]) ) / begin_time_increment[1] + 1} "+\
                            "--ntasks-per-node=${begin_time_increment[1] / 2}"}
        }
        else {
            clusterOptions {"--nodes=${1} "+\
                            "--ntasks-per-node=${params.max_cpus_per_node}"}
        }
    }

    beforeScript "module use ${params.module_dir}; module load vcstools/${params.vcstools_version}; module load mwa-voltage/master; module load mpi4py"
    
    input:
    val data_type
    val begin_time_increment

    when:
    data_type == '11'

    script:
    """
    srun --export=all recombine.py -o ${params.obsid} -s ${begin_time_increment[0]} -w ${params.scratch_basedir}/${params.obsid} -e recombine
    checks.py -m recombine -o ${params.obsid} -w ${params.scratch_basedir}/${params.obsid}/combined/ -b ${begin_time_increment[0]} -i ${begin_time_increment[1]}
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
}
