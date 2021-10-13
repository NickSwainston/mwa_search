#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

params.obsid = null
params.calid = null
params.pointings = null
params.pointing_file = null
params.bestprof_pointings = null

params.begin = 0
params.end = 0
params.all = false

params.summed = true
params.channels = null
params.vcstools_version = 'master'
params.mwa_search_version = 'master'

params.didir = "${params.scratch_basedir}/${params.obsid}/cal/${params.calid}/rts"
params.out_dir = "${params.search_dir}/${params.obsid}_candidates"

params.dm_min = 1
params.dm_max = 250
params.dm_min_step = 0.02
params.max_dms_per_job = 5000
params.zmax = 0

params.no_combined_check = false


if ( params.max_dms_per_job != 5000 ) {
    // If using non default max_dms_per_job then use a make the groupTuple size sudo infinite
    total_dm_jobs = 10000
}
// If doing an acceleration search, lower the number of DMs per job so the jobs don't time out
else if ( params.zmax == 0 ) {
    // Periodic search defaults
    total_dm_jobs = 6
}
else {
    // Accel search defaults
    total_dm_jobs = 24
    params.max_dms_per_job = 128
}

if ( params.bestprof_pointings ) {
    bestprof_files = Channel.fromPath("${params.bestprof_pointings}/*.bestprof").collect()
}
else {
    bestprof_files = Channel.from(" ")
}

params.help = false
if ( params.help ) {
    help = """mwa_search_pipeline.nf: A pipeline that will beamform and perform a pulsar search
             |                        in the entire FOV.
             |Required argurments:
             |  --obsid     Observation ID you want to process [no default]
             |  --calid     Observation ID of calibrator you want to process [no default]
             |  --pointings A comma sepertated list of pointings with the RA and Dec seperated
             |              by _ in the format HH:MM:SS_+DD:MM:SS, e.g.
             |              "19:23:48.53_-20:31:52.95,19:23:40.00_-20:31:50.00" [default: None]
             |  --pointing_file
             |              A file containing pointings with the RA and Dec seperated by _
             |              in the format HH:MM:SS_+DD:MM:SS on each line, e.g.
             |              "19:23:48.53_-20:31:52.95\\n19:23:40.00_-20:31:50.00" [default: None]
             |  --begin     First GPS time to process [no default]
             |  --end       Last GPS time to process [no default]
             |  --all       Use entire observation span. Use instead of -b & -e. [default: false]
             |  --summed    Add this flag if you the beamformer output as summed polarisations
             |              (only Stokes I). This reduces the data size by a factor of 4.
             |              [default: False]
             |
             |Optional arguments:
             |  --dm_min    Minimum DM to search over [default: 1]
             |  --dm_max    Maximum DM to search over [default: 250]
             |  --dm_min_step
             |              Minimum DM step size (Delta DM) [default: 0.1]
             |  --out_dir   Output directory for the candidates files
             |              [default: ${params.search_dir}/<obsid>_candidates]
             |  --ipfb      Perform an the inverse PFB to produce high time resolution beamformed
             |              vdif files [default: false]
             |  --publish_fits
             |              Publish to the fits directory (/group on Galaxy).
             |  --publish_fits_scratch
             |              Publish to the scratch fits directory (/astro on Galaxy).
             |  --vcstools_version
             |              The vcstools module version to use [default: master]
             |  --mwa_search_version
             |              The mwa_search module bersion to use [default: master]
             |  --no_combined_check
             |              Don't check if all the combined files are available [default: false]
             |  -w          The Nextflow work directory. Delete the directory once the processs
             |              is finished [default: ${workDir}]""".stripMargin()
    println(help)
    exit(0)
}

process bestprof_pointings {
    input:
    val pointings
    file bestprof_files

    output:
    file "${params.obsid}_DM_pointing.csv"

    """
    #!/usr/bin/env python

    import glob
    import csv
    from vcstools.pointing_utils import format_ra_dec

    dm_pointings = []
    if "${params.bestprof_pointings}" == "null":
        pointings = ["${pointings.join('", "')}"]
        for p in pointings:
            ra, dec = format_ra_dec([[p.split("_")[0], p.split("_")[1]]])[0]
            dm_pointings.append(["{}_{}".format(ra, dec), "Blind", "None"])
    else:
        bestprof_files = glob.glob("*.bestprof")
        for bfile_loc in bestprof_files:
            pointing = bfile_loc.split("${params.obsid}_")[-1].split("_DM")[0]
            ra, dec = format_ra_dec([[pointing.split("_")[0], pointing.split("_")[1]]])[0]
            with open(bfile_loc,"r") as bestprof:
                lines = bestprof.readlines()
                dm = lines[14][22:-1]
                period = lines[15][22:-1]
                period, period_uncer = period.split('  +/- ')
            dm_pointings.append(["{}_{}".format(ra, dec), "dm_{}".format(dm), period])

    with open("${params.obsid}_DM_pointing.csv", "w") as outfile:
        spamwriter = csv.writer(outfile, delimiter=',')
        for dm_point in dm_pointings:
            spamwriter.writerow(dm_point)
    """
}

process follow_up_fold {
    label 'cpu'
    time "6h"
    publishDir params.out_dir, mode: 'copy'
    errorStrategy 'retry'
    maxRetries 1

    when:
    params.bestprof_pointings != null

    input:
    tuple file(fits_files), val(dm), val(period)

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
    """
    # Set up the prepfold options to match the ML candidate profiler
    temp_period=${Float.valueOf(period)/1000}
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

    prepfold -ncpus $task.cpus -o follow_up_${params.obsid}_P${period.replaceAll(~/\s/,"")}_DM${dm} -n \$nbins -dm ${dm} -p \$period -noxwin -noclip -nsub 256 \
-npart \$ntimechunk -dmstep \$dmstep -pstep 1 -pdstep 2 -npfact \$period_search_n -ndmfact \$ndmfact -runavg *.fits
    """
}

if ( params.pointing_file ) {
    pointings = Channel
        .fromPath(params.pointing_file)
        .splitCsv()
        .collect()
}
else if ( params.pointings ) {
    pointings = Channel
        .from(params.pointings.split(","))
        .collect()
}
else if ( params.bestprof_pointings ) {
    pointings = Channel.from("null")
}
else {
    println "No pointings given. Either use --pointing_file or --pointings. Exiting"
    exit(1)
}

include { pre_beamform; beamform } from './beamform_module'
include { pulsar_search } from './pulsar_search_module'
include { classifier }   from './classifier_module'

workflow {
    bestprof_pointings( pointings,
                        bestprof_files )
    pre_beamform()
    beamform( pre_beamform.out[0],\
              pre_beamform.out[1],\
              pre_beamform.out[2],\
              bestprof_pointings.out.splitCsv().map{ it -> it[0] }.flatten().unique().collate( params.max_pointings ) )
    follow_up_fold( beamform.out[1].map{ it -> [ it[0].getBaseName().split("/")[-1].split("_ch")[0], it ] }.cross(
                    bestprof_pointings.out.splitCsv().map{ it -> ["${params.obsid}_"+it[0], it[1], it[2]]}.map{ it -> [it[0].toString(), [it[1], it[2]]] }).\
                    map{ it -> [it[0][1], it[1][1][0].split("_")[-1], it[1][1][1]] } )
    pulsar_search( beamform.out[1].map{ it -> [ it[0].getBaseName().split("/")[-1].split("_ch")[0], it ] }.concat(
                   bestprof_pointings.out.splitCsv().map{ it -> ["${params.obsid}_"+it[0], it[1]]}).map{ it -> [it[0].toString(), it[1]] }.\
                   groupTuple( size: 2 ).map{ it -> [it[1][1]+"_"+it[0], it[1][0]] } )
    classifier( pulsar_search.out[1].flatten().collate( 120 ) )
}
