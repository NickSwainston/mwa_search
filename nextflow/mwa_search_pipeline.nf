#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.help = false
if ( params.help ) {
    help = """mwa_search_pipeline.nf: A pipeline that will beamform and perform a pulsar search
             |                        in the entire FOV.
             |Required argurments:
             |  --obsid     Observation ID you want to process [no default]
             |  --calid     Observation ID of calibrator you want to process [no default]
             |  --begin     First GPS time to process [no default]
             |  --end       Last GPS time to process [no default]
             |  --all       Use entire observation span. Use instead of -b & -e. [default: ${params.all}]
             |
             |Pointing arguments (one is required):
             |  --pointings A comma-separated list of pointings with the RA and Dec separated
             |              by _ in the format HH:MM:SS_+DD:MM:SS, e.g.
             |              "19:23:48.53_-20:31:52.95,19:23:40.00_-20:31:50.00" [default: None]
             |  --pointing_file
             |              A file containing pointings with the RA and Dec separated by _
             |              in the format HH:MM:SS_+DD:MM:SS on each line, e.g.
             |              "19:23:48.53_-20:31:52.95\\n19:23:40.00_-20:31:50.00" [default: None]
             |  --bestprof_pointings
             |              A directory that contains bestprof files that you would like to
             |              follow up. The pipeline will beamform on their pointings, prepfold
             |              on their DM and period, and perform a blind search. [default: None]
             |
             |Beamforming types arguments (optional):
             |  --summed   Sum the Stoke paramters [default: ${params.summed}]
             |  --incoh    Also produce an incoherent beam [default: ${params.incoh}]
             |  --ipfb     Also produce a high time resolution Inverse Polyphase Filter Bank beam
             |             [default: ${params.ipfb}]
             |  --offringa Use offringa calibration solution instead of RTS [default: ${params.offringa}]
             |
             |Dedispersion arguments (optional):
             |  --dm_min    Minimum DM to search over [default: ${params.dm_min}]
             |  --dm_max    Maximum DM to search over [default: ${params.dm_max}]
             |  --dm_min_step
             |              Minimum DM step size (Delta DM) [default: ${params.dm_min_step }]
             |  --dm_max_step
             |              Maximum DM step size (Delta DM) [default: ${params.dm_max_step }]
             |  --max_dms_per_job
             |              Maximum number of DM steps a single job will procces.
             |              Lowering this will reduce memory usage and increase parellelisation.
             |              [default: ${params.max_dms_per_job}]
             |
             |Pulsar search arguments (optional):
             |  --sp        Perform only a single pulse search [default: ${params.sp }]
             |  --cand      Label given to output files [default: ${params.cand }]
             |  --nharm     Number of harmonics to search [default: ${params.nharm }]
             |  --min_period
             |              Min period to search for in sec (ANTF min = 0.0013)
             |              [default: ${params.min_period }]
             |  --max_period
             |              Max period to search for in sec (ANTF max = 23.5)
             |              [default: ${params.max_period }]
             |  --zmax      Maximum acceleration to search (0 will do a simpler periodic search).
             |              I recomend you use 200 and set --max_dms_per_job 32
             |              [default: ${params.zmax }]
             |
             |Other arguments (optional):
             |  --out_dir   Output directory for the candidates files
             |              [default: ${params.out_dir}]
             |  --publish_fits
             |              Publish to the fits files to the vcs subdirectory.
             |  --vcstools_version
             |              The vcstools module version to use [default: ${params.vcstools_version}]
             |  --mwa_search_version
             |              The mwa_search module bersion to use [default: ${params.mwa_search_version}]
             |  --no_combined_check
             |              Don't check if all the combined files are available [default: ${params.no_combined_check}]
             |  -w          The Nextflow work directory. Delete the directory once the processs
             |              is finished [default: ${workDir}]""".stripMargin()
    println(help)
    exit(0)
}

// Command line argument parsing

// if ( params.bestprof_pointings ) {
//     bestprof_files = Channel.fromPath("${params.bestprof_pointings}/*.bestprof").collect()
// }
// else {
//     bestprof_files = Channel.empty()
// }
bestprof_files = Channel.empty()
if ( params.pointing_file ) {
    // Grab the pointings from a pointing file
    pointings = Channel
        .fromPath(params.pointing_file)
        .splitCsv()
}
else if ( params.pointings ) {
    // Grab the pointings from the command line
    pointings = Channel
        .of(params.pointings.split(","))
}
else if ( params.bestprof_pointings ) {
    // Grab the pointings from the bestprof files
    bestprof_files = Channel.fromPath("${params.bestprof_pointings}/*.bestprof").collect()
}
else {
    println "No pointings given. Either use --pointing_file or --pointings. Exiting"
    exit(1)
}


process bestprof_pointings {
    input:
    path bestprof_files

    output:
    path "${params.obsid}_DM_pointing.csv"

    """
    #!/usr/bin/env python

    import glob
    import csv
    from vcstools.pointing_utils import format_ra_dec

    dm_pointings = []
    bestprof_files = ["${bestprof_files.join('", "')}"]
    for bfile_loc in bestprof_files:
        pointing = bfile_loc.split("${params.obsid}_")[-1].split("_DM")[0]
        ra, dec = format_ra_dec([[pointing.split("_")[0], pointing.split("_")[1]]])[0]
        with open(bfile_loc,"r") as bestprof:
            lines = bestprof.readlines()
            dm = lines[14][22:-1]
            period = lines[15][22:-1]
            period, period_uncer = period.split('  +/- ')
        dm_pointings.append([f"{ra}_{dec}", f"dm_{dm}_{ra}_{dec}", dm, period])

    with open("${params.obsid}_DM_pointing.csv", "w") as outfile:
        spamwriter = csv.writer(outfile, delimiter=',')
        for dm_point in dm_pointings:
            spamwriter.writerow(dm_point)
    """
}

process follow_up_fold {
    label 'cpu'
    label 'presto'

    time "6h"
    publishDir params.out_dir, mode: 'copy'
    errorStrategy 'retry'
    maxRetries 1

    when:
    params.bestprof_pointings != null

    input:
    tuple path(fits_files), val(dm), val(period)

    output:
    path "*pfd*"

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

    prepfold -ncpus ${task.cpus} -o follow_up_${params.obsid}_P${period.replaceAll(~/\s/,"")}_DM${dm} -n \$nbins -dm ${dm} -p \$period -noxwin -noclip -nsub 256 \
-npart \$ntimechunk -dmstep \$dmstep -pstep 1 -pdstep 2 -npfact \$period_search_n -ndmfact \$ndmfact -runavg *.fits
    """
}

include { pre_beamform; beamform } from './beamform_module'
include { pulsar_search          } from './pulsar_search_module'
include { classifier             } from './classifier_module'

workflow {
    // Parse pointing input
    if ( params.pointing_file ) {
        // Grab the pointings from a pointing file
        pointings = Channel.fromPath(params.pointing_file).splitCsv()
    }
    else if ( params.pointings ) {
        // Grab the pointings from the command line
        pointings = Channel.of(params.pointings.split(","))
    }
    else if ( params.bestprof_pointings ) {
        // Grab the pointings from the bestprof files
        bestprof_files = Channel.fromPath("${params.bestprof_pointings}/*.bestprof").collect()
        // If given bestprof files, extract pointings from them
        bestprof_pointings( bestprof_files )
        bestprof_data = bestprof_pointings.out.splitCsv()
        // Only grab the pointings
        pointings = bestprof_data.map{ point, name, dm, period -> point }
    }
    else {
        println "No pointings given. Either use --pointing_file or --pointings. Exiting"
        exit(1)
    }

    // Beamform
    pre_beamform()
    beamform(
        pre_beamform.out.utc_beg_end_dur,
        pre_beamform.out.channels,
        pointings
    )

    // If doing follow up also do a fold
    if ( params.bestprof_pointings ) {
        // Grab name that includes DM so the ddplan will only use nearby DMs
        bestprof_data_fits = bestprof_data.concat( beamform.out ).groupTuple()
        // Do a follow up fold if bestprofs given
        follow_up_fold(
            bestprof_data_fits.map{ pointing, name_fits, dm, period -> [ name_fits[1], dm[0], period[0] ] }
        )
        name_fits = bestprof_data_fits.map{ pointing, name_fits, dm, period -> name_fits }
    }
    else {
        // Make a name for each fits file
        name_fits = beamform.out.map{ pointing, fits -> [ "${params.cand}_${params.obsid}_${pointing}".toString(), fits ] }
    }

    // Perform pulsar search
    pulsar_search( name_fits )
    classifier( pulsar_search.out )
}
