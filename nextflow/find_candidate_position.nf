#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

params.obsid = null
params.calid = null
params.pointings = null
params.pointing_file = null
params.begin = 0
params.end = 0
params.all = false

params.pointing_grid = null
params.fwhm_deg = null
params.fraction = 0.8
params.loops = 1

params.summed = true
params.channels = null
params.vcstools_version = 'master'
params.mwa_search_version = 'master'

params.bins = 128
params.subint = 60
params.nchan = 48
params.pulsar = 0
params.period = 0.90004
params.dm = 23.123

params.no_pdmp = false
params.fwhm_ra = "None"
params.fwhm_dec = "None"

if ( params.no_pdmp ) {
    tuple_size = 4
}
else {
    tuple_size = 7
}

params.help = false
if ( params.help ) {
    help = """mwa_search_pipeline.nf: A pipeline that will beamform and perform a pulsar search
             |                        in the entire FOV.
             |Required argurments:
             |  --obsid     Observation ID you want to process [no default]
             |  --calid     Observation ID of calibrator you want to process [no default]
             |  --begin     First GPS time to process [no default]
             |  --end       Last GPS time to process [no default]
             |  --all       Use entire observation span. Use instead of -b & -e. [default: false]
             |  --summed    Add this flag if you the beamformer output as summed polarisations
             |              (only Stokes I). This reduces the data size by a factor of 4.
             |              [default: False]
             |
             |Required pointing arguments:
             |  --pointings A space sepertated list of pointings with the RA and Dec seperated
             |              by _ in the format HH:MM:SS_+DD:MM:SS, e.g.
             |              "19:23:48.53_-20:31:52.95 19:23:40.00_-20:31:50.00" [default: None]
             |  --pointing_file
             |              A file containing pointings with the RA and Dec seperated by _
             |              in the format HH:MM:SS_+DD:MM:SS on each line, e.g.
             |              "19:23:48.53_-20:31:52.95\\n19:23:40.00_-20:31:50.00" [default: None]
             |
             |Pointing grid arguments:
             | --pointing_grid
             |              Pointing which grid.py will make a loop of pointings around eg.
             |              "19:23:48.53_-20:31:52.95" [default: None]
             | --fwhm_deg   The FWHM of the observation in degrees (used by grid.py) [default: 0.021]
             | --fraction   The fraction of the FWHM to space the grid by [default: 0.8]
             | --loops      The number of loops of beamd to surround the centre pointing [default: 1]
             |
             |Presto and dspsr options:
             | --bins       Number of bins to use [default: 128]
             | --period     Period in seconds to fold on [default: 0.90004]
             | --dm         The dispersion measure to use [default: 23.123]
             | --subint     The number of subints to use in pmdp [default: 60]
             | --nchan      The number of subchans to use in pmdp [default: 48]
             |
             |Optional arguments:
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

if ( params.pulsar == 0 ) {
    eph_command = "-p ${params.period} -dm ${params.dm}"
    psrcat_command = ""
}
else {
    eph_command = "-par ${params.pulsar}.eph"
    psrcat_command = "psrcat -e ${params.pulsar} | grep -v UNITS > ${params.pulsar}.eph"
}


include { pre_beamform; beamform } from './beamform_module'
include { fwhm_calc } from './data_processing_pipeline'

params.didir = "${params.scratch_basedir}/${params.obsid}/cal/${params.calid}/rts"
params.out_dir = "${params.search_dir}/${params.obsid}_candidate_follow_up"
params.final_dir = "${params.search_dir}/psr2_J0024-1932/${params.obsid}"

if ( params.pointing_file ) {
    pointings = Channel
        .fromPath(params.pointing_file)
        .splitCsv()
        .collect()
        .flatten()
        .collate( params.max_pointings )
}
else if ( params.pointings ) {
    pointings = Channel
        .from(params.pointings.split(","))
        .collect()
        .flatten()
        .collate( params.max_pointings )
}
else if ( params.pointing_grid ) {
    pointing_grid = Channel.from(params.pointing_grid)
}
else if ( ! params.pulsar ) {
    println "No pointings given. Either use --pointing_file, --pointings or --pointing_grid. Exiting"
    exit(1)
}

if ( params.no_pdmp ) {
    input_sn_option = " -b "
}
else {
    input_sn_option = " -p "
}


process get_pulsar_ra_dec {
    output:
    file 'pulsar_ra_dec.txt'

    """
    #!/usr/bin/env python3

    import csv
    from vcstools.catalogue_utils import get_psrcat_ra_dec
    from vcstools.pointing_utils import format_ra_dec

    pulsar_list = ["${params.pulsar.split(",").join('","')}"]
    pulsar_ra_dec = get_psrcat_ra_dec(pulsar_list=pulsar_list)
    pulsar_ra_dec = format_ra_dec(pulsar_ra_dec, ra_col = 1, dec_col = 2)
    pointing = []
    for prd in pulsar_ra_dec:
        pointing.append("{}_{}".format(prd[1], prd[2]))

    with open("pulsar_ra_dec.txt", "w") as outfile:
        spamwriter = csv.writer(outfile, delimiter=',')
        for prd in pulsar_ra_dec:
            spamwriter.writerow([prd[0], "{}_{}".format(prd[1], prd[2])])
    """
}

process grid {
    input:
    tuple val(pulsar), val(pointings), val(fwhm)

    output:
    file "*txt"

    """
    grid.py -o $params.obsid -d $fwhm -f $params.fraction -p $pointings -l $params.loops --label ${pulsar}
    """

}

process prepfold {
    label 'cpu'
    time '3h'
    publishDir params.out_dir, mode: 'copy'

    input:
    tuple val(pointing), file(fits_files), val(pulsar)

    output:
    file "*bestprof"
    file "*ps"

    if ( "$HOSTNAME".startsWith("farnarkle") ) {
        beforeScript "module use ${params.presto_module_dir}; module load presto/${params.presto_module}"
    }
    else if ( "$HOSTNAME".startsWith("x86") || "$HOSTNAME".startsWith("garrawarla") || "$HOSTNAME".startsWith("galaxy") ) {
        container = "file:///${params.containerDir}/presto/presto.sif"
    }
    else {
        container = "nickswainston/presto:realfft_docker"
    }

    //no mask command currently
    """
    if [ ${params.pulsar} == 0 ]; then
        eph_command="-p ${params.period} -dm ${params.dm}"
    else
        eph_command="-par ${pulsar}.eph"
        psrcat -e ${pulsar} | grep -v UNITS > ${pulsar}.eph
    fi
    prepfold -ncpus ${task.cpus} -o ${params.obsid}_${pointing}_pos -n ${params.bins} \${eph_command} -noxwin -noclip -nsub 256 -npart 120 \
-dmstep 1 -pstep 1 -pdstep 2 -npfact 1 -ndmfact 1 -runavg ${params.obsid}*.fits
    """
}

process pdmp {
    label 'cpu'
    time '8h'
    publishDir params.out_dir, mode: 'copy'

    when:
    params.no_pdmp == false

    input:
    tuple val(pointings), file(bestprof), file(fits_files)

    output:
    file "*ps"
    file "*posn"
    file "*ar"

    if ( "$HOSTNAME".startsWith("farnarkle") ) {
        beforeScript "module use ${params.presto_module_dir}; module load dspsr/master"
    }
    else if ( "$HOSTNAME".startsWith("x86") || "$HOSTNAME".startsWith("garrawarla") || "$HOSTNAME".startsWith("galaxy") ) {
        container = "file:///${params.containerDir}/dspsr/dspsr.sif"
    }
    else {
        container = "nickswainston/dspsr_docker"
    }

    //may need to add some channel names
    """
    DM=\$(grep DM *.bestprof | tr -s ' ' | cut -d ' ' -f 5)
    echo "DM: \$DM"
    period=\$(grep P_topo *.bestprof | tr -s ' ' | cut -d ' ' -f 5)
    period="\$(echo "scale=10;\${period}/1000"  |bc)"
    echo "period: \$period"
    samples="\$(grep "Data Folded" *.bestprof | tr -s ' ' | cut -d ' ' -f 5)"
    #One subint per 30 seconds
    subint=\$(python -c "print('{:d}'.format(int(\$samples/300000)))")
    if [ \$subint -lt 30 ]; then subint=30; fi
    echo "subint: \$subint"
    name=${params.obsid}_${pointings}_pos_${bestprof.baseName.split("pos_")[1].split(".pfd")[0]}
    dspsr -t $task.cpus -b ${params.bins} -c \${period} -D \${DM} -L \${subint} -e subint -cont -U 4000 ${params.obsid}*.fits
    psradd *.subint -o \${name}.ar
    pam --setnchn ${params.nchan} -m \${name}.ar
    pdmp -g \${name}_pdmp.ps/cps \${name}.ar
    mv pdmp.posn \${name}_pdmp.posn
    """
}

process bestgridpos {
    publishDir params.out_dir, mode: 'copy'

    input:
    tuple val(pulsar), file(posn_or_bestprof), val(fwhm), val(orig_pointing)

    output:
    file "*predicted_pos.txt"
    file "*png"
    file "*orig_best_SN.txt"

    """
    if [[ ${params.fwhm_ra} == None || ${params.fwhm_dec} == None ]]; then
        fwhm_option="-fr ${fwhm} -fd ${fwhm}"
    else
        fwhm_option="-fr ${params.fwhm_ra} -fd ${params.fwhm_dec}"
    fi
    bestgridpos.py -o ${params.obsid} -O ${params.calid} ${input_sn_option} ./ -w \$fwhm_option --orig_pointing ${[orig_pointing].flatten().findAll{ it != null }.join(" ")} --label ${pulsar}
    """
}

process format_output {
    publishDir params.final_dir, mode: 'copy'
    echo true

    input:
    tuple file(orig_best_file), file(posn_or_bestprof)

    output:
    file "*orig_best_predicted_sn.csv"

    """
    #!/usr/bin/env python3

    import csv

    # process input bestprof or posn files
    if "${params.no_pdmp}" == "true":
        with open("${posn_or_bestprof.baseName}.bestprof" , "r") as bestprof:
            lines = bestprof.readlines()
            ra, dec = lines[0].split("=")[-1].split("_")[1:3]
            sn = float(lines[13].split("~")[-1].split(" ")[0])
    if "${params.no_pdmp}" == "false":
        with open("${posn_or_bestprof.baseName}.posn" , "r") as pdmp:
            lines = pdmp.readlines()
            sn = float(lines[0].split()[3])
            ra, dec = lines[0].split()[9].split("_")[1:3]
            dec = dec

    # read input csv
    with open("${orig_best_file}", "r") as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',')
        csv_output = []
        for row in spamreader:
            csv_output.append(row)
        orig_point, orig_sn = csv_output[0]
        best_point, best_sn = csv_output[1]

    # output finial file
    with open("{}_{}_${orig_best_file.baseName.split("_")[0]}_orig_best_predicted_sn.csv".format("${params.obsid}", "${params.calid}"), "w") as outfile:
        spamwriter = csv.writer(outfile, delimiter=',')
        spamwriter.writerow([orig_point, orig_sn])
        if sn > float(best_sn):
            spamwriter.writerow(["{}_{}".format(ra, dec), sn])
        else:
            spamwriter.writerow([best_point, best_sn])

    print("${orig_best_file.baseName}".split("_")[0])
    print("Original  {} {}".format(orig_point, orig_sn))
    if sn > float(best_sn):
        print("Best      {} {}".format("{}_{}".format(ra, dec), sn))
    else:
        print("Best      {} {}".format(best_point, best_sn))
    """
}

process publish_best_pointing {
    publishDir params.final_dir, mode: 'copy'

    input:
    file fits

    output:
    file '*' includeInputs true

    """
    echo outputing ${fits}
    """
}

workflow find_pos {
    take:
        pointings
        pre_beamform_1
        pre_beamform_2
        pre_beamform_3
        fwhm
        orig_pointing
        pulsar_pointings
    main:
        beamform( pre_beamform_1,
                  pre_beamform_2,
                  pre_beamform_3,
                  pointings )
        prepfold( // combine pulsar names with fits files
                  pulsar_pointings.map{ it -> [ it[1], it[0] ] }.concat(beamform.out[3]).\
                  // group by pointing
                  groupTuple().map{ it -> [ it[0], it[1][1], it[1][0] ] } )
        pdmp( // combine bestprof and fits files
              prepfold.out[0].map{ it -> [it.baseName.split("_pos")[0].split("${params.obsid}_")[1], it ] }.concat(beamform.out[3])
              // group by pointing
              groupTuple().map{ it -> [ it[0], it[1][0], it[1][1] ] } )

        if ( params.no_pdmp ) {
            if ( params.pulsar == 0 ) {
                bestprof_or_pdmp = prepfold.out[0].flatten().map{ it -> [ "cand", it ] }.groupTuple()
            }
            else {
                bestprof_or_pdmp = prepfold.out[0].flatten().map{ it -> [ "J"+it.baseName.split(".pfd")[0].split("PSR_")[1], it ] }.groupTuple()
            }
        }
        else {
            if ( params.pulsar == 0 ) {
                bestprof_or_pdmp = pdmp.out[1].flatten().map{ it -> [ "cand", it ] }.groupTuple()
            }
            else {
                bestprof_or_pdmp = pdmp.out[1].flatten().map{ it -> [ "J"+it.baseName.split("_pdmp")[0].split("PSR_")[1], it ] }.groupTuple()
            }
        }
        bestgridpos( bestprof_or_pdmp.combine(fwhm).combine(orig_pointing.toList()) )
        //tuple val(pulsar), file(posn_or_bestprof), val(fwhm), val(orig_pointing)
    emit:
        bestgridpos.out[0] // label and new pointing
        bestgridpos.out[2] // orig and best pointing SN file
        beamform.out[1] // fits files of first grid
        prepfold.out[0].concat(prepfold.out[1]).flatten().map{ it -> [ it.baseName.split("_pos")[0], it ]} // presto outputs
        pdmp.out[0].concat(pdmp.out[1], pdmp.out[2]).flatten().map{ it -> [ it.baseName.split("_pos")[0], it ]} // dspsr outputs

}

workflow {
    pre_beamform()
    fwhm_calc( pre_beamform.out[1].map{ it -> it[0] }.collect() )

    // work out intial pointings
    if ( params.pointing_grid ) {
        grid( pointing_grid.map{ it -> [ "cand", it ] }.combine( fwhm_calc.out.splitCsv().flatten() ) )
        pointings = grid.out.splitCsv().map{ it -> it[1] }.toSortedList().flatten().collate( params.max_pointings )
        orig_pointing = pointings.flatten().first()
    }
    else if ( params.pulsar ) {
        get_pulsar_ra_dec()
        grid( get_pulsar_ra_dec.out.splitCsv().combine(fwhm_calc.out.splitCsv().flatten()) )
        pointings = grid.out.splitCsv().map{ it -> it[1] }.toSortedList().flatten().collate( params.max_pointings )
        //pulsar_pointings = grid.out.splitCsv()
        orig_pointing = get_pulsar_ra_dec.out.splitCsv().map{ it -> it[1] }.collect()
    }

    find_pos( pointings,
              pre_beamform.out[0],
              pre_beamform.out[1],
              pre_beamform.out[2],
              fwhm_calc.out.splitCsv().flatten(),
              orig_pointing,
              grid.out.splitCsv() )
    beamform( pre_beamform.out[0],
              pre_beamform.out[1],
              pre_beamform.out[2],
              find_pos.out[0].splitCsv().map{ it -> it[1] }.toSortedList().flatten().collate( params.max_pointings ) )
    prepfold( // combine pulsar names with fits files
              find_pos.out[0].splitCsv().map{ it -> [ it[1].replaceAll(~/\s/,""), it[0] ] }.concat(beamform.out[3]).\
              // group by pointing
              groupTuple().map{ it -> [ it[0], it[1][1], it[1][0] ] } )
    pdmp( // combine bestprof and fits files
          prepfold.out[0].map{ it -> [it.baseName.split("_pos")[0].split("${params.obsid}_")[1], it ] }.concat(beamform.out[3])
          // group by pointing
          groupTuple().map{ it -> [ it[0], it[1][0], it[1][1] ] } )

    // Work out the best pointing
    if ( params.no_pdmp ) {
        if ( params.pulsar == 0 ) {
            bestprof_or_pdmp = prepfold.out[0].flatten().map{ it -> [ "cand", it ] }.groupTuple()
        }
        else {
            bestprof_or_pdmp = prepfold.out[0].flatten().map{ it -> [ "J"+it.baseName.split(".pfd")[0].split("PSR_")[1], it ] }.groupTuple()
        }
    }
    else {
        if ( params.pulsar == 0 ) {
            bestprof_or_pdmp = pdmp.out[1].flatten().map{ it -> [ "cand", it ] }.groupTuple()
        }
        else {
            bestprof_or_pdmp = pdmp.out[1].flatten().map{ it -> [ "J"+it.baseName.split("_pdmp")[0].split("PSR_")[1], it ] }.groupTuple()
        }
    }
    format_output( find_pos.out[1].map{ it -> [ it.baseName.split("_${params.obsid}")[0], it ] }.concat(bestprof_or_pdmp).\
                   groupTuple().map{ it -> it[1] } )

    // Find the best pointing fits file
    publish_best_pointing( // The pointing we want
                           format_output.out.splitCsv( skip: 1 ).map{ it -> [ ("${params.obsid}_" + it[0]).toString(), it[0] ]}.\
                           // All fits files
                           concat(beamform.out[1].concat(find_pos.out[2]).flatten().map{ it -> [ it.baseName.split("_ch")[0], it ] },
                           // Add the presto outputs
                                   find_pos.out[3], prepfold.out[0].concat(prepfold.out[1]).flatten().map{ it -> [ it.baseName.split("_pos")[0], it ]},
                           // Add the dspsr outputs
                                   find_pos.out[4], pdmp.out[0].concat(pdmp.out[1], pdmp.out[2]).flatten().map{ it -> [ it.baseName.split("_pos")[0], it ]}).\
                           // Filter the pointing
                           groupTuple( size: tuple_size ).map{ it -> it[1].tail() } )
}
