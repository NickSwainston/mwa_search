#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

params.obsid = null
params.calid = null

params.begin = null
params.end = null
params.all = false

params.summed = false
params.channels = null
params.vcstools_version = 'master'
params.mwa_search_version = 'master'

params.didir = "${params.scratch_basedir}/${params.obsid}/cal/${params.calid}/rts"

params.no_combined_check = false

params.help = false
if ( params.help ) {
    help = """beamform.nf: A pipeline that will beamform and splice on all input pointings.
             |Required argurments:
             |  --obsid     Observation ID you want to process [no default]
             |  --calid     Observation ID of calibrator you want to process [no default]
             |  --begin     First GPS time to process [no default]
             |  --end       Last GPS time to process [no default]
             |  --all       Use entire observation span. Use instead of -b & -e. [default: false]
             |
             |Optional arguments:
             |  --summed    Add this flag if you the beamformer output as summed polarisations
             |              (only Stokes I). This reduces the data size by a factor of 4.
             |              [default: False]
             |  --ipfb      Perform an the inverse PFB to produce high time resolution beamformed
             |              vdif files [default: false]
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
//Work out total obs time
if ( params.all ) {
    // an estimation since there's no easy way to make this work
    println("Benchmarks don't need the full observation")
    exit(0)
}
else {
    obs_length = params.end - params.begin + 1
}

range = Channel.from( ['001', '002', '003', '004', '005', '006',\
                       '007', '008', '009', '010', '011', '012',\
                       '013', '014', '015', '016', '017', '018',\
                       '019', '020', '021', '022', '023', '024'] )



bench_max_pointings = 20

//Required temp SSD mem required for gpu jobs
temp_mem = (int) (0.0012 * obs_length * bench_max_pointings + 1)
temp_mem_single = (int) (0.0024 * obs_length + 2)
if ( ! params.summed ) {
    temp_mem = temp_mem * 4
    temp_mem_single = temp_mem_single *4
}

if ( params.summed ) {
    bf_out = " -p -s "
}
else {
    bf_out = " -p "
}

process make_pointings {
    output:
    file "*txt"

    """
    #!/usr/bin/env python

    pointing_list_list = []
    arcsec = 0
    p_ra, p_dec = ["00:00:00.00", "00:00:00.00"]
    dec_deg, dec_min, dec_sec = p_dec.split(":")
    for pn in range(1, ${bench_max_pointings} + 1):
        temp_list = []
        for n in range(1, pn + 1):
            out_dec_sec = float(dec_sec) + arcsec
            out_dec_min = int(dec_min)
            out_dec_deg = int(dec_deg)
            if out_dec_sec > 59:
                out_dec_min = out_dec_min + out_dec_sec // 60
                out_dec_sec = out_dec_sec % 60
            if out_dec_min > 59:
                out_dec_deg = out_dec_deg + out_dec_min // 60
                out_dec_min = out_dec_min % 60
            temp_list.append("{0}_{1:02d}:{2:d}:{3:05.2f}".format(p_ra,
                             out_dec_deg, int(out_dec_min), out_dec_sec))
            arcsec += 1
        pointing_list_list.append(temp_list)
    for pointing_list in pointing_list_list:
        with open('pointings_n{}.txt'.format(len(pointing_list)),'w') as out_file:
            for pointing in pointing_list:
                out_file.write("{}\\n".format(pointing))
    """
}

process make_beam {
    label 'gpu'
    time '24h'
    errorStrategy 'retry'
    maxRetries 1
    maxForks 24

    if ( "$HOSTNAME".startsWith("farnarkle") ) {
        clusterOptions = "--gres=gpu:1  --tmp=${temp_mem}GB"
        scratch '$JOBFS'
        beforeScript "module use ${params.module_dir}; module load vcstools/${params.vcstools_version}"
    }
    else if ( "$HOSTNAME".startsWith("x86") ) {
        clusterOptions = "--gres=gpu:1"
        scratch '/ssd'
        //container = "file:///${params.containerDir}/vcstools/vcstools_${params.vcstools_version}.sif"
        beforeScript "module use ${params.module_dir}; module load vcstools/${params.vcstools_version}"
    }
    else if ( "$HOSTNAME".startsWith("garrawarla") ) {
        clusterOptions = "--gres=gpu:1  --tmp=${temp_mem}GB"
        scratch '/nvmetmp'
        container = "file:///${params.containerDir}/vcstools/vcstools_${params.vcstools_version}.sif"
    }
    else if ( "$HOSTNAME".startsWith("galaxy") ) {
        beforeScript "module use ${params.module_dir}; module load vcstools/${params.vcstools_version}"
    }
    else {
        container = "cirapulsarsandtransients/vcstools:${params.vcstools_version}"
    }

    input:
    each channel_pair
    val utc
    each point
    tuple val(begin), val(end)

    output:
    file "make_beam*txt"

    """
    make_beam -o $params.obsid -b $begin -e $end -a 128 -n 128 \
-f ${channel_pair[0]} -J ${params.didir}/DI_JonesMatrices_node${channel_pair[1]}.dat \
-d ${params.scratch_basedir}/${params.obsid}/combined -P ${point.join(",")} \
-r 10000 -m ${params.scratch_basedir}/${params.obsid}/${params.obsid}_metafits_ppds.fits \
${bf_out} -t 6000 -z $utc &> make_beam_${channel_pair[1]}_n${point.size()}_output.txt
    rm */*fits
    """
}


process make_beam_ipfb {
    label 'gpu'
    time '24h'
    errorStrategy 'retry'
    maxRetries 1
    maxForks 24
    
    if ( "$HOSTNAME".startsWith("farnarkle") ) {
        clusterOptions = "--gres=gpu:1  --tmp=${temp_mem_single}GB"
        scratch '$JOBFS'
        beforeScript "module use ${params.module_dir}; module load vcstools/${params.vcstools_version}"
    }
    else if ( "$HOSTNAME".startsWith("x86") ) {
        clusterOptions = "--gres=gpu:1"
        scratch '/ssd'
        container = "file:///${params.containerDir}/vcstools/vcstools_${params.vcstools_version}.sif"
    }
    else if ( "$HOSTNAME".startsWith("garrawarla") ) {
        clusterOptions = "--gres=gpu:1  --tmp=${temp_mem_single}GB"
        scratch '/nvmetmp'
        container = "file:///${params.containerDir}/vcstools/vcstools_${params.vcstools_version}.sif"
    }
    else if ( "$HOSTNAME".startsWith("galaxy") ) {
        beforeScript "module use ${params.module_dir}; module load vcstools/${params.vcstools_version}"
    }
    else {
        container = "cirapulsarsandtransients/vcstools:${params.vcstools_version}"
    }

    input:
    each channel_pair
    val utc
    each point
    tuple val(begin), val(end)

    output:
    file "make_beam*txt"
    
    """
    make_beam -o $params.obsid -b $begin -e $end -a 128 -n 128 \
-f ${channel_pair[0]} -J ${params.didir}/DI_JonesMatrices_node${channel_pair[1]}.dat \
-d ${params.scratch_basedir}/${params.obsid}/combined -P ${point} \
-r 10000 -m ${params.scratch_basedir}/${params.obsid}/${params.obsid}_metafits_ppds.fits \
-p -v -t 6000 -z $utc &> make_beam_${channel_pair[1]}_IPFB_output.txt
    rm */*fits
    """
}

process make_beam_single {
    label 'gpu'
    time '24h'
    errorStrategy 'retry'
    maxRetries 1
    maxForks 24

    if ( "$HOSTNAME".startsWith("galaxy") ) {
        beforeScript "module use ${params.module_dir}; module load vcstools/nswainston"
    }
    else if ( "$HOSTNAME".startsWith("farnarkle") ) {
        clusterOptions = "--gres=gpu:1  --tmp=${temp_mem_single}GB"
        scratch '$JOBFS'
        container = "file:///${params.containerDir}/vcstools/vcstools_single-pixel_legacy.sif"
    }
    else if ( "$HOSTNAME".startsWith("x86") ) {
        clusterOptions = "--gres=gpu:1"
        scratch '/ssd'
        container = "file:///${params.containerDir}/vcstools/vcstools_single-pixel_legacy.sif"
    }
    else if ( "$HOSTNAME".startsWith("garrawarla") ) {
        clusterOptions = "--gres=gpu:1  --tmp=${temp_mem_single}GB"
        scratch '/nvmetmp'
        container = "file:///${params.containerDir}/vcstools/vcstools_single-pixel_legacy.sif"
    }

    input:
    each channel_pair
    val utc
    val point
    tuple val(begin), val(end)

    output:
    file "make_beam*txt"

    """
    make_beam -o $params.obsid -b $begin -e $end -a 128 -n 128 \
-f ${channel_pair[0]} -J ${params.didir}/DI_JonesMatrices_node${channel_pair[1]}.dat \
-d ${params.scratch_basedir}/${params.obsid}/combined -R ${point.split("_")[0]} -D ${point.split("_")[1]} \
-r 10000 -m ${params.scratch_basedir}/${params.obsid}/${params.obsid}_metafits_ppds.fits \
-p -z $utc &> make_beam_${channel_pair[1]}_single-pixel_output.txt
    rm *fits
    """
}

process calc_beamformer_benchmarks {
    input:
    file files

    output:
    stdout()
    
    """
    calc_beamformer_benchmarks.py --max_pointing_num ${bench_max_pointings}
    """
}

include { pre_beamform } from './beamform_module'

single_pointing = Channel.of("00:00:00.00_00:00:00.00")

workflow {
    pre_beamform()
    make_pointings()
    make_beam( pre_beamform.out[1].flatten().merge(range),\
               pre_beamform.out[2],\
               make_pointings.out.flatten().map{ it -> it.splitCsv().collect().flatten() },\
               pre_beamform.out[0] )
    make_beam_single( pre_beamform.out[1].flatten().merge(range),\
                      pre_beamform.out[2],\
                      single_pointing,\
                      pre_beamform.out[0] )
    make_beam_ipfb( pre_beamform.out[1].flatten().merge(range),\
                    pre_beamform.out[2],\
                    single_pointing,\
                    pre_beamform.out[0] )
    calc_beamformer_benchmarks( make_beam.out.concat(make_beam_single.out, make_beam_ipfb.out).collect() )
    calc_beamformer_benchmarks.out.view()
}
