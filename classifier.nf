#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.cand_dir = null
params.out_dir = "${params.cand_dir}"

params.help = false
if ( params.help ) {
    help = """classifier.nf: A pipeline that will use the LOTAAS classifier to sort candidates
             |               into positive and negative detections.
             |Required argurments:
             |  --cand_dir  The directory containing the pfd files [no default]
             |
             |Optional arguments:
             |  --out_dir   The output directory [default: cand_dir]""".stripMargin()
    println(help)
    exit(0)
}

include { classifier } from './classifier_module'

workflow {
    classifier( Channel.fromPath("${params.cand_dir}/*pfd").toList( ))
}
