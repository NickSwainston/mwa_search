nextflow.preview.dsl = 2
include classifier from './classifier_module'

params.cand_dir = null
params.out_dir = "${params.cand_dir}"

workflow {
    classifier( Channel.fromPath("${params.cand_dir}/*pfd").toList( )
    publish:
        classifier.out[0] to: params.out_dir
        classifier.out[1] to: params.out_dir
}
