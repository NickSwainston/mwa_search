params.cand_dir = null
params.out_dir = "${params.cand_dir}"

process feature_extract {
    
    output:
    file '*.arff' into fex_out
    
    """
    PulsarFeatureLab.py -d ${params.cand_dir} -f ${params.out_dir} -t 6 -c 3 --meta --arff
    """

}

process classifier {

    input:
    file fex_out

    output:
    file "${params.out_dir}/*" into models

    """
    java -jar `which LOTAASClassifier.jar` -m ${params.out_dir}/V1.3.1model5.model -p ${fex_out} -a 5 -d
    """

}

models.view()