params.cand_dir = null
params.out_dir = "${params.cand_dir}"

cand_dirs = Channel.fromPath("${params.cand_dir}/*pfd")
    .toList()

process feature_extract {
    container = "cirapulsarsandtransients/pulsarfeaturelab:V1.3.2"
    input:
    file cand_files from cand_dirs

    output:
    file '*.arff' into fex_out
    
    """
    ls
    cat `ls ./ | head -n 1`
    python `which PulsarFeatureLab.py` -d `pwd` -f feature_extraction.arff -t 6 -c 3 --meta --arff
    """

}

process classify {
    input:
    path fex_out

    output:
    file "feature_extraction*" into classifier_files
    

    """
    REALPATH=`realpath ${fex_out}`
    for i in {1..5}; do
        java -jar `which LOTAASClassifier.jar` -m ${LOTAAS_MLC_MODEL_DIR}/V1.3.1model\${i}.model -p `realpath ${fex_out}` -a 1 -d
        mv \${REALPATH%arff}positive feature_extraction_m\${i}.positive
        mv \${REALPATH%arff}negative feature_extraction_m\${i}.negative
    done
    """
}

process sort_detections {
    input:
    file classifier_files

    """
    python /group/mwaops/k_smith/repos/mwa_search/mwa_search/scripts/LOTAAS_wrapper.py --out_dir ${params.out_dir}
    """
}
