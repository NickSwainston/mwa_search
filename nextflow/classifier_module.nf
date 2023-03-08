

process feature_extract {
    label 'cpu'
    label 'lofar_feature_lab'

    time '1h'
    errorStrategy 'retry'
    maxRetries 1

    input:
    tuple path(pfd), path(bestprof), path(ps), path(png)

    output:
    tuple path(pfd), path(bestprof), path(ps), path(png), path("*.arff")

    """
    ls
    cat `ls ./ | head -n 1`
    python `which PulsarFeatureLab.py` -d `pwd` -f feature_extraction.arff -t 6 -c 3 --meta --arff
    """

}

process classify {
    label 'lofar_ml'

    input:
    tuple path(pfd), path(bestprof), path(ps), path(png), path(fex_out)

    output:
    tuple path(pfd), path(bestprof), path(ps), path(png), path(fex_out), path("feature_extraction*")

    """
    REALPATH=`realpath ${fex_out}`
    for i in {1..5}; do
        java -jar `which LOTAASClassifier.jar` -m \${LOTAAS_MLC_MODEL_DIR}/V1.3.1model\${i}.model -p `realpath ${fex_out}` -a 1 -d
        if [ -f "\${REALPATH%arff}positive" ]; then
            mv \${REALPATH%arff}positive feature_extraction_m\${i}.positive
        fi
        if [ -f "\${REALPATH%arff}negative" ]; then
            mv \${REALPATH%arff}negative feature_extraction_m\${i}.negative
        fi
    done
    """
}

process sort_detections {
    label 'lofar_feature_lab'

    publishDir params.out_dir, mode: 'copy', enabled: params.publish_all_classifer_cands

    input:
    tuple path(pfd), path(bestprof), path(ps), path(png), path(fex_out), path(classifier_files)

    output:
    path "positive_detections/*", optional: true, emit: positive
    path "negative_detections/*", optional: true, emit: negative

    """
    LOTAAS_wrapper.py
    if [ -f LOTAAS_positive_detections.txt ]; then
        mkdir positive_detections
        for i in \$(cat LOTAAS_positive_detections.txt); do
            mv \${i}* positive_detections
        done
    fi
    if [ -f LOTAAS_negative_detections.txt ]; then
        mkdir negative_detections
        for i in \$(cat LOTAAS_negative_detections.txt); do
            mv \${i}* negative_detections
        done
    fi
    """
}


workflow classifier {
    take:
        presto_candiates
    main:
        // Collate into groups of 30 candidates
        collated_cands = presto_candiates.transpose().collate( 30 ).map{ it.transpose() }
        feature_extract( collated_cands )
        classify( feature_extract.out )
        sort_detections( classify.out )
    emit:
        positive = sort_detections.out.positive
        negative = sort_detections.out.negative
}
