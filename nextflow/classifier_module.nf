nextflow.preview.dsl = 2


params.out_dir = "${params.search_dir}/${params.obsid}_candidates"

process feature_extract {
    //container = "cirapulsarsandtransients/pulsarfeaturelab:V1.3.2"
    label 'cpu'
    errorStrategy 'retry'
    maxRetries 1
    
    input:
    file pfd_files

    output:
    file "*.arff"
    file "*pfd*" includeInputs true
    
    if ( "$HOSTNAME".startsWith("x86") ) {
        container = 'lofar_pulsar_ml.sif'
    }
        //container = 'nickswainston/lofar_pulsar_ml'
    else {
        beforeScript "module use $params.module_dir; module load PulsarFeatureLab/V1.3.2"
    }

    """
    ls
    cat `ls ./ | head -n 1`
    python `which PulsarFeatureLab.py` -d `pwd` -f feature_extraction.arff -t 6 -c 3 --meta --arff
    """

}

process classify {
    input:
    path fex_out
    file pfd_files

    output:
    file "feature_extraction*"
    file "*pfd*" includeInputs true
    
    if ( "$HOSTNAME".startsWith("x86") ) {
        container = 'lofar_pulsar_ml.sif'
    }
        //container = 'nickswainston/lofar_pulsar_ml'
    else {
        beforeScript "module use $params.module_dir; module load LOTAASClassifier/master"
    }

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
    publishDir params.out_dir, mode: 'copy'

    input:
    file classifier_files
    file pfd_files

    output:
    file "positive_detections/*" optional true
    file "negative_detections/*" optional true

    if ( "$HOSTNAME".startsWith("x86") ) {
        container = 'lofar_pulsar_ml.sif'
    }
        //container = 'nickswainston/lofar_pulsar_ml'
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
        pfd_files
    main:
        feature_extract( pfd_files )
        classify( feature_extract.out[0],\
                  feature_extract.out[1] )
        sort_detections( classify.out[0],\
                         classify.out[1] )//pfd_files )
    emit:
        sort_detections.out[0]
        sort_detections.out[1]
}
