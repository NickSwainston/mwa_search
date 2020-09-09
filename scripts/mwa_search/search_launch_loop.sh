#!/bin/bash -l 

if [ "$1" != "" ]; then
    OBSID=${1}
else
    read -p "Observation ID: " obsid
    OBSID=${obsid}
fi


if [ "$2" != "" ]; then
    CALID=${2}
else
    read -p "Calibration ID: " calid
    CALID=${calid}
fi

if [ "$3" != "" ]; then
    NAME=${3}
else
    read -p "Name: " name
    NAME=${name}
fi

if [ "$4" != "" ]; then
    BEGIN=${4}
else
    read -p "Begin: " begin
    BEGIN=${begin}
fi

if [ "$5" != "" ]; then
    END=${5}
else
    read -p "End: " end
    END=${end}
fi

cd /fred/oz125/nswainst/pulsar_search
grid.py -o $OBSID -a -b $BEGIN -e $END -d 0.3 -f 0.9 -n 1080 --out_file_name SMART_${NAME}_grid

for SMART_job in $(ls SMART_${NAME}_grid*txt); do
    mkdir -p ${SMART_job%.txt}
    cd ${SMART_job%.txt}
    if [ ! -f "${SMART_job%.txt}_done" ]; then
        mwa_search_pipeline.nf --obsid $OBSID --calid $CALID --pointing_file ../${SMART_job} --begin $BEGIN --end $END -resume \
            --vcstools_version devel --mwa_search_version devel --summed -with-report ${SMART_job%.txt}.html -w ${SMART_job%.txt}_work --out_dir ${SMART_job%.txt}_cands
        errorcode=$?
        echo "Errorcode: $errorcode"
        if [ "$errorcode" != "0" ]; then
            echo "Error in ${SMART_job%.txt}, exiting"
            break
        else
            echo "${SMART_job%.txt} done"
            touch ${SMART_job%.txt}_done
            #rsync --copy-links -zru ${SMART_job%.txt}_cands pulsar-desktop:~/SMART_cand_sorting/${OBSID}; rm -rf ${SMART_job%.txt}_work &
        fi
    else
        echo "${SMART_job%.txt} already finished so skipping"
    fi
    cd ..
done