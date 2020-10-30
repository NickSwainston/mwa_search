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

grid.py -o $OBSID -a -b $BEGIN -e $END -d 0.3 -f 0.9 -n 1080 --out_file_name SMART_${NAME}_grid

n_grids=$(ls SMART_${NAME}_grid*txt | wc -l)
n_done=0

echo "Looping over $n_grids files"

while [ $n_done -lt $n_grids ]; do
    n_done=0
    for SMART_job in $(ls SMART_${NAME}_grid*txt); do
        mkdir -p ${SMART_job%.txt}
        cd ${SMART_job%.txt}
        if [ -f "${SMART_job%.txt}_done" ]; then
            n_done=$(expr $n_done + 1)
            if [ -d "${SMART_job%.txt}_cands" ]; then
                echo "Syncing ${SMART_job%.txt}_cands"
                rsync --copy-links -zru ${SMART_job%.txt}_cands prometheus:/data/nswainston/SMART_cand_sorting/${OBSID}
                if [ $? != 0 ]; then
                    echo "rsync error exiting"
                    exit
                fi
                echo "Syncing ${SMART_job%.txt}_cands done"
                echo "Deleting ${SMART_job%.txt}_cands"
                rm -rf ${SMART_job%.txt}_cands
                echo "Deleting ${SMART_job%.txt}_cands done"
            fi
            if [ -d "${SMART_job%.txt}_work" ]; then
                echo "Deleting ${SMART_job%.txt}_work"
                rm -rf ${SMART_job%.txt}_work
                echo "Deleting ${SMART_job%.txt}_work done"
            fi
        fi
        # Even if it isn't done start using rysnc update
        if [ -d "${SMART_job%.txt}_cands" ]; then
            echo "Syncing ${SMART_job%.txt}_cands"
            rsync --copy-links -zru ${SMART_job%.txt}_cands prometheus:/data/nswainston/SMART_cand_sorting/${OBSID}
            echo "Syncing ${SMART_job%.txt}_cands done"
        fi
        cd ..
    done
    sleep 600
done
