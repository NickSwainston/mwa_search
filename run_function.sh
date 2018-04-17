#!/bin/bash
# run command and add relevant data to the job database
# 1st parameter is command to be run (e.g. wsclean)
# 2nd parameter is parameters to that command (e.g. "-j $ncpus")
# 3rd parameter is bs_id
# 4th parameter is DM file int [optional]
if [ "$1" == "rfifind" ] || [ "$1" == "prepsubband" ]; then
if [ -z "$4" ]; then
    rownum=`blindsearch_database.py -m s -c $1 -a "$2" -b $3 -n 1` 
else
    rownum=`blindsearch_database.py -m s -c $1 -a "$2" -b $3 -n 1 -d $4`
fi
$1 $2
echo $1 $2
errcode=$?
blindsearch_database.py -m e -c $1 -r $rownum --errorcode $errcode
echo "blindsearch_database.py -m e -c $1 -r $rownum --errorcode $errcode"
if [ "$errcode" != "0" ]; then
    exit $errcode
fi
else
if [ -z "$4" ]; then
    echo `date +%Y-%m-%d" "%H:%M:%S`",$1,$2,$3,1" >> ${1}_temp_database_file.csv
else
    echo `date +%Y-%m-%d" "%H:%M:%S`",$1,$2,$3,1,$4" >> ${1}_temp_database_file.csv
fi
$1 $2
echo $1 $2
errcode=$?
echo `date +%Y-%m-%d" "%H:%M:%S`",$errcode" >> ${1}_temp_database_file.csv
if [ "$errcode" != "0" ]; then
    exit $errcode
fi
fi
