#!/usr/bin/env bash
#changes the PATH and PYTHON PATH so it use the scripts in the main directory
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
echo $DIR
export PATH="$( echo $PATH| tr : '\n' |grep -v ${DIR}/bin/ | paste -s -d: )"
export PATH=${PATH}:${DIR}/
export PYTHONPATH="$( echo $PATH| tr : '\n' |grep -v ${DIR}/bin/ | paste -s -d: )"
export PYTHONPATH=${PATH}:${DIR}/

