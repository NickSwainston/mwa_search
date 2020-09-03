#!/usr/bin/env bash

module load git/2.18.0

# Find directory of this build script
ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
ROOT_DIR=$ROOT

cd $ROOT
# VERSION is... the version of the software. This could be hard-coded, or you could determine it from something else.
# The following git branches, if possible.
if [ $# -eq 0 ]
then
    if ! git branch | grep "\\*" | grep "detached"; then
        VERSION=$(git branch | grep "\\*" | cut -d" " -f2-)
    elif git describe --tags &> /dev/null; then
        VERSION=$(git describe --tags)
    else
        VERSION=$(git branch | grep "\\*" | grep "detached" | rev | cut -d" " -f1 | rev | sed "s|)||")
    fi
else
    VERSION=$1
fi
cd -

PACKAGE=mwa_search

if [[ $ROOT == *${PACKAGE}/${PACKAGE} ]]; then
    # This is the prefereed directory structure
    ROOT=${ROOT%"/${PACKAGE}"}
fi

if [ ! -d $ROOT/$VERSION ]; then
    mkdir $ROOT/$VERSION
fi

# Updates to the project version are explained in CHANGELOG.md
python3 setup.py install --prefix=$ROOT/$VERSION/ --single-version-externally-managed --record=record.txt
rm -r build
rm record.txt
find $ROOT/$VERSION/ -type d -exec chmod 775 {} \;
