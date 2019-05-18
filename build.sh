#!/usr/bin/env bash

# Find directory of this build script
ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
ROOT_DIR=$ROOT

cd $ROOT
# VERSION is... the version of the software. This could be hard-coded, or you could determine it from something else.
# The following git branches, if possible.
if ! git branch | grep "\\*" | grep "detached"; then
    VERSION=$(git branch | grep "\\*" | cut -d" " -f2-)
elif git describe --tags &> /dev/null; then
    VERSION=$(git describe --tags)
else
    VERSION=$(git branch | grep "\\*" | grep "detached" | rev | cut -d" " -f1 | rev | sed "s|)||")
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
PROJECT_VERSION="1.0.0"
VERSION_GIT=$(git log -1 --format=%h)
SEARCH_VERSION="${PROJECT_VERSION}_${VERSION_GIT}"
# Creates version.py to track git version in scripts
echo '__version__ = ''"'${SEARCH_VERSION}'"' > ${ROOT}/${VERSION}/version.py

cp ${ROOT_DIR}/*py $ROOT/$VERSION/

# Fix permissions.
find "$ROOT/$VERSION" -user "$USER" -type d -exec chmod u+rwx,g+rwx,o+rx,o-w {} \;
find "$ROOT/$VERSION" -user "$USER" -type f -exec chmod u+rwx,g+rwx,o+rx,o-w {} \;
