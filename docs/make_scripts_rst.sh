# Loops through all the executable scripts and prints all their help text to scripts.rst so it can be used the the sphinx documentation

# mwa_search ------------------------------------------------------------------
echo ".. _mwa_search_scripts:

==================
mwa_search scripts
==================

The following scripts are provided as part of the mwa_search package:
" > mwa_search_scripts.rst

for script in $(ls ../scripts/mwa_search/*py); do
    base_name=${script##*/}
    echo $base_name
    header_dashes=$(perl -E "say '-' x ${#base_name}")

    # Split help into usage, decription and arguments
    python $script -h | awk -v RS= '{print > ("help-" NR ".txt")}'
    usage=$(cat help-1.txt | awk '{print "  " $0}')
    description=$(cat help-2.txt | tr '\n' ' ')
    arguments=$(cat help-3.txt | awk '{print "  " $0}')
    rm help-*txt

    # Print the help and pipe it to a file
    echo ".. _${base_name//_/-}-label:

${base_name}
${header_dashes}

${description}::


${usage}

${arguments}

"  >> mwa_search_scripts.rst
done

# dpp ------------------------------------------------------------------
echo ".. _dpp_scripts:

===========
dpp scripts
===========

The following scripts are provided as part of the dpp package:
" > dpp_scripts.rst

for script in $(ls ../scripts/dpp/*py); do
    base_name=${script##*/}
    echo $base_name
    header_dashes=$(perl -E "say '-' x ${#base_name}")

    # Split help into usage, decription and arguments
    python $script -h | awk -v RS= '{print > ("help-" NR ".txt")}'
    usage=$(cat help-1.txt | awk '{print "  " $0}')
    description=$(cat help-2.txt | tr '\n' ' ')
    arguments=$(cat help-3.txt | awk '{print "  " $0}')
    rm help-*txt

    # Print the help and pipe it to a file
    echo ".. _${base_name//_/-}-label:

${base_name}
${header_dashes}

${description}::


${usage}

${arguments}

"  >> dpp_scripts.rst
done

# nextflow ------------------------------------------------------------------
echo ".. _nextflow_scripts:

================
nextflow scripts
================

The following scripts are provided as part of the nextflow package:
" > nextflow_scripts.rst



for script in $(ls ../nextflow/*nf | grep -v module); do
    base_name=${script##*/}
    echo $base_name
    header_dashes=$(perl -E "say '-' x ${#base_name}")

    # Split help into usage, decription and arguments
    nlines=$(nextflow run $script --help | wc -l)
    description=$(nextflow run $script --help | tail -n $((nlines - 2)))

    # Print the help and pipe it to a file
    echo ".. _${base_name//_/-}-label:

${base_name}
${header_dashes}

${description}

"  >> nextflow_scripts.rst
done

# plotting ------------------------------------------------------------------
echo ".. _plotting_scripts:

================
plotting scripts
================

The following scripts are provided as part of the plotting package:
" > plotting_scripts.rst



for script in $(ls ../scripts/plotting/*py | grep -v module); do
    base_name=${script##*/}
    echo $base_name
    header_dashes=$(perl -E "say '-' x ${#base_name}")

    # Split help into usage, decription and arguments
    python $script -h | awk -v RS= '{print > ("help-" NR ".txt")}'
    usage=$(cat help-1.txt | awk '{print "  " $0}')
    description=$(cat help-2.txt | tr '\n' ' ')
    arguments=$(cat help-3.txt | awk '{print "  " $0}')
    rm help-*txt

    # Print the help and pipe it to a file
    echo ".. _${base_name//_/-}-label:

${base_name}
${header_dashes}

${description}::


${usage}

${arguments}

"  >> plotting_scripts.rst
done