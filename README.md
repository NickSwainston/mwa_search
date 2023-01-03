# The MWA pulsar search pipeline (deprecated)
This repository is no longer being actively developed. For an up to date version see https://github.com/CIRA-Pulsars-and-Transients-Group/mwa_search.

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/eedca9f0fca94e7cb67b45059eee1da3)](https://www.codacy.com/app/NickSwainston/blindsearch_scripts?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=NickSwainston/blindsearch_scripts&amp;utm_campaign=Badge_Grade)

This repository was written by Nick Swainston to automate pulsar searching using the PRESTO software suite. An explanation of the search procedure can be found on the wiki of the GitHub page. The pipeline uses Nextflow to manage all the required jobs for both beamforming and searching.

## Prerequisites

Requires the [PRESTO](https://github.com/scottransom/presto) software suite, [Nextflow](https://www.nextflow.io/) and [_vcstools_](https://github.com/CIRA-Pulsars-and-Transients-Group/vcstools).

## Installing

On Swinburne's Ozstar supercomputer you can load the module using
```
module use /fred/oz125/software/modulefiles
module load vcstools/master
module load mwa_search/master
```

On Pawsey's Galaxy supercomputer you can load the software on the module using
```
module use /group/mwa/software/modulefiles
module load vcstools
module load mwa_search
```

If you want to install this pipeline on your supercomputer I will likely need to assistant the installation and the writing of your config file. If you would like to try yourself, I recommend installing the _vcstools_ beamformer [docker image](https://hub.docker.com/repository/docker/cirapulsarsandtransients/vcstools) and then installing the python scripts with the setup.py If you would like to attempt it yourself do the following. Make a directory called mwa\_search and then move into it. Then clone the repository and move into the directory it creates. Run the build script using 
```
./build.sh
```
This will move all the python scripts to a directory called master. Then create a module that does the following.
```
export PATH=${PATH}:<your_install_directory>/master
export PYTHONPATH=${PYTHONPATH}:<your_install_directory>/master
```
where \<your\_install_directory\> is the directory where you ran the git clone command, and \<search\_directory\> is where you would like your search pipeline products (make sure this directory exists).

You will also need to edit _config.py_ in _vcstools_ to comply with the modules and directory structure of your supercomputer.

## Developing
If you create a new branch of the git repo then when you use the _build.sh_ script it will make a directory based on your branch name which can be used to test changes to the code without disrupting currently running versions. _mwa\_search\_pipeline.nf_ has an option --mwa_search_version which can use a different module version (which you will have to create) and used to test it, You can then submit a pull request to the GitHub.

## Common Use Cases
All Nextflow scripts have a --help option to explain all the available arguments.

Beamforming a single pointing:
```
beamform.nf --obsid 1272486104 --calid 1272449208 --pointings 19:32:14.05_+10:59:33.38 --all --publish_fits
```

Performing a pulsar search on some fits files:
```
pulsar_search.nf --obsid 1265470568 --pointing 07:42:49.05_-28:22:43.76 --fits_file_dir ./
```

Use the full search pipeline (includes beamforming) on a large number of pointings (listed in the pointings file):
```
mwa_search_pipeline.nf --obsid 1255444104 --calid 1255443816 --pointings 00:34:08.87_-07:21:53.40 --begin 1255444106 --end 1255444705 --summed --pointing_file pointing_file.txt
```
