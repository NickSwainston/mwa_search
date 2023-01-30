# The MWA pulsar search pipeline
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/eedca9f0fca94e7cb67b45059eee1da3)](https://www.codacy.com/app/NickSwainston/blindsearch_scripts?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=NickSwainston/blindsearch_scripts&amp;utm_campaign=Badge_Grade)

This repository was written by Nick Swainston to automate pulsar searching using the PRESTO software suite. An explanation of the search procedure can be found on the wiki of the GitHub page. The pipeline uses Nextflow to manage all the required jobs for both beamforming and searching.

## Documentation

Documentation for `mwa_search` is hosted at [this ReadTheDocs link](https://mwa-search-cira.readthedocs.io/en/latest/).
Source code for this documentation is in the [docs][docs] folder.

## Prerequisites

To run the pipelines contained in the `nextflow` directory requires [Nextflow](https://www.nextflow.io/).
The dependancies of the pipeline are containerised so you will need container software such as Docker installed. If you would like an explanation of the depenancies see the dependancy section below.

## Installing

The repository's scripts can be installed using either:
```
pip install .
```
or
```
python setup.py install.
```

On Swinburne's Ozstar supercomputer, the pipeline is already installed so you can load the module using
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

If you want to install this pipeline on your supercomputer you will need to edit the `nextflow.config` based on your cluster.
To do this, copy one of the `if ( hostname.startsWith("<cluster>") ) {` sections of the config
and edit to describe your clusters' directories structure and dependancy installation.
I will likely need to assistant the installation and the writing of your config file so feel free to make a GitHub issue to ask for assistance.


## Dependancies
The following is all the software we use in the `mwa_search_pipeline.nf`. The following will describe the version of the software we use, the location of Docker images and any changes we have made to the software.

### PRESTO
[PRESTO](https://github.com/scottransom/presto) pulsar search software suite.
We use [this fork](https://github.com/NickSwainston/presto) which includes out custom `ACCEL_sift.py` script.
We use version [v4.0_7ec3c83](https://hub.docker.com/repository/docker/nickswainston/presto/general)


### vcstools
[vcstools](https://github.com/CIRA-Pulsars-and-Transients-Group/vcstools).

I recommend installing the _vcstools_ beamformer [docker image](https://hub.docker.com/repository/docker/cirapulsarsandtransients/vcstools) and then installing the python scripts with the setup.py If you would like to attempt it yourself do the following.

## Developing
If you create a new branch of the git repo then when you use the _build.sh_ script it will make a directory based on your branch name which can be used to test changes to the code without disrupting currently running versions.
_mwa\_search\_pipeline.nf_ has an option --mwa_search_version which can use a different module version (which you will have to create) and used to test it.
You can then submit a pull request to the GitHub.

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
