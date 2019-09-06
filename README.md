# The MWA pulsar search pipeline
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/eedca9f0fca94e7cb67b45059eee1da3)](https://www.codacy.com/app/NickSwainston/blindsearch_scripts?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=NickSwainston/blindsearch_scripts&amp;utm_campaign=Badge_Grade)
This repository was written by Nick Swainston to automate pulsar searching using the PRESTO software suite. An explanation of the search procedure can be found on the wiki of the GitHub page. The pipeline writes batch scripts that are submitted to SLURM job queues, records any errors and reruns jobs that didn't complete successfully. It can also be used as a beamform wrapper for the process\_vcs.py script from the vcs\_tools repository. 

## Prerequisites

Requires the PRESTO software suite and _vcstools_ (private ICRAR repo).

## Installing

On Swinburne's Ozstar supercomputer you can load the module using
```
module use /fred/oz125/software/modulefiles
module load mwa_search/master
```

On Pawsey's Galaxy supercomputer you can load the software on the module using
```
module use /group/mwa/software/modulefiles
module load mwa_search
```

If you want to install this pipeline on your own supercomputer it may be easiest if I assist you as the private _vcstools_ scripts must be installed first. If you would like to attempt it yourself do the following. Make a directory called mwa\_search and then move into it. Then clone the repository and move into the directory it creates. Run the build script using 
```
./build.sh
```
This will move all the python scripts to a directory called master. Then make a module that does the following.
```
export PATH=${PATH}:<your_install_directory>/master
export PYTHONPATH=${PYTHONPATH}:<your_install_directory>/master
export SEARCH_WORK_DIR=<search_directory>
export SEARCH_DB=<where_you_would_like_you_database>.db
```
where \<your\_install_directory\> is the directory where you ran the git clone command and \<search\_directory\> is where you would like your search pipeline products (make sure this directory exists).

then run
```bash
touch <where_you_would_like_you_database>.db
python init_search_database.py
```
You will also need to edit _config.py_ in _vcstools_ to comply with the modules and directory structure of your supercomputer

## Developing
If you create a new branch of the git repo then when you use the _build.sh_ script it will make a directory based on your branch name which can be used to test changes to the code without disrupting currently running versions. _mwa\search\pipeline.py_ has an option -v which can use a different module version (which you will have to create) and used to test it, You can then submit a pull request to the GitHub.

## Common Use Cases

#TODO
