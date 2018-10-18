# The MWA pulsar search pipeline

This repository was written by Nick Swainston to automate pulsar searching using the PRESTO software suit. It creates batch scripts that are launched onto SLURM job queues and records any errors and reruns jobs that didn't complete succesfully. It can also be used as a beamform wrapper for the process\_vcs.py script from the vcs\_tools repository. 

### Prerequisites

The PRESTO software suite and the vcs\_tools repoistory must be installed on your machine for these scripts to run

### Installing

To install you must put the bin directory into your PATH and PYTHONPATH. So add these to your .bashrc replacing \<your\_install_directory\> with the directory that you clone the repository

```
export PATH=${PATH}:<your_install_directory>/blindsearch_scripts/bin/
export PYTHONPATH=${PYTHONPATH}:<your_install_directory>/blindsearch_scripts/bin/
```

If you would like to use my blindsearch database add:

```
export CMD_BS_DB_DEF_FILE=/group/mwaops/nswainston/blindsearch/.blindsearch_database_defensive.db
```

If you would like your own databse:

```
export CMD_BS_DB_DEF_FILE=<where_you_would_like_you_database>.db
```

then run

```
touch <where_you_would_like_you_database>.db
python init_blindsearch_database.py
```

