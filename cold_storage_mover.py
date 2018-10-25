#! /usr/bin/env python

import argparse
import os
from blindsearch_pipeline import your_slurm_queue_check
from job_submit import submit_slurm

def tar_job_wrapper(pointing_file, default_work_dir, obsid=1117643248):
    with open(pointing_file) as f:
        pointing_list = f.readlines()
    print pointing_list
    for pn, pointing in enumerate(pointing_list):
        if pointing.endswith("\n"):
            pointing = pointing[:-1]
        print pointing
        #only lets 5 jobs at a time to not flood the transfer (could probably handle more
        your_slurm_queue_check(max_queue=5, grep='tar')
        fits_dir='/group/mwaops/vcs/{0}/pointings/{1}/'.format(obsid,pointing)
        commands = []
        commands.append('cd {}'.format(fits_dir))
        commands.append("srun -n 1 -c 10 tar zcvvf - {0}_00*.fits | ssh hsm 'cat - > /project/mwaops/nswainston/yogesh_low_DM_candiate/temp_{1}.tar.gz'".format(obsid, pn))
        commands.append('echo $?')
        #TODO and a move command
        submit_slurm('tar_{0}_{1}'.format(pointing,pn), commands,
                 batch_dir="{}tar_batch".format(default_work_dir),
                 slurm_kwargs={"time": "5:00:00", "partition": "workq"},
                 submit=True, export='ALL')


#srun -n 1 -c 10 tar zcvvf - 1117643248_00*.fits | ssh hsm 'cat - > /project/mwaops/nswainston/yogesh_low_DM_candiate/temp_${2}.tar.gz
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    Used to put pointings into cold storage on galaxy.
    """, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-p','--pointing_file',type=str,help='A text file with a pointing on each line in the format HH:MM:SS_+DD:MM:SS. Such a file can be made using the command: for i in $(blindsearch_database.py -m vc --all | grep low | grep -v None); do if [[ $i = 17:* ]]; then echo $i; fi; done >> pointings_to_tar.txt')
    args=parser.parse_args() 
    
    default_work_dir = os.environ['BLINDSEARCH_WORK_DIR']
    tar_job_wrapper(args.pointing_file, default_work_dir)
