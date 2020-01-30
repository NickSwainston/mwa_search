#! /usr/bin/env python3

import argparse
import os
from mwa_search_pipeline import your_slurm_queue_check
from job_submit import submit_slurm

def tar_job_wrapper(hsm_work_dir, file_list, remove=True):
    print(file_list)
    for fn, fits_dir in enumerate(file_list):
        if fits_dir.endswith("\n"):
            fits_dir = fits_dir[:-1]
        dir_name = fits_dir.split("/")[-1]
        #only lets 5 jobs at a time to not flood the transfer (could probably handle more
        if os.path.exists(fits_dir):
            print(dir_name)
            your_slurm_queue_check(max_queue=5, grep='tar', queue='workq')
            commands = []
            commands.append('cd {}'.format(fits_dir))
            commands.append("srun -n 1 -c 1 tar zcvvf - {0} | ssh hpc-hsm.pawsey.org.au 'cat - > {1}/temp_{2}.tar.gz'".format(fits_dir, hsm_work_dir, fn))
            commands.append('errorcode=$?')
            commands.append('if [ $errorcode == "0" ]; then')
            commands.append('   echo mv temp_{0}.tar.gz {1}.tar.gz'.format(fn, dir_name))
            commands.append('   ssh hpc-hsm.pawsey.org.au "mv {0}/temp_{1}.tar.gz {0}/{2}.tar.gz"'.format(hsm_work_dir, fn, dir_name))
            if remove:
                commands.append('   cd ..')
                commands.append('   rm -r {}'.format(dir_name))
            commands.append('fi')
            #TODO and a move command
            submit_slurm('tar_{0}_{1}'.format(dir_name,fn), commands,
                     batch_dir="./",
                     slurm_kwargs={"time": "5:00:00"},
                     queue='copyq',
                     submit=True, export='ALL')
        else:
           print('{} does not exist'.format(dir_name))


#srun -n 1 -c 10 tar zcvvf - 1117643248_00*.fits | ssh hsm 'cat - > /project/mwaops/nswainston/yogesh_low_DM_candiate/temp_${2}.tar.gz
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    Used to put the fits files of MWA pointings into cold storage on galaxy. You must have an sshkey for hpc-hsm.pawsey.org.au
    """)
    parser.add_argument('-w', '--work_dir', type=str, help='The working directory on Pawsey\'s HSM cold storage that you would like to store your files in (please make this foulder before running this script). For example /project/mwaops/nswainston/yogesh_low_DM_candiate/.')
    parser.add_argument('-r', '--no_remove', action="store_true", help='By default this script will remove the directory after the file has succesfully been transfered to cold storage. If you use this option it will NOT delete the foulder after transfering it.')
    default_options = parser.add_argument_group('Default Directory Structure Options', "Options to use if the fits files you would like to transfer are in the directory format /group/mwaops/vcs/<obsid>/pointings/<pointing_ra_dec>/")
    default_options.add_argument('-p','--pointing_file',type=str,help='A text file with a pointing on each line in the format HH:MM:SS_+DD:MM:SS. Such a file can be made using the command: for i in $(search_database.py -m vc --all | grep low | grep -v None); do if [[ $i = 17:* ]]; then echo $i; fi; done >> pointings_to_tar.txt')
    default_options.add_argument("-o", '--obsid', type=int, help='The obsid of the MWA observation')
    full_options = parser.add_argument_group('Full Path Options', 'If you do not have the defualt directory structure please use these options.')
    full_options.add_argument('-d', '--dir_file', type=str, help='The text file that contains full path to a foulder full of fits files on each line. For example /group/mwaops/xuemy/pol_census/1118509648/ which contains fits files')

    args=parser.parse_args()

    #parseing the options
    if not args.work_dir:
        print("No --work_dir option given, please add one. Exiting")
        exit()
    else:
        if args.work_dir.endswith("/"):
            args.work_dir = args.wor_kdir[:-1]

    if args.pointing_file and args.dir_file:
        print("Please us only --pointing_dir or --dir_file, not both. Exiting")
        exit()
    elif args.pointing_file:
        if not args.obsid:
            print("--pointing_file requires the option --obsid but none is given. Exiting")
            exit()
        with open(args.pointing_file) as f:
            temp_list = f.readlines()
        fits_foulder_list = []
        for temp in temp_list:
            fits_foulder_list.append('/group/mwaops/vcs/{0}/pointings/{1}'.format(args.obsid,temp))
    elif args.dir_file:
        with open(args.dir_file) as f:
            fits_foulder_list = f.readlines()
    else:
        print("Please us either --pointing_dir or --dir_file. Exiting")
        exit()

    #TODO add a step that checks if you have an sshkey (that actually works)
    """
    try:
        output = subprocess.check_output('ssh -oNumberOfPasswordPrompts=0 hpc-hsm.pawsey.org.au "echo hello"', shell=True)
    except subprocess.CalledProcessError as error:
        print("Error sshing into hpc-hsm.pawsey.org.au. Pleas emake sure you've set up an sshkey. Exiting")
        print("Full error: " + error.output)
        exit()
    """

    tar_job_wrapper(args.work_dir, fits_foulder_list, remove=(not args.no_remove))
