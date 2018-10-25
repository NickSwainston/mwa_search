#! /usr/bin/env python

import subprocess
import os
import argparse
import textwrap
import urllib
import urllib2
import json
import pipes
import glob
from time import sleep
import datetime
import numpy as np

#vcstools imports
import blindsearch_database
import mwa_metadb_utils as meta
import process_vcs as pvcs
from job_submit import submit_slurm

#python /group/mwaops/nswainston/bin/blindsearch_pipeline.py -o 1133329792 -p 19:45:14.00_-31:47:36.00
#python /group/mwaops/nswainston/bin/blindsearch_pipeline.py -o 1150234552 -p 00:34:08.8703_-07:21:53.409 --pulsar J0034-0721
#python /group/mwaops/nswainston/bin/blindsearch_pipeline.py -o 1099414416 -p 05:34:32_+22:00:53 --pulsar J0534+2200

#1163853320 47 tuck data

def chunks(l, n):
    # For item i in a range that is a length of l,
    for i in range(0, len(l), n):
        # Create an index range for l of n items:
        yield l[i:i+n]

def exists_remote(host, path):
    """Test if a file exists at path on a host accessible with SSH."""
    status = subprocess.call(
        ['ssh', host, 'test -f {}'.format(pipes.quote(path))])
    if status == 0:
        return True
    if status == 1:
        return False
    raise Exception('SSH failed')


def your_slurm_queue_check(max_queue = 80, pbs = False, queue = 'workq', grep=None):
    """
    Checks if you have over 100 jobs on the queue, if so waits until your queue clears
    """
    if pbs:
        submit_line = 'qstat -u $USER '
    else:
        submit_line = 'squeue -u $USER --partition={0}'.format(queue)
    if grep is not None:
        submit_line += ' | grep "{0}"'.format(grep)
    submit_line += ' | wc -l'
    submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
    q_num = ""
    for line in submit_cmd.stdout:
            q_num += line
    print "q: " + str(int(q_num))
    while (int(q_num) > max_queue ):
        print "waiting 5 mins for queue to clear"
        sleep(300)
        submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
        q_num = ""
        for line in submit_cmd.stdout:
                q_num += line
    return


def job_setup_headers(pbs=False, script_test=False, n_omp_threads=1):
    batch_line ='ncpus={0}\n'.format(n_omp_threads) +\
                'export OMP_NUM_THREADS={0}\n'.format(n_omp_threads) +\
                'export CMD_VCS_DB_FILE=/astro/mwaops/vcs/.vcs.db\n' + \
                'export CMD_BS_DB_DEF_FILE=/group/mwaops/nswainston/blindsearch/.blindsearch_database_defensive.db\n'
    if script_test:
        batch_line += """export PATH="$( echo $PATH| tr : '\n' |grep -v /group/mwaops/nswainston/code/blindsearch_scripts/bin/ | paste -s -d: )"\n
                      export PATH=${PATH}:/group/mwaops/nswainston/code/blindsearch_scripts/\n"""
    else:
        batch_line += 'export PATH=${PATH}:/group/mwaops/nswainston/code/blindsearch_scripts/bin/\n'
    return batch_line


def add_database_function(pbs=False, script_test=False, n_omp_threads=1):
    batch_line ='function run\n' +\
                '{\n' +\
                '    # run command and add relevant data to the job database\n' +\
                '    # 1st parameter is command\n' +\
                '    # 2nd parameter is argument\n' +\
                '    # 3rd parameter is bsd_row_num\n' +\
                '    # 4th parameter is script_row_num\n' +\
                '    # 5th parameter is attempt_num\n' +\
                '    echo blindsearch_database.py -m s -c $1 -b $3 -r $4 -a $5 \n' +\
                '    blindsearch_database.py -m s -c $1 -b $3 -r $4 -a $5 \n' +\
                '    srun --export=ALL -n 1 -c $ncpus $1 $2\n' +\
                '    echo $1 $2\n' +\
                '    errcode=$?\n' +\
                '    blindsearch_database.py -m e -c $1 -b $3 -r $4 -a $5 --errorcode $errcode\n' +\
                '    if [ "$errcode" != "0" ]; then\n' +\
                '        exit $errcode\n' +\
                '    fi\n' +\
                '}\n' +\
                '\n'
    batch_line += job_setup_headers(pbs=pbs, script_test=script_test, n_omp_threads=n_omp_threads)
    return batch_line


def add_temp_database_function(job_num, attempt_num, pbs=False, n_omp_threads=1, threads=True):
    #srun removed because all quick jobs are run with a single srun command within a bash script

    batch_line ='function run\n' +\
                '{\n' +\
                '    # run command and add relevant data to the job database\n' +\
                '    # 1st parameter is command to be run (e.g. wsclean)\n' +\
                '    # 2nd parameter is parameters to that command (e.g. "-j $ncpus")\n' +\
                '    # 3rd parameter is bsd_row_num\n' +\
                '    # 4th parameter is DM file int [optional]\n' +\
                '    echo `date +%Y-%m-%d" "%H:%M:%S`",$3,$4,$5" >> ${1}_temp_database_file_'+\
                        str(attempt_num) + '_' + str(job_num) + '.csv\n' +\
                '    $1 $2\n' +\
                '    echo $1 $2\n' +\
                '    errcode=$?\n' +\
                '    echo `date +%Y-%m-%d" "%H:%M:%S`",$errcode" >> ${1}_temp_database_file_'+\
                        str(attempt_num) + '_' + str(job_num) + '.csv\n' +\
                '    if [ "$errcode" != "0" ]; then\n' +\
                '        exit $errcode\n' +\
                '    fi\n' +\
                '}\n' +\
                'ncpus={0}\n'.format(n_omp_threads)
    return batch_line
    
    
def numout_calc(fits_dir, obsid):
    """
    Find the number of time samples of all fits files in the given directory
    """
    dirlist = glob.glob("{0}/{1}_00*.fits".format(fits_dir, obsid))
    numout = 0 
    for d in dirlist:
        submit_line = 'readfile ' + d
        print submit_line
        submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
        for line in submit_cmd.stdout:
            if 'Time per file (sec) =' in line:
                subint = int(line.split('=')[-1])
        numout += int(subint * 1e4)

    if numout%2:
        numout += 1        
    return numout
    
def get_pulsar_dm_p(pulsar):
    #Gets the ra and dec from the output of PSRCAT
    cmd = ['psrcat', '-c', 'dm', pulsar]
    output = subprocess.Popen(cmd,stdout=subprocess.PIPE).communicate()[0]
    temp = []
    lines = output.split('\n')
    for l in lines[4:-1]: 
        columns = l.split()

        if len(columns) > 1:
            dm = columns[1]
    cmd = ['psrcat', '-c', 'p0', pulsar]
    output = subprocess.Popen(cmd,stdout=subprocess.PIPE).communicate()[0]
    temp = []
    lines = output.split('\n')
    for l in lines[4:-1]: 
        columns = l.split()
        if len(columns) > 1:
            p = columns[1]
    return [dm, p]
    
def dm_i_to_file(dm_i):
    if 0 <= dm_i <= 3:
        dm_file = "DM_00" + str(dm_i*2) + "-00" + str((dm_i+1)*2)
    elif dm_i == 4:
        dm_file = "DM_008-010"
    elif 5 <= dm_i <= 48:
        dm_file = "DM_0" + str(dm_i*2) + "-0" + str((dm_i+1)*2)
    elif dm_i == 49:
        dm_file = "DM_098-100"
    elif 50 <= dm_i <= 149:
        dm_file = "DM_" + str(dm_i*2) + "-" + str((dm_i+1)*2)
    else:
        print dm_i
    return dm_file

def process_vcs_wrapper(obsid, begin, end, pointing, args, DI_dir,
                        pointing_dir,relaunch_script,
                        search=False, bsd_row_num=None, nice=100,
                        pulsar_check=None, cal_id=None, vdif=False):
    """
    Does some basic checks and formating before 
    if args.pulsar_file:
        code_comment += using beamforming from process_vcs.py
    """
    #TODO add script_test
    #check queue
    your_slurm_queue_check(queue = 'gpuq')

    #check for incoh file which is used to predict if you have used rfifind
    if not os.path.exists('/group/mwaops/vcs/{0}/incoh'.format(obsid)):
        bf_formats = " -p -i"
    else:
        bf_formats = " -p"
    if vdif:
        bf_formats += " -u"


    #set up and launch beamfroming
    data_dir = '/astro/mwaops/vcs/{0}'.format(obsid)
    product_dir = '/group/mwaops/vcs/{0}'.format(obsid)
    metafits_dir = "{0}/{1}_metafits_ppds.fits".format(data_dir, obsid)
    pvcs.ensure_metafits(data_dir, obsid, metafits_dir)
    job_id_list = pvcs.coherent_beam(obsid, begin, end, data_dir, product_dir,
                  "{0}/batch".format(product_dir), 
                  metafits_dir, 128, pointing, args,
                  bf_formats=bf_formats, DI_dir=DI_dir, calibration_type="rts", nice=nice)
    
    pointing = "{0}_{1}".format(pointing[0],pointing[1])
    dependant_splice_batch(obsid, pointing, product_dir, pointing_dir, job_id_list,
                           bsd_row_num=bsd_row_num, pulsar_check=pulsar_check, 
                           relaunch_script=relaunch_script, cal_id=cal_id)
    return


def dependant_splice_batch(obsid, pointing, product_dir, pointing_dir, job_id_list,
                            bsd_row_num=None, pulsar_check=None, pbs=False, 
                            relaunch_script="echo no relaunch", cal_id=None):
    """
    Launches a script that splices the beamformed files and, where approriate,
    launches the blindsearch pipeline or folds on known pulsars.
    """
    #create a split wrapper dependancy
    splice_wrapper_batch = 'splice_wrapper_{0}_{1}'.format(obsid, pointing)
    commands = []
    commands.append(job_setup_headers(pbs=pbs))
    if bsd_row_num is not None:
        #record beamforming processing time
        commands.append('blindsearch_database.py -m b -b {0} -f mb_{1}'.format(bsd_row_num, pointing))
    if pulsar_check is not None:
        #check_known_pulsars.py uses this to check if it was detected and if so upload it
        commands.append('splice_wrapper.py -o {0} -w {1} -d'.format(obsid, pointing_dir))
        commands.append('cd {0}'.format(pointing_dir))
        for pulsar in pulsar_check:
            commands.append("prepfold -o {0} -noxwin -runavg -noclip -psr {1} -nsub 256 {2}/1*fits".\
                            format(obsid, pulsar, pointing_dir))
            commands.append('chi=`sed "13q;d" {0}_PSR_{1}.pfd.bestprof`'.format(obsid,pulsar))
            commands.append('chi=${chi#*=}')
            commands.append('if [ ${chi%.*} -ge 5 ]; then')
            commands.append('submit_to_database.py -o {0} --cal_id {1} -p {2} --bestprof {0}_PSR_{2}.pfd.bestprof --ppps {0}_PSR_{2}.pfd.ps'.format(obsid, cal_id, pulsar))
            #move files for mengyao to analyse
            pol_census_dir = '/group/mwaops/xuemy/pol_census/{0}/pointing/{1}'.format(obsid,pulsar)
            pol_census_fold_dir = '/group/mwaops/xuemy/pol_census/{0}/pfold/'.format(obsid)
            vcs_pointing_dir = '/group/mwaops/vcs/{0}/pointings/{1}'.format(obsid, pointing)
            if not os.path.exists(pol_census_fold_dir):
               os.makedirs(pol_census_fold_dir)
            commands.append('mv {0}/{1}_PSR_{2}.pfd* {3}'.format(vcs_pointing_dir, 
                                            obsid, pulsar, pol_census_fold_dir))
            if not os.path.exists('/group/mwaops/xuemy/pol_census/{0}/pointing'.format(obsid)):
                os.makedirs('/group/mwaops/xuemy/pol_census/{0}/pointing'.format(obsid))
            commands.append('mv {0} {1}'.format(vcs_pointing_dir, pol_census_dir))
            commands.append('echo "{0}_PSR_{1}.pfd.png is over 5"'.format(obsid,pulsar))
        commands.append("fi")
    elif os.path.exists('/group/mwaops/vcs/{0}/incoh'.format(obsid)):
        #run normally
        commands.append('splice_wrapper.py -o {0} -w {1} -d'.format(obsid, pointing_dir))
    else:
        #make a incoh pointing
        commands.append('mkdir /group/mwaops/vcs/{0}/incoh'.format(obsid))
        commands.append('splice_wrapper.py -o {0} -w {1} -i -d'.format(obsid, pointing_dir))
        commands.append('mv /group/mwaops/vcs/{0}/pointings/{1}/*incoh* /group/mwaops/vcs/{0}/incoh/'.
                            format(obsid, pointing))
        #run rfi job
        if bsd_row_num is None:
            commands.append('{0} -m r -p {1}'.format(relaunch_script, pointing))
        else:
            commands.append('{0} -m r -r {1} -p {2}'.format(relaunch_script, bsd_row_num, pointing))
    #add relaunch script
    if bsd_row_num is None:
        commands.append('{0} -m b -p {1}'.format(relaunch_script, pointing))
    else:
        commands.append('{0} -m b -r {1} -p {2}'.format(relaunch_script, bsd_row_num, pointing))
    
    if os.path.exists("{0}/batch/".format(product_dir)):
        batch_dir = "{0}/batch/".format(product_dir)
    else:
        batch_dir = product_dir
    submit_slurm(splice_wrapper_batch, commands,
                 batch_dir=batch_dir,
                 slurm_kwargs={"time": "5:00:00", "partition": "workq"},
                 submit=True, depend=job_id_list, depend_type='afterany',
                 module_list=["presto/master"])
    return


def beamform(pointing_list, obsid, begin, end, DI_dir, 
             work_dir='/group/nswainston/blindsearch/', relaunch=False,             
             relaunch_script=None, code_comment=None, dm_max=4,
             search=False, bsd_row_num_input=None, incoh=False, 
             pbs=False, pulsar=None, args=None, script_test=False,
             fits_dir_base=None, pulsar_check=None, cal_id=None,
             vdif=False):
    
    for n, line in enumerate(pointing_list):
        if line.startswith("#"):
            continue
        print "Checking pointing {0} out of {1}".format(n+1, len(pointing_list))
        if incoh:
            pointing = "incoh"
        elif fits_dir_base is not None:
            pointing = fits_dir_base.split("/")[-1]
        else:
            ra, dec = line.split(" ")
            if dec.endswith("\n"):
                dec = dec[:-1]
            pointing = ra + "_" + dec
        
        
        #fits dir parsing
        if fits_dir_base is None:
            if pbs:
                fits_dir = '/lustre/projects/p125_astro/DATA/'
                product_dir = fits_dir
            else:
                if incoh:
                    fits_dir='/group/mwaops/vcs/{0}/incoh/'.format(obsid)
                else:
                    fits_dir='/group/mwaops/vcs/{0}/pointings/{1}/'.format(obsid,pointing)
                product_dir = '/group/mwaops/vcs/{0}'.format(obsid)
        else:
            #if pulsar is None:
            #    fits_dir = '{0}/{1}/'.format(fits_dir_base, pointing)
            #else:
            #    fits_dir = '{0}/{1}/'.format(fits_dir_base, pulsar)
            fits_dir = fits_dir_base + "/"
            product_dir = '/group/mwaops/vcs/{0}'.format(obsid)
        #Check if pointing in cold storage
        try :
            exists_remote_check = exists_remote("hsm",
                    "/project/mwaops/nswainston/yogesh_low_DM_candiate/{0}_pointing.tar.gz".\
                    format(pointing))
            if exists_remote_check and len(pointing_list) > 1:
                print "The pointing is in cold storage so assumed it is analysised so not reprocessing"
                continue
        except:
            print "Connection to cold storage failed. Will only check for local files"

        
        #Go through some file checks (true is something is missing)
        path_check = False
        missing_file_check = False
        unspliced_check = False
        missing_chan_list = []
        if os.path.exists(fits_dir):
            #first check is there's already spliced files
            #does check if they have the same start time
            expected_file_num = int( (end-begin)/200 ) + 2
            for fnc in range(1,expected_file_num):
                if not glob.glob(fits_dir+obsid+"_*"+str(fnc)+".fits"):
                    missing_file_check = True
            if fits_dir_base is not None:
                #assumes that maybe they didn't use the whole obs and that's ok
                if len(glob.glob(fits_dir+obsid+"_*.fits")) > 0:
                    missing_file_check = False
            
            if missing_file_check:
                #check if we have any unspliced files
                #there are some so going to resubmit jobs
                beam_meta_data = meta.getmeta(service='obs', params={'obs_id':obsid})
                channels = beam_meta_data[u'rfstreams'][u"0"][u'frequencies']
                
                job_id_list =[]
                unspliced_check = False
                for ch in channels:
                    for ne in range(1,expected_file_num):
                        if not glob.glob(fits_dir+"*_"+obsid+"_ch*"+str(ch)+"_00*"+\
                                str(ne)+".fits"):
                            unspliced_check = True
                            if ch not in missing_chan_list:
                                missing_chan_list.append(ch)
        else:
            path_check = True

        
        if (not path_check and not missing_file_check and not unspliced_check):
            #everything is ok so start blind search database recording
            if search and ((relaunch and len(pointing_list) > 1) or len(pointing_list) == 1):
                #start the blind search database recording
                bsd_row_num = blindsearch_database.database_blindsearch_start(obsid,
                                          pointing, "{0} {1}".format(code_comment,n))
            else:
                bsd_row_num = bsd_row_num_input
        else:
            #need to fix some files then search
            if search:
                #start the blind search database recording
                bsd_row_num = blindsearch_database.database_blindsearch_start(obsid,
                                          pointing, "{0} {1}".format(code_comment,n))
            else:
                bsd_row_num = bsd_row_num_input



        #work out what needs to be done
        if path_check or len(missing_chan_list) == 24:
            # do beamforming
            print "No pointing directory or files for {0} starting beamforming".format(pointing)
            process_vcs_wrapper(obsid, begin, end, [ra,dec], args, DI_dir,
                                fits_dir, relaunch_script, vdif=vdif, 
                                pulsar_check=pulsar_check, cal_id=cal_id,
                                search=search, bsd_row_num=bsd_row_num)
        elif missing_file_check and not unspliced_check:
            #splice files
            print "Splicing the files in {0}".format(pointing)
            dependant_splice_batch(obsid, pointing, product_dir, fits_dir, None,
                           bsd_row_num=bsd_row_num, pulsar_check=pulsar_check, 
                           relaunch_script=relaunch_script, cal_id=cal_id)
 
        elif unspliced_check:
            #resubmit any channels that are incomplete
            print "Some channels missing, resubmitting make beam scripts for {0}".format(pointing)
            if len(pointing_list) > 1:
                your_slurm_queue_check(max_queue = 50, queue='gpuq')
                    
            job_id_list = []
            for ch in missing_chan_list:   
                if os.path.exists("/group/mwaops/vcs/{0}/batch/mb_{1}_ch{2}.batch".\
                                  format(obsid, pointing, ch)):
                    #remove all files to prevent errors
                    remove_files = glob.glob("{}*_{}_ch{:03d}_*fits".format(fits_dir, obsid, ch))
                    for rf in remove_files:
                        os.remove(rf)
                    
                    #resubmit missing channels
                    submit_line = "sbatch /group/mwaops/vcs/{0}/batch/mb_{1}_ch{2}.batch".\
                                  format(obsid, pointing, ch)
                    submit_cmd = subprocess.Popen(submit_line,shell=True,\
                                                    stdout=subprocess.PIPE)
                    print submit_line
                    for line in submit_cmd.stdout:
                        print line,
                        if "Submitted" in line:
                            temp = line.split()
                            job_id_list.append(temp[3])

                else:
                    print "ERROR no batch file found"
            dependant_splice_batch(obsid, pointing, product_dir, fits_dir, job_id_list,
                           bsd_row_num=bsd_row_num, pulsar_check=pulsar_check, 
                           relaunch_script=relaunch_script, cal_id=cal_id)

        else:
            #All files there so the check has succeded and going to start the pipeline
            if search and ((relaunch and len(pointing_list) > 1) or len(pointing_list) == 1):
                print "Fits files available, begining pipeline for {0}".format(pointing)
                sub_dir = pointing + "/" + obsid + "/"
                if len(pointing_list) > 1:
                    your_slurm_queue_check(max_queue = 50)
                prepdata(obsid, pointing, "{0} -p {1}".format(relaunch_script, pointing),
                         work_dir=work_dir, pbs=pbs,
                         bsd_row_num=bsd_row_num, pulsar=pulsar,
                         fits_dir=fits_dir, dm_max=dm_max, script_test=script_test)
            #remove any extra unspliced files
            for fr in glob.glob(fits_dir+"*_"+obsid+"_*.fits"):
                os.remove(fr)

            #Sending off the splice wrapper just for the folding
            #splice_wrapper should fail
            if pulsar_check is not None:
                dependant_splice_batch(obsid, pointing, product_dir, fits_dir, None,
                                       bsd_row_num=None, pulsar_check=pulsar_check, 
                                       cal_id=cal_id)
    return
 

#-------------------------------------------------------------------------------------------------------------
def rfifind(obsid, pointing, sub_dir, relaunch_script,
            work_dir='/group/mwaops/nswainston/blindsearch/', pbs=False,
            n_omp_threads = 20,
            bsd_row_num=None, pulsar=None,
            fits_dir=None, script_test=False):

    if fits_dir == None:
        fits_dir='/group/mwaops/vcs/{0}/pointings/{1}/'.format(obsid,pointing)
           
    #Calculates -numout for prepsubbands
    numout = numout_calc(fits_dir, obsid)
    
    
    rfi_batch = str(bsd_row_num) + '_rfi_{0}'.format(obsid)
    commands = []
    commands.append("source /group/mwaops/PULSAR/psrBash.profile")
    commands.append("ncpus={0}".format(n_omp_threads))
    commands.append("export OMP_NUM_THREADS={0}".format(n_omp_threads))
    if not os.path.exists("{0}/rfi_masks/{1}_rfifind.mask".format(work_dir, obsid)):
        #if there is not already a rfi mask make one
        commands.append('rfifind -ncpus $ncpus -noclip -time 12.0 ' + '-o ' + str(obsid) +\
                        ' -zapchan 0:19,108:127,128:147,236:255,256:275,364:383,384:403,492:511,512:531,620:639,640:659,748:767,768:787,876:895,896:915,1004:1023,1024:1043,1132:1151,1152:1171,1260:1279,1280:1299,1388:1407,1408:1427,1516:1535,1536:1555,1644:1663,1664:1683,1772:1791,1792:1811,1900:1919,1920:1939,2028:2047,2048:2067,2156:2175,2176:2195,2284:2303,2304:2323,2412:2431,2432:2451,2540:2559,2560:2579,2668:2687,2688:2707,2796:2815,2816:2835,2924:2943,2944:2963,3052:3071 ' + fits_dir +\
                        '*incoh*.fits')
        commands.append('mv {0}_rfifind.mask {1}/rfi_masks/'.format(obsid,work_dir))
        commands.append('blindsearch_database.py -c rfifind -m p -b ' +str(bsd_row_num))
        commands.append("prepdata -ncpus $ncpus -dm 0 " "-numout " + str(numout) + " -o " +\
                        str(obsid) + "_DM0.00 " + fits_dir + "*incoh*.fits")
        commands.append('mv {0}_rfifind.mask {1}/rfi_masks/'.format(obsid,work_dir))
        commands.append('mv {0}_DM0.00.dat {1}/rfi_masks/'.format(obsid,work_dir))
        commands.append('mv {0}_DM0.00.inf {1}/rfi_masks/'.format(obsid,work_dir))
    #commands.append("{0} -m p -r {1} -s {2}/{3}".format(relaunch_script, bsd_row_num, pointing, obsid))

    submit_slurm(rfi_batch, commands,
                 batch_dir="{0}{1}/{2}/batch".format(work_dir,pointing,obsid),
                 slurm_kwargs={"time": "2:00:00", "partition": "workq"},#4 hours
                 submit=True, module_list=["presto/master"])
 
    return


#-------------------------------------------------------------------------------------------------------------
def prepdata(obsid, pointing, relaunch_script,
             work_dir='/group/mwaops/nswainston/blindsearch/', pbs=False,
             bsd_row_num=None, pulsar=None,
             n_omp_threads = 20,
             fits_dir=None, dm_max=4, script_test=False):
   
    if fits_dir == None:
        fits_dir='/group/mwaops/vcs/{0}/pointings/{1}/'.format(obsid,pointing)

    #Set up some directories and move to it
    if not os.path.exists(work_dir + pointing):
            os.mkdir(work_dir + pointing)
    if not os.path.exists(work_dir + "/rfi_masks"):
            os.mkdir(work_dir + "/rfi_masks")
    if not os.path.exists("{0}{1}/{2}".format(work_dir, pointing, obsid)): 
            os.mkdir("{0}{1}/{2}".format(work_dir, pointing, obsid))
    if not pulsar == None:
        if not os.path.exists("{0}{1}/{2}/{3}".format(work_dir, pointing, obsid, pulsar)): 
            os.mkdir("{0}{1}/{2}/{3}".format(work_dir, pointing, obsid, pulsar))
        os.chdir("{0}{1}/{2}/{3}".format(work_dir, pointing, obsid, pulsar))
        sub_dir = "{0}/{1}/{2}/".format(pointing, obsid, pulsar)
    else:
        os.chdir(work_dir + pointing + "/" + obsid)
        sub_dir = "{0}/{1}/".format(pointing, obsid)
    
    if not os.path.exists("{0}{1}/batch".format(work_dir, sub_dir)): 
            os.mkdir("{0}{1}/batch".format(work_dir, sub_dir))
    

    
    if not pulsar == None:
        os.chdir("{0}{1}/{2}/{3}".format(work_dir, pointing, obsid, pulsar))
        sub_dir = "{0}/{1}/{2}/".format(pointing, obsid, pulsar)
    else:
        os.chdir("{0}{1}/{2}".format(work_dir, pointing, obsid))
        sub_dir = "{0}/{1}/".format(pointing, obsid)
    
    #Get the centre freq channel and then run DDplan.py to work out the most effective DMs
    print "Obtaining metadata from http://mwa-metadata01.pawsey.org.au/metadata/ for OBS ID: " + str(obsid)
    beam_meta_data = meta.getmeta(service='obs', params={'obs_id':obsid})
    channels = beam_meta_data[u'rfstreams'][u"0"][u'frequencies']
    minfreq = float(min(channels))
    maxfreq = float(max(channels))
    centrefreq = 1.28 * (minfreq + (maxfreq-minfreq)/2) #in MHz
    
    
    if not pulsar == None:
        dm, p = get_pulsar_dm_p(pulsar)
        output = subprocess.Popen(['DDplan.py','-l',str(float(dm) - 1.),
                                   '-d', str(float(dm) + 1.),
                                   '-f',str(centrefreq),
                                   '-b','30.7200067160534',
                                   '-t','0.0001',
                                   '-n','3072',
                                   '-o','dm_temp'],stdout=subprocess.PIPE).communicate()
    else:
        output = subprocess.Popen(['DDplan.py','-l','1',
                                   '-d',str(dm_max),
                                   '-f',str(centrefreq),
                                   '-b','30.7200067160534',
                                   '-t','0.0001',
                                   '-n','3072',
                                   '-o','dm_temp'],stdout=subprocess.PIPE).communicate()
    subprocess.check_call("\n", shell=True)
    os.remove('dm_temp.eps')
    dm_list = []
    print output[0]
    lines = output[0].split('\n')
    for l in lines[13:-4]: 
        columns = l.split()
        dm_list.append(columns)
    
    #dm_list = [['1.500','3.500','0.01','4','200','1']]
          
    #Calculates -numout for prepsubbands
    numout = numout_calc(fits_dir, obsid)

    #Calculates the expected procesing time (conservative) in seconds
    #the number of DMs doesn't appear to greatly affect the processing time so not included in the calculation
    expe_proc_time = 1.7*10**11 / numout 
    
    print dm_list    

    #Create a list of all the commands needed
    dms_per_job = 1024
    commands_list = []
    for dm_line in dm_list:
        dm_start = dm_line[0]
        dm_end = float(dm_line[2]) * float(dm_line[4]) + float(dm_start)
        while ( (dm_end - float(dm_start)) / float(dm_line[2])) > float(dms_per_job) :
            #dedisperse for only 1024 steps
            commands_list.append('-ncpus $ncpus -lodm ' + str(dm_start) +\
                            " -dmstep " + str(dm_line[2]) + " -numdms "+str(dms_per_job)+\
                            " -numout " + str(numout) + " -o " + str(obsid) + " " + fits_dir +\
                            str(obsid) + '_00*.fits')
            dm_start = str(float(dm_start) + (float(dms_per_job) * float(dm_line[2])))
        steps = int((dm_end - float(dm_start)) / float(dm_line[2]))
        #last loop to get the <1024 steps
        commands_list.append('-ncpus $ncpus -lodm ' + str(dm_start) +\
                            " -dmstep " + str(dm_line[2]) + " -numdms "+str(steps)+\
                            " -numout " + str(numout) + " -o " + str(obsid) + " " + fits_dir +\
                            str(obsid) + '_00*.fits')
    #Puts all the expected jobs on the databse
    #blindsearcg_database_script_id_list
    bdsil = blindsearch_database.database_script_list(bsd_row_num, 'prepsubband', commands_list, 
                         n_omp_threads, expe_proc_time)
    #Submit a bunch some prepsubbands to create our .dat files
    job_id_list = []
    for cli, prepdata_command in enumerate(commands_list):
        DM_batch = str(bsd_row_num) + '_DM_' + str(cli)
        commands = []
        commands.append(add_database_function(pbs=pbs, script_test=script_test,
                                              n_omp_threads=n_omp_threads))
        commands.append('cd {0}{1}/{2}'.format(work_dir, pointing, obsid))
        commands.append('run prepsubband "{0}" "{1}" "{2}" "1"'.format(prepdata_command,
                                bsd_row_num, cli))
        
        job_id = submit_slurm(DM_batch, commands,
                         batch_dir="{0}{1}/{2}/batch".format(work_dir,pointing,obsid),
                         slurm_kwargs={"time": "3:00:00", "partition": "workq"},#4 hours
                         submit=True, module_list=["presto/master"])
        job_id_list.append(job_id)    
        


    #make a job that simply restarts this program when all prepsubband jobs are complete
    print "Waiting 5 sec to make sure to the dependancy script works"
    sleep(5)
    job_id_str = ""
    for i in job_id_list:
        job_id_str += ":" + str(i)
    
    DM_depend_batch = str(bsd_row_num) + '_dep_prepsubbands'
    commands = []
    commands.append(add_database_function(pbs=pbs, script_test=script_test, n_omp_threads=n_omp_threads))
    commands.append("{0} -m c --table Prepdata -r {1} --attempt 1".\
                    format(relaunch_script, bsd_row_num))
    
    submit_slurm(DM_depend_batch, commands,
                 batch_dir="{0}{1}/{2}/batch".format(work_dir,pointing,obsid),
                 slurm_kwargs={"time": "20:00", "partition": "workq"},#4 hours
                 submit=True, depend=job_id_str[1:], module_list=["presto/master"])
    return
                
#-------------------------------------------------------------------------------------------------------------
def sort_fft(obsid, pointing, sub_dir, relaunch_script,
             work_dir='/group/mwaops/nswainston/blindsearch/', pbs=False,
             bsd_row_num=None, pulsar=None,
             fits_dir=None, dm_max=4, n_omp_threads = 20, script_test=False):

    if fits_dir == None:
        fits_dir='/group/mwaops/vcs/{0}/pointings/{1}/'.format(obsid,pointing)

    #Makes 90 files to make this all a bit more managable and sorts the files.
    os.chdir(work_dir + "/" + sub_dir)
    if not os.path.exists("over_3_png"):
        os.mkdir("over_3_png")
    if not os.path.exists("other_png"):
        os.mkdir("other_png")
    if not os.path.exists("cand_files"):
        os.mkdir("cand_files")
        
    DIR=work_dir + sub_dir
    os.chdir(DIR)
    all_files = os.listdir(DIR+ "/")
    

    #Set up jobs on database
    expe_proc_time = 180. #giving it a generous 60 seconds as these jobs don't take long at all
    commands_list = []
    dat_files = glob.glob(work_dir + sub_dir + "/*dat")
    
    fft_command = ''
    for fi, dat in enumerate(dat_files):
        fft_command += ' ' + str(dat.split("/")[-1])
        if fi%16 == 15:
            #only 16 ffts can be done at once
            commands_list.append(fft_command)
            fft_command = ''
    if fi%16 != 15:
        commands_list.append(fft_command)
    blindsearch_database.database_script_list(bsd_row_num, 'realfft', commands_list,
                                 n_omp_threads, expe_proc_time)

    #Send off jobs
    error_check('FFT', 0, bsd_row_num, relaunch_script,
                obsid, pointing, pbs=pbs, script_test=script_test, bash_job=True,
                work_dir=work_dir, total_job_time=3600)

    print "Sent off fft jobs"
    return
                
                
#-------------------------------------------------------------------------------------------------------------
def accel(obsid, pointing, sub_dir, relaunch_script,
          work_dir='/group/mwaops/nswainston/blindsearch/', pbs=False,
          bsd_row_num=None, pulsar=None, n_omp_threads=20, script_test=False):
    
    DIR=work_dir + sub_dir
    os.chdir(DIR)

    dir_files = glob.glob("*fft")
    commands_list = []
    for fft_file in dir_files:
        #calc processing time
        expe_proc_time = 300.
        commands_list.append('-ncpus $ncpus -zmax 0 -flo 0.75 -fhi 500 '+\
                              '-numharm 8 {}'.format(fft_file))
    
    blindsearch_database.database_script_list(bsd_row_num, 'accelsearch', commands_list,
                                 n_omp_threads, expe_proc_time)

    #Send off jobs
    error_check('Accel', 0, bsd_row_num, relaunch_script,
                obsid, pointing, pbs=pbs, script_test=script_test, bash_job=True,
                work_dir=work_dir, n_omp_threads=n_omp_threads)

    print "Sent off accel jobs"
    return
       
#-------------------------------------------------------------------------------------------------------------
def fold(obsid, pointing, sub_dir, relaunch_script,
         work_dir='/group/mwaops/nswainston/blindsearch/', pbs=False,
         bsd_row_num=None, pulsar=None, fits_dir=None, n_omp_threads = 20, script_test=False):
    from math import floor
    
    if fits_dir == None:
        fits_dir='/group/mwaops/vcs/{0}/pointings/{1}/'.format(obsid,pointing)
    
    DIR=work_dir + str(sub_dir)
    os.chdir(DIR)
    if not os.path.exists("presto_profiles"):
        os.mkdir("presto_profiles")
    
    #run accel_sift.py to find candidates with DM structure
    if pulsar == None: 
        if pbs:
            submit_line = 'python ~/My-Scripts/ACCEL_sift.py ' + dm_file_orig
            file_loc = 'cand_files/cands_'+dm_file_orig+'.txt'
        else:
            submit_line = 'ACCEL_sift.py {0}/'.format(sub_dir)
            file_loc = '{0}cand_files/cands_{1}.txt'.format(work_dir, sub_dir.replace("/","_"))
    else:
        if pbs:
            submit_line = 'python ~/My-Scripts/ACCEL_sift.py .'
        else:
            submit_line = 'python /group/mwaops/nswainston/bin/ACCEL_sift.py .'
        file_loc = 'ACCEL_sift_cands.txt'
            
    print fits_dir   
    #calcs sn_min for candidates
    numout = numout_calc(fits_dir, obsid)
    from math import sqrt,log
    sn_min = ( sqrt(log(numout)) - 0.88 ) / 0.47
    
    cand_list = [] 
    #TODO temporary fix because ACCEL_Sift doesn't work on gstar    
    if pbs:
        #don't do accelsift
        all_files = os.listdir(DIR+ "/*/")#I may of broke it
        for f in all_files:
            if f.endswith("ACCEL_200") or f.endswith("ACCEL_0"):
                with open(DIR+ "/*/" +f,'rb') as accel_cand_file:
                    lines = accel_cand_file.readlines()
                    for l in lines[3:]:
                        l = l.split()
                        if len(l) == 0:
                            break
                        if float(l[1]) > sn_min:
                            cand_list.append([f,l[0],l[1],f[13:-8],l[5]])
    else:
        os.chdir(work_dir)
        submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
        print submit_cmd.communicate()[0],
        if os.path.exists(file_loc):
            #creat a list of s/n>7 cands
            #[[accel_file,cand,s/n,DM,period(ms)]]
            with open(file_loc,"rb") as sifted_cands:
                lines = sifted_cands.readlines()
                for l in lines:
                    if l.startswith(obsid):
                        cand_line = l.split()
                        if float(cand_line[2]) > sn_min:
                            cand_list.append([cand_line[0].split(':')[0],cand_line[0].split(':')[1],\
                                              cand_line[2],cand_line[1],cand_line[7]])
                        #print cand_line[2]
        #print len(cand_list)
        os.chdir(DIR + "presto_profiles")
    
    #cand_list = [accel_file_name, cand_num, SN, DM, period(ms)]
    print "Sending off jobs with fft sn greater than {}".format(sn_min)
    if len(cand_list) > 0:
        #sort by DM
        from operator import itemgetter
        cand_list.sort(key=itemgetter(3))
        print "Number of cands in this file: " + str(len(cand_list))
        
        expe_proc_time = 5400.
        commands_list = []

        for i,c in enumerate(cand_list):
            accel_file_name, cand_num, SN, cand_DM, period = c
            #through some stuffing around sort the fold into 100 folds per job
            #the fold option using .dat files which is quicker but inaccurate
            #fold_command = 'run "prepfold" "-ncpus $ncpus -n 128 -nsub 128 '+\
            #           "-accelcand "+c[1]+" -accelfile "+c[0]+".cand  -o " +\
            #           c[0][:-8] + " " + c[0][:-8] + '.dat" "'+str(bsd_row_num)+'" "'+str(dm_i)+'"' 
           
            #the fold options that uses .fits files which is slower but more accurate
            commands_list.append('-n 128 -noxwin -noclip -o {0}_{1}_{2} -p {3} -dm {4} -nosearch {5}*.fits'.format(accel_file_name, cand_num, pointing, float(period)/1000.,cand_DM, fits_dir, bsd_row_num))
        blindsearch_database.database_script_list(bsd_row_num, 'prepfold', commands_list,
                                                  n_omp_threads, expe_proc_time)

        #Send off jobs
        error_check('Fold', 0, bsd_row_num, relaunch_script,
                    obsid, pointing, pbs=pbs, script_test=script_test,
                    work_dir=work_dir, n_omp_threads=n_omp_threads)

    

    return



def wrap_up(obsid, pointing, 
            work_dir='/group/mwaops/nswainston/blindsearch/', pbs=False,
            bsd_row_num=None, script_test=False):

    #put all significant candidates into the over directory


    #update final database logs

    #TODO make this the clean up script            
    wrap_batch = str(bsd_row_num) + "_wrap_up" 
    commands = []
    commands.append(add_database_function(pbs=pbs, script_test=script_test, n_omp_threads=n_omp_threads))
    commands.append('cd {0}{1}/{2}/presto_profiles'.format(work_dir, pointing, obsid))
    commands.append('echo "Searching for a profile with a sigma greater than 3"')
    commands.append('count=0')
    commands.append('total=`ls *.ps | wc -l`')
    commands.append('for i in $(ls *.ps); do')
    commands.append('   if (( $count % 100 == 0 )); then')
    commands.append('       echo "$count / $total searched"')
    commands.append('   fi')
    commands.append('   chi=`sed "13q;d" ${i%.ps}.bestprof`')
    commands.append('   chi=${chi#*=}')
    commands.append('   if [ ${chi%.*} -ge 2 ]; then')
    commands.append('       mv "${i%.ps}".png ../over_3_png/"${i%.ps}".png')
    commands.append('       echo "${i%.ps}.png is over 3"')
    commands.append('   fi')
    commands.append('   count=$(($count+1))')
    commands.append('done')
    commands.append('over_sn=`ls ../over_3_png/*.png | wc -l`')
    commands.append('blindsearch_database.py -m w -b ' +str(bsd_row_num) +' --cand_val "$total $over_sn 0"')
    submit_slurm(wrap_batch, commands,
                 batch_dir="{0}{1}/{2}/batch".format(work_dir,pointing,obsid),
                 slurm_kwargs={"time": "2:50:00", "partition": "workq"},#4 hours
                 submit=True, module_list=["presto/master"])
    return


def error_check(table, attempt_num, bsd_row_num, relaunch_script,
                obsid, pointing, pbs=False, script_test=False, bash_job=False,
                work_dir='/group/mwaops/nswainston/blindsearch/', n_omp_threads=1,
                total_job_time=18000.):
    """
    Checkes the database for any jobs that didn't complete (or didn't even start)
    and reruns any errors before continuing to the next step
    """
    subdir = ''
    if table == 'Prepdata':
        next_mode = 's'
        cur_mode = 'p'
    elif table == 'FFT':
        n_omp_threads = 1 #fft is not parrelised
        total_job_time = 6000
        threads = False
        bash_job = True
        next_mode = 'a'
        cur_mode = 's'
    elif table == 'Accel':
        threads = True
        bash_job = True
        next_mode = 'f'
        cur_mode = 'a'
    elif table == 'Fold':
        subdir = 'presto_profiles'
        next_mode = 'w'
        cur_mode = 'f'
    
    print table, bsd_row_num, attempt_num
    total_job_time_str = datetime.timedelta(seconds=total_job_time)
    if attempt_num == 0:
        error_data = blindsearch_database.database_script_check(table, bsd_row_num, 1)
    else:
        error_data = blindsearch_database.database_script_check(table, bsd_row_num, attempt_num)
    print "Number of jobs needed to run {}".format(len(error_data))
    if len(error_data) == 0:
        print "{0} -m {1}".format(relaunch_script, next_mode)
        cmd=subprocess.Popen("{0} -m {1}".format(relaunch_script, next_mode),
                             shell=True,stdout=subprocess.PIPE)
        for line in cmd.stdout:
            print line,
    else:
        print error_data
        presto_command = error_data[0][0]
        if attempt_num != 0 :
            command_list = []
            for er in error_data:
               command_list.append(er[1])
            bdsil = blindsearch_database.database_script_list(bsd_row_num, presto_command,
                      command_list, n_omp_threads, error_data[0][2], attempt=(attempt_num+1)) 

        #submit jobs
        check_job_num = 1
        job_id_list = []

        processing_time = 0.0
        check_batch = "{0}_{1}_check_a{2}_{3}".format(bsd_row_num, presto_command,
                            attempt_num+1, check_job_num)
        commands = []
        bash_commands = []

        if bash_job:
            commands.append(job_setup_headers(pbs=pbs, script_test=script_test,
                                              n_omp_threads=n_omp_threads))
        else:
            commands.append(add_database_function(pbs=pbs, script_test=script_test, 
                                                  n_omp_threads=n_omp_threads))
        commands.append('cd {0}{1}/{2}/{3}'.format(work_dir, pointing, obsid, subdir))
        for ei, er in enumerate(error_data):
            bash_commands.append('run "{0}" "{1}" "{2}" "{3}" "{4}"'.format(presto_command, 
                                er[1], bsd_row_num, ei, (attempt_num+1)))
            processing_time += float(er[2])

            # if next step will go over the expected time limit then create a new script
            if (processing_time + float(er[2])) > total_job_time:
                if bash_job:
                    #if it's a bash job make a file for it to run
                    with open('{0}{1}/{2}/{3}.bash'.format(work_dir,
                              pointing,obsid,check_batch), "w") as srun_bash:
                        srun_bash.write(add_temp_database_function(check_job_num, 
                                        attempt_num + 1, pbs=pbs, n_omp_threads=n_omp_threads,
                                        threads=threads))
                        for sc in bash_commands:
                            srun_bash.write("{}\n".format(sc))
                    
                    commands.append("srun --export=ALL -n 1 -c {0} bash {1}.bash".\
                                    format(n_omp_threads, check_batch))
                    commands.append('blindsearch_database.py -c {0} -m m --file_location {0}_temp_database_file_{1}_{2}.csv\n'.format(presto_command, attempt_num + 1, check_job_num))
                else:
                    commands = commands + bash_commands
                
                job_id = submit_slurm(check_batch, commands,
                         batch_dir="{0}{1}/{2}/batch".format(work_dir,pointing,obsid),
                         slurm_kwargs={"time": total_job_time_str , "partition": "workq"},#4 hours
                         submit=True, module_list=["presto/master"])
                job_id_list.append(job_id)
                
                check_job_num += 1
                processing_time = 0.0
                check_batch = "{0}_{1}_check_a{2}_{3}".format(bsd_row_num, presto_command,
                        attempt_num + 1, check_job_num)
                commands = []
                bash_commands = []

                if bash_job:
                    commands.append(job_setup_headers(pbs=pbs, script_test=script_test,
                                                      n_omp_threads=n_omp_threads))
                else:
                    commands.append(add_database_function(pbs=pbs, script_test=script_test, 
                                                          n_omp_threads=n_omp_threads))
                commands.append('cd {0}{1}/{2}/{3}'.format(work_dir, pointing, obsid, subdir))
        
        #Check there is a command to run
        if len(bash_commands) > 0:
            if bash_job:
                #if it's a bash job make a file for it to run
                with open('{0}{1}/{2}/{3}.bash'.format(work_dir,
                          pointing,obsid,check_batch), "w") as srun_bash:
                    srun_bash.write(add_temp_database_function(check_job_num,
                                    attempt_num + 1, pbs=pbs, n_omp_threads=n_omp_threads,
                                    threads=threads))
                    for sc in bash_commands:
                        srun_bash.write("{}\n".format(sc))
                
                commands.append("srun --export=ALL -n 1 -c {0} bash {1}.bash".\
                                format(n_omp_threads, check_batch))
                commands.append('blindsearch_database.py -c {0} -m m --file_location {0}_temp_database_file_{1}_{2}.csv\n'.format(presto_command, attempt_num + 1, check_job_num))
            else:
                commands = commands + bash_commands
                
            print "submit at end {}".format(check_batch)
            job_id = submit_slurm(check_batch, commands,
                     batch_dir="{0}{1}/{2}/batch".format(work_dir,pointing,obsid),
                     slurm_kwargs={"time": total_job_time_str , "partition": "workq"},#4 hours
                     submit=True, module_list=["presto/master"])
            job_id_list.append(job_id)
        
        #Dependancy job
        print "Waiting 5 sec to make sure to the dependancy script works"
        sleep(5)
        
        check_depend_batch = '{0}_dep_{1}_check_a{2}'.format(bsd_row_num,
                                presto_command, attempt_num +1)
        commands = []
        commands.append(job_setup_headers(pbs=pbs, script_test=script_test,
                                               n_omp_threads=n_omp_threads))
        commands.append("{0} --attempt {1} -m c --table {2}".format(relaunch_script,
                                              attempt_num + 1, table))
        
        submit_slurm(check_depend_batch, commands,
                     batch_dir="{0}{1}/{2}/batch".format(work_dir,pointing,obsid),
                     slurm_kwargs={"time": "20:00", "partition": "workq"},
                     submit=True, depend=job_id_list, depend_type="afterany", 
                     module_list=["presto/master"])
    return




if __name__ == "__main__":
    mode_options = ["b", "r", "p", "t", "a", "f", "c", "w"]
    default_work_dir = os.environ['BLINDSEARCH_WORK_DIR']
    parser = argparse.ArgumentParser(description="""
    Used to automate mass beamforming of MWA data and pulsar searches using the galaxy supercomputer (Ozstar coming soon).
    """, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-o','--observation',type=str,help='The observation ID of the MWA observation')
    parser.add_argument('--pbs',action="store_true",help="PBS queue mode. Not fully implemented yet.")
    parser.add_argument('-t','--test',action="store_true",help="Uses the scripts that haven't been put in bin/ yet. Used to test them before using the build script")

    pointing_options = parser.add_argument_group('Pointing Options')
    pointing_options.add_argument('-p','--pointing',type=str,help='The pointing of the fits file to be searched. The RA and Dec seperated by _ in the format HH:MM:SS_+DD:MM:SS.')
    pointing_options.add_argument("--pulsar_file", default=None, help='Location of a file containting the pointings to be processed. Each line contains a pointing with an RA and Dec seperated by a space in the format "HH:MM:SS +DD:MM:SS". Can be made using grid.py.')
    pointing_options.add_argument('-i','--incoh', action="store_true", help='Uses the incoh fits file location') 
    
    beamform_options = parser.add_argument_group('Beamforming Options')
    beamform_options.add_argument("--DI_dir", default=None, help="Directory containing either Direction Independent Jones Matrices (as created by the RTS) or calibration_solution.bin as created by Andre Offringa's tools.[no default]")
    beamform_options.add_argument('-O','--cal_obs', type=int, help="Observation ID of calibrator you want to process. Used to work out default DI_dir directory.")
    beamform_options.add_argument("-b", "--begin", type=int, help="First GPS time to process [no default]")
    beamform_options.add_argument("-e", "--end", type=int, help="Last GPS time to process [no default]")
    beamform_options.add_argument("-a", "--all", action="store_true",  help="Perform on entire observation span. Use instead of -b & -e.")
    beamform_options.add_argument("--fits_dir", default=None, help="The directory containing the fits files. Only required if the fits files aren't in the standard vcs location.")

    blindsearch_options = parser.add_argument_group('Blindsearch Options')     
    blindsearch_options.add_argument("--search", action="store_true",  help="Continue with the blindsearch pipeline after a successful beamforming check. Default False")
    blindsearch_options.add_argument("--relaunch", action="store_true",  help="Will rerun the blindsearch pipeline even when no beamforming is required. Otherwise assumes that if the beamforming is complete then a search has already been performed..")
    blindsearch_options.add_argument("--attempt", type=int, help="The number of attempts that a script has made. Default is 1.", default=1)
    blindsearch_options.add_argument('-w','--work_dir',type=str,help='Work directory. Default: {}'.format(default_work_dir))
    blindsearch_options.add_argument('-s','--sub_dir',type=str,help='Used by the program to keep track of the sub directory its using')
    blindsearch_options.add_argument('-r','--bsd_row_num',type=int,help='Database row reference number for keeping track of the scripts.')
    blindsearch_options.add_argument('--dm_max',type=int, default = 250,help='DM max searched. Default 250')
    blindsearch_options.add_argument('--pulsar',type=str,help="Used to search for a known pulsar by inputing it's Jname. The code then looks within 1 DM and 15%% of the pulsar's period.")
    blindsearch_options.add_argument('-m','--mode',type=str,help=textwrap.dedent('''Modes used by the pipeline do indicate which step in the process it is up to. The default is beamform mode.
"b" Checks the fits file directory to decide if all files are there or if splicing or beamforming is required.
"r" Uses PRESTO's rfifind to create a RFI mask.
"p" Uses PRESTO's prepdata to dedisperses the the fits files into .dat files. 
"t" Uses PRESTO's realfft to do a fast fourier transform on the .dat files.
"a" Uses PRESTO's accelsearch to search for periodic signals.
"f" Uses PRESTO's prepfold to fold on any significant candidates.
"w" Wraps up the pipeline by moving significant pulse profiles to it's own file for inspection.
"c" Check script that runs after each step to make sure no jobs have failed.'''), default="b")
    blindsearch_options.add_argument("--table", type=str, help="The table name for the blindsearch database to check, similar to the commands used.")
    args=parser.parse_args()
   
    #argument parsing
    if args.mode not in mode_options:
        print "Unrecognised mode [{0}]. Please use a mode from {1}. Exiting".format(args.mode, mode_options)
        quit()

    if not args.observation:
        print "No observation id given. Please supply one using -o. Exiting"
        quit()

    if not (args.pointing or args.pulsar_file or args.incoh):
        print "No pointing option supplied. Please either use -p or --pulsar_file. Exiting"
        quit()

    
    if args.mode == "b":
        if not args.DI_dir:
            if args.cal_obs:
                args.DI_dir = "/group/mwaops/vcs/{0}/cal/{1}/rts/".format(args.obsid, args.cal_obs)
                print "No DI_dir given so assuming {0} is the directory".format(args.DI_dir)
            else:
                print "Please use either --DI_dir to explicitly or use --cal_obs if the DI Jones matrices are in the standard directory. Exiting."
                quit()
        #check begining end times
        if args.all and (args.begin or args.end):
            print "Please specify EITHER (-b,-e) OR -a"
            quit()
        elif args.all:
            args.begin, args.end = meta.obs_max_min(args.observation)
        elif not (args.all or args.begin or args.end):
            print "Please specify the beginning and end times using (-b,-e) OR -a. Exiting"
            quit()

    
    #Default parsing
    obsid = args.observation
    if args.incoh:
        pointing = 'incoh'
        pointing_list = ['incoh']
    else:
        pointing = args.pointing

    if args.work_dir:
        work_dir = args.work_dir
    elif args.pbs:
        work_dir = '/home/nswainst/blindsearch/'
    else:
        work_dir = '/group/mwaops/nswainston/blindsearch/'

    if args.sub_dir:
        sub_dir = args.sub_dir
    elif not args.mode == 'b':
        sub_dir = pointing + "/" + obsid + "/"
    if args.bsd_row_num:
        bsd_row_num = args.bsd_row_num

    if args.pointing:
        if args.pbs:
            fits_dir = '/lustre/projects/p125_astro/DATA/'
        else:
            if args.incoh:
                fits_dir='/group/mwaops/vcs/{0}/incoh/'.format(obsid)
            else:
                fits_dir='/group/mwaops/vcs/{0}/pointings/{1}/'.format(obsid,args.pointing)

    if args.pbs:
        n_omp_threads = 8
    else:
        n_omp_threads = 20
   

    #Base launch of this code (everything except mode and dmfile int)
    relaunch_script = "blindsearch_pipeline.py -o " + str(obsid) + " -w " + work_dir 

    if not args.mode =='b':
        relaunch_script +=  " -s " + str(sub_dir)
    if args.DI_dir:
        relaunch_script += " --DI_dir="+args.DI_dir
    if args.bsd_row_num:
        relaunch_script += " -r " + str(args.bsd_row_num)
    if not args.pulsar == None:
        relaunch_script += " --pulsar " + pulsar
    if args.pbs:
        relaunch_script += " --pbs "
    if args.pointing:
        relaunch_script += " -p " + str(pointing)
    if args.begin and args.end:
        relaunch_script += " -b " + str(args.begin) + " -e " + str(args.end)
    if args.search and args.mode == 'b':
        relaunch_script += " --search"
    if args.dm_max:
        relaunch_script += " --dm_max " + str(args.dm_max)
    if args.incoh:
        relaunch_script +=  " --incoh "
    if args.test:
        relaunch_script +=  " --test"
    if args.fits_dir:
        relaunch_script +=  " --fits_dir " + str(args.fits_dir)


    #work out start and stop times for beamforming
    if args.mode == "b" or args.mode == None:
        if args.pulsar_file:
            with open(args.pulsar_file) as f:
                pointing_list = f.readlines()
        elif args.pointing:
            pointing_list = [args.pointing.replace("_"," ")]
        elif args.fits_dir:
            #gives it a pointing name based on the last dir in the fits 
            #directory as that's normaly a pulsar id or pointing
            if args.fits_dir.endswith("/"):
                args.fits_dir = args.fits_dir[:-1]
            pointing_list = [args.fits_dir.split("/")[-1]]
           
        #If in search mode start up the database entry
        if args.search and not args.bsd_row_num:
            code_comment = raw_input("Please write a comment describing the purpose of this blindsearch. eg testing: ")
            if args.pulsar_file:
                code_comment += " (using: {0}) ".format(args.pulsar_file)
        else:
            code_comment = None
        beamform(pointing_list, obsid, args.begin, args.end, args.DI_dir,
                 work_dir=work_dir, relaunch=args.relaunch, dm_max=args.dm_max,
                 relaunch_script=relaunch_script, code_comment=code_comment,
                 search=args.search, bsd_row_num_input=args.bsd_row_num, incoh=args.incoh,
                 pbs=args.pbs, pulsar=args.pulsar, args=args, script_test=args.test,
                 fits_dir_base=args.fits_dir)
    elif args.mode == "c":
        error_check(args.table, args.attempt, args.bsd_row_num, relaunch_script,
                    obsid, pointing, pbs=args.pbs, script_test=args.test,
                    n_omp_threads=n_omp_threads)
        
    elif args.mode == "r":
        rfifind(obsid, pointing, sub_dir, relaunch_script,
                work_dir=work_dir, pbs=args.pbs, bsd_row_num=args.bsd_row_num, pulsar=args.pulsar,
                n_omp_threads=n_omp_threads, script_test=args.test)
    elif args.mode == "p":
        prepdata(obsid, pointing, relaunch_script,
                 work_dir=work_dir, pbs=args.pbs, bsd_row_num=args.bsd_row_num, pulsar=args.pulsar,
                 n_omp_threads=n_omp_threads, fits_dir=args.fits_dir, dm_max=args.dm_max, 
                 script_test=args.test)
    elif args.mode == "t":
        sort_fft(obsid, pointing, sub_dir, relaunch_script,
                 work_dir=work_dir, pbs=args.pbs, bsd_row_num=args.bsd_row_num, pulsar=args.pulsar,
                 n_omp_threads=n_omp_threads, dm_max=args.dm_max, script_test=args.test)
    elif args.mode == "a":
        accel(obsid, pointing, sub_dir, relaunch_script,
              work_dir=work_dir, pbs=args.pbs, bsd_row_num=args.bsd_row_num, pulsar=args.pulsar,
              n_omp_threads=n_omp_threads, script_test=args.test)
    elif args.mode == "f":
        fold(obsid, pointing, sub_dir, relaunch_script,
             work_dir=work_dir, pbs=args.pbs, bsd_row_num=args.bsd_row_num, pulsar=args.pulsar,
             n_omp_threads=n_omp_threads, fits_dir=args.fits_dir, script_test=args.test)
    elif args.mode == "w":
        wrap_up(obsid, pointing, 
                work_dir=work_dir, pbs=args.pbs,
                bsd_row_num=args.bsd_row_num, script_test=args.test)

            
        
    #blindsearch_pipeline.py -o 1166459712 -p 06:30:00.0_-28:34:00.0
