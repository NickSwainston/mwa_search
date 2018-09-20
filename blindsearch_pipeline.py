#! /usr/bin/env python

import subprocess
import os
import argparse
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


def your_slurm_queue_check(max_queue = 80, pbs = False, queue = 'workq'):
    """
    Checks if you have over 100 jobs on the queue, if so waits until your queue clears
    """
    if pbs:
        submit_line = 'qstat -u $USER | wc -l'
    else:
        submit_line = 'squeue -u $USER --partition={0}| wc -l'.format(queue)
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

def add_database_function(pbs):
    batch_line ='function run\n' +\
                '{\n' +\
                '    # run command and add relevant data to the job database\n' +\
                '    # 1st parameter is command to be run (e.g. wsclean)\n' +\
                '    # 2nd parameter is parameters to that command (e.g. "-j $ncpus")\n' +\
                '    # 3rd parameter is bsd_row_num\n' +\
                '    # 4th parameter is DM file int [optional]\n' +\
                '    if [ -z "$4" ]; then\n' +\
                '        rownum=`blindsearch_database.py -m s -c $1 -a "$2" -b $3 -n 1` \n' +\
                '    else\n' +\
                '        rownum=`blindsearch_database.py -m s -c $1 -a "$2" -b $3 -n 1 -d $4`\n' +\
                '    fi\n' +\
                '    srun -n 1 -c $ncpus $1 $2\n' +\
                '    echo $1 $2\n' +\
                '    errcode=$?\n' +\
                '    blindsearch_database.py -m e -c $1 -r $rownum --errorcode $errcode\n' +\
                '    echo "blindsearch_database.py -m e -c $1 -r $rownum --errorcode $errcode"\n' +\
                '    if [ "$errcode" != "0" ]; then\n' +\
                '        exit $errcode\n' +\
                '    fi\n' +\
                '}\n'
    return batch_line


def add_temp_database_function(pbs, threads=True):
    batch_line ='function run\n' +\
                '{\n' +\
                '    # run command and add relevant data to the job database\n' +\
                '    # 1st parameter is command to be run (e.g. wsclean)\n' +\
                '    # 2nd parameter is parameters to that command (e.g. "-j $ncpus")\n' +\
                '    # 3rd parameter is bsd_row_num\n' +\
                '    # 4th parameter is DM file int [optional]\n' +\
                '    if [ -z "$4" ]; then\n' +\
                '        echo `date +%Y-%m-%d" "%H:%M:%S`",$1,$2,$3,1" >> ${1}_temp_database_file.csv\n' +\
                '    else\n' +\
                '        echo `date +%Y-%m-%d" "%H:%M:%S`",$1,$2,$3,1,$4" >> ${1}_temp_database_file.csv\n' +\
                '    fi\n'
    #if threads:
    #    batch_line+='    srun -n 1 -c $ncpus $1 $2\n'
    #else:
    #    batch_line+='    srun -n 1 $1 $2\n'
    batch_line+='    $1 $2\n' +\
                '    echo $1 $2\n' +\
                '    errcode=$?\n' +\
                '    echo `date +%Y-%m-%d" "%H:%M:%S`",$errcode" >> ${1}_temp_database_file.csv\n' +\
                '    if [ "$errcode" != "0" ]; then\n' +\
                '        exit $errcode\n' +\
                '    fi\n' +\
                '}\n'
    return batch_line
    
    
def numout_calc(DIR):
    #DIR = '/astro/mwaops/vcs/1133329792/pointings/19:45:14.00_-31:47:36.00/'

    dirlist =[]
    for filen in os.listdir(DIR):
        if filen.endswith(".fits") and filen.startswith("1"):
            dirlist.append(filen)

    numout = 0 
    for d in dirlist:
        readfile_out = []
        submit_line = 'readfile ' + DIR + '/' + d
        submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
        for line in submit_cmd.stdout:
            temp = line.split('=')
            readfile_out.append(temp)
         
        array = []       
        for line in readfile_out:
            temp = []
            for r in line:
                rr = r
                rr = rr.replace(' ','')
                rr = rr.replace('\n','')
                temp.append(rr)
            array.append(temp)
            
        #print array
             
          
        for a in array:
            if a[0] == 'Timeperfile(sec)':
                subint = int(a[1])
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

def process_vcs_wrapper(obsid, begin, end, pointing, args, DI_dir,pointing_dir,\
                        check, search, bsd_row_num, relaunch_script, nice = 100):
    """
    Does some basic checks and formating before 
    if args.pulsar_file:
        code_comment += using beamforming from process_vcs.py
    """
    #check queue
    your_slurm_queue_check(queue = 'gpuq')

    #check for incoh file which is used to predict if you have used rfifind
    if not os.path.exists('/group/mwaops/vcs/{0}/incoh'.format(obsid)):
        bf_formats = " -p -i"
    else:
        bf_formats = " -p"


    #set up and launch beamfroming
    data_dir = '/astro/mwaops/vcs/{0}'.format(obsid)
    product_dir = '/group/mwaops/vcs/{0}'.format(obsid)
    job_id_list = pvcs.coherent_beam(obsid, begin, end, data_dir, product_dir,
                  "{0}/batch".format(product_dir), 
                  "{0}/{1}_metafits_ppds.fits".format(data_dir, obsid), 128, pointing, args,
                  bf_formats=bf_formats, DI_dir=DI_dir, calibration_type="rts", nice=nice)
    
    #get a job dependancy string
    job_id_str = ""
    for i in job_id_list:
        job_id_str += ":" + str(i)
    
    pointing = "{0}_{1}".format(pointing[0],pointing[1])
    #create a split wrapper dependancy
    splice_wrapper_batch = 'splice_wrapper_{0}_{1}'.format(obsid, pointing)
    commands = []
    for f in glob.glob("{0}/batch/mb_{1}*.batch".format(product_dir,pointing)):
        commands.append('blindsearch_database.py -m b -b ' +str(bsd_row_num) + " -f " + str(f)[:-6])
    commands.append('blindsearch_database.py -c make_beam -m p -b ' +str(bsd_row_num))
    if os.path.exists('/group/mwaops/vcs/{0}/incoh'.format(obsid)):
        commands.append('splice_wrapper.py -o {0} -w {1} -d'.format(obsid, pointing_dir))
        commands.append('{0} -m b -r {1} -p {2}'.format(relaunch_script, bsd_row_num, pointing))
    else:
        commands.append('mkdir /group/mwaops/vcs/{0}/incoh'.format(obsid))
        commands.append('splice_wrapper.py -o {0} -w {1} -i -d'.format(obsid, pointing_dir))
        commands.append('mv /group/mwaops/vcs/{0}/pointings/{1}/*incoh* /group/mwaops/vcs/{0}/incoh/'.
                            format(obsid, pointing))
        commands.append('{0} -m r -r {1} -p {2}'.format(relaunch_script, bsd_row_num, pointing))
        commands.append('{0} -m b -r {1} -p {2}'.format(relaunch_script, bsd_row_num, pointing))
    submit_slurm(splice_wrapper_batch, commands,
                 batch_dir="{0}/batch".format(product_dir),
                 slurm_kwargs={"time": "1:00:00", "partition": "workq"},
                 submit=True, depend=job_id_str[1:])
 
    return


def beamform(pointing_list, obsid, begin, end, DI_dir, 
             work_dir='/group/nswainston/blindsearch/', relaunch=False,             
             relaunch_script=None, code_comment=None, dm_max=4,
             search=False, bsd_row_num=None, incoh=False, 
             pbs=False, pulsar=None, check=0, args=None):
    
    for n, line in enumerate(pointing_list):
        if line.startswith("#"):
            continue
        print "Checking pointing {0} out of {1}".format(n+1, len(pointing_list))
        if incoh:
            pointing = "incoh"
        else:
            ra, dec = line.split(" ")
            if dec.endswith("\n"):
                dec = dec[:-1]
            pointing = ra + "_" + dec
        
        #fits dir parsing
        if pbs:
            fits_dir = '/lustre/projects/p125_astro/DATA/'
        else:
            if incoh:
                fits_dir='/group/mwaops/vcs/{0}/incoh/'.format(obsid)
            else:
                fits_dir='/group/mwaops/vcs/{0}/pointings/{1}/'.format(obsid,pointing)

        #Check if pointing in cold storage
        try :
            exists_remote_check = exists_remote("hsm",
                    "/project/mwaops/nswainston/yogesh_low_DM_candiate/{0}_pointing.tar.gz".\
                    format(pointing))
            if exists_remote_check:
                print "The pointing is in cold storage so assumed it is analysised so not reprocessing"
                continue
        except:
            print "Connection to cold storage failed. Will only check for local files"
        
        if os.path.exists(fits_dir):
            #first check is there's already spliced files
            expected_file_num = int( (end-begin)/200 ) + 2
            
            missing_file_check = False
            for fnc in range(1,expected_file_num):
                if not glob.glob(fits_dir+obsid+"_*"+str(fnc)+".fits"):
                    missing_file_check = True
            
            
            if missing_file_check:
                #check if we have any unspliced files
                print fits_dir+"*_"+obsid+"_*.fits"
                if glob.glob(fits_dir+"*_"+obsid+"_*.fits"):
                    #there are some so going to resubmit jobs
                    beam_meta_data = meta.getmeta(service='obs', params={'obs_id':obsid})
                    channels = beam_meta_data[u'rfstreams'][u"0"][u'frequencies']
                    
                    job_id_list =[]
                    unspliced_check = False
                    for ch in channels:
                        channel_check = False
                        for n in range(1,expected_file_num):
                            if not glob.glob(fits_dir+"*_"+obsid+"_ch*"+str(ch)+"_00*"+\
                                    str(n)+".fits"):
                                channel_check = True
                                unspliced_check = True
                        if channel_check:
                            #missing some files for that channel so resumbit script
                            if os.path.exists("/group/mwaops/vcs/"+obsid+"/batch/mb_"+ra+"_"+\
                                    dec+"_ch"+str(ch)+".batch"):
                                #delete files of that channel
                                files_list = os.listdir(fits_dir)
                                for f in files_list:
                                    if ("ch0"+str(ch) in f) or ("ch"+str(ch) in f):
                                        os.remove(fits_dir + f)
                                submit_line = "sbatch /group/mwaops/vcs/"+obsid+"/batch/mb_"+ra+\
                                              "_"+dec+"_ch"+str(ch)+".batch"
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

                    if unspliced_check:
                        #splice wraps them when they're done
                        sleep(1)
                        job_id_str = ""
                        for j in job_id_list:
                            job_id_str += ":" + str(j)
                        splice_wrapper_batch = 'splice_wrapper_{0}_{1}'.format(obsid, pointing)
                        commands = []
                        commands.append('splice_wrapper.py -o {0} -w {1} -d'.format(obsid, fits_dir))
                        if search:
                            #if args.bsd_row_num:
                            bsd_row_num = blindsearch_database.database_blindsearch_start(obsid,
                                                 pointing, "{0} {1}".format(code_comment,n))
                            commands.append('{0} -m b -r {1} -p {2}'.format(relaunch_script,\
                                                              bsd_row_num, pointing))
                        submit_slurm(splice_wrapper_batch, commands,
                                     batch_dir='/group/mwaops/vcs/{0}/batch'.format(obsid),
                                     slurm_kwargs={"time": "1:00:00", "partition": "workq"},
                                     submit=True, depend=job_id_str[1:])
 
                    else:
                        #If only unspliced files then splice
                        splice_wrapper_batch = 'splice_wrapper_{0}_{1}'.format(obsid, pointing)
                        commands = []
                        commands.append('splice_wrapper.py -o {0} -w {1} -d'.format(obsid, fits_dir))
                        if search:
                            if not bsd_row_num:
                                bsd_row_num = blindsearch_database.database_blindsearch_start(obsid,
                                                  pointing, "{0} {1}".format(code_comment,n))
                            commands.append('{0} -m b -r {1} -p {2}'.format(relaunch_script,\
                                                            bsd_row_num, pointing))
                        submit_slurm(splice_wrapper_batch, commands,
                                     batch_dir='/group/mwaops/vcs/{0}/batch'.format(obsid),
                                     slurm_kwargs={"time": "1:00:00", "partition": "workq"},
                                     submit=True)
 
 
                else:
                    print "No files in "+ra+"_"+dec+" starting beamforming"
                    if search:
                        if not bsd_row_num:
                            bsd_row_num = blindsearch_database.database_blindsearch_start(obsid,
                                                    pointing, "{0} {1}".format(code_comment,n))
                    process_vcs_wrapper(obsid, begin, end, [ra,dec], args, DI_dir,\
                                     fits_dir, check, search, bsd_row_num, relaunch_script)
                              
            else:
                #All files there so the check has succeded and going to start the pipeline
                if search and not relaunch:
                    sub_dir = pointing + "/" + obsid + "/"
                    if not bsd_row_num:
                        bsd_row_num = blindsearch_database.database_blindsearch_start(obsid,
                                            pointing, "{0} {1}".format(code_comment,n))
                    if len(pointing_list) > 1:
                        your_slurm_queue_check(max_queue = 300)
                    prepdata(obsid, pointing, "{0} -p {1}".format(relaunch_script, pointing),
                             work_dir=work_dir, pbs=pbs,
                             bsd_row_num=bsd_row_num, pulsar=pulsar,
                             fits_dir=fits_dir, dm_max=dm_max)
                #remove any extra unspliced files
                for fr in glob.glob(fits_dir+"*_"+obsid+"_*.fits"):
                    os.remove(fr)
        else:
            # do beamforming
            print "No pointing directory for "+ra+"_"+dec+" starting beamforming"
            if search:
                bsd_row_num = blindsearch_database.database_blindsearch_start(obsid,
                                               pointing, "{0} {1}".format(code_comment,n))
            else:
                bsd_row_num = None
            process_vcs_wrapper(obsid, begin, end, [ra,dec], args, DI_dir,\
                                fits_dir, check, search, bsd_row_num, relaunch_script)
    return
 

#-------------------------------------------------------------------------------------------------------------
def rfifind(obsid, pointing, sub_dir, relaunch_script,
            work_dir='/group/mwaops/nswainston/blindsearch/', pbs=False,
            n_omp_threads = 20,
            bsd_row_num=None, pulsar=None,
            fits_dir=None):

    if fits_dir == None:
        fits_dir='/group/mwaops/vcs/{0}/pointings/{1}/'.format(obsid,pointing)
           
    #Calculates -numout for prepsubbands
    numout = numout_calc(fits_dir)
    
    
    rfi_batch = str(bsd_row_num) + '_rfi_{0}'.format(obsid)
    commands = []
    commands.append(add_database_function(pbs))
    commands.append("source /group/mwaops/PULSAR/psrBash.profile")
    commands.append("ncpus={0}".format(n_omp_threads))
    commands.append("export OMP_NUM_THREADS={0}".format(n_omp_threads))
    if not os.path.exists("{0}/rfi_masks/{1}_rfifind.mask".format(work_dir, obsid)):
        #if there is not already a rfi mask make one
        commands.append('run "rfifind" "-ncpus $ncpus -noclip -time 12.0 ' + '-o ' + str(obsid) +\
                        ' -zapchan 0:19,108:127,128:147,236:255,256:275,364:383,384:403,492:511,512:531,620:639,640:659,748:767,768:787,876:895,896:915,1004:1023,1024:1043,1132:1151,1152:1171,1260:1279,1280:1299,1388:1407,1408:1427,1516:1535,1536:1555,1644:1663,1664:1683,1772:1791,1792:1811,1900:1919,1920:1939,2028:2047,2048:2067,2156:2175,2176:2195,2284:2303,2304:2323,2412:2431,2432:2451,2540:2559,2560:2579,2668:2687,2688:2707,2796:2815,2816:2835,2924:2943,2944:2963,3052:3071 ' + fits_dir +\
                        '*incoh*.fits" '+ str(bsd_row_num))
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
                 submit=True, module_list=["presto/master"], export="ALL")
 
    return


#-------------------------------------------------------------------------------------------------------------
def prepdata(obsid, pointing, relaunch_script,
             work_dir='/group/mwaops/nswainston/blindsearch/', pbs=False,
             bsd_row_num=None, pulsar=None,
             n_omp_threads = 20,
             fits_dir=None, dm_max=4):
   
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
        output = subprocess.Popen(['DDplan.py','-l',str(float(dm) - 1.),'-d',str(float(dm) + 1.),'-f',str(centrefreq),'-b','30.7200067160534','-t','0.0001','-n','3072','-o','dm_temp'],stdout=subprocess.PIPE).communicate()
    else:
        output = subprocess.Popen(['DDplan.py','-l','1','-d',str(dm_max),'-f',str(centrefreq),'-b','30.7200067160534','-t','0.0001','-n','3072','-o','dm_temp'],stdout=subprocess.PIPE).communicate()
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
    numout = numout_calc(fits_dir)
    
    print dm_list    
    #Submit a bunch some prepsubbands to create our .dat files
    job_id_list = []
    dms_per_job = 1024
    for dm_line in dm_list:
        dm_start = dm_line[0]
        dm_end = float(dm_line[2]) * float(dm_line[4]) + float(dm_start)
        #DM_batch = 'DM_' + dm_start + '.batch'
        #commands = []
        #commands.append(add_database_function(pbs))
        #commands.append("source /group/mwaops/PULSAR/psrBash.profile")
        #commands.append("ncpus={0}".format(n_omp_threads))
        #commands.append("export OMP_NUM_THREADS={0}".format(n_omp_threads))
        #commands.append('run "prepsubband" "-ncpus $ncpus -lodm '+str(dm_start) +\
        #                " -dmstep " + str(dm_line[2]) + " -numdms " + str(dm_line[4])+ " -numout " +\
        #                str(numout) + " -o " + str(obsid) + " " + fits_dir + str(obsid) + \
        #                "/pointings/" + str(pointing) + "/" + str(obsid) + '*.fits" ' +str(bsd_row_num))
        
        
        while ( (dm_end - float(dm_start)) / float(dm_line[2])) > float(dms_per_job) :
            #dedisperse for only 1024 steps
            DM_batch = str(bsd_row_num) + '_DM_' + dm_start 
            prepdata_command = 'run "prepsubband" "-ncpus $ncpus -lodm ' + str(dm_start) +\
                            " -dmstep " + str(dm_line[2]) + " -numdms "+str(dms_per_job)+\
                            " -numout " + str(numout) + " -o " + str(obsid) + " " + fits_dir +\
                            '1*.fits" '+str(bsd_row_num)
            commands = []
            commands.append(add_database_function(pbs))
            commands.append("source /group/mwaops/PULSAR/psrBash.profile")
            commands.append("ncpus={0}".format(n_omp_threads))
            commands.append("export OMP_NUM_THREADS={0}".format(n_omp_threads))
            commands.append(prepdata_command)
            
            job_id = submit_slurm(DM_batch, commands,
                         batch_dir="{0}{1}/{2}/batch".format(work_dir,pointing,obsid),
                         slurm_kwargs={"time": "3:00:00", "partition": "workq"},#4 hours
                         submit=True, module_list=["presto/master"], export="ALL")
            job_id_list.append(job_id)
           
            dm_start = str(float(dm_start) + (float(dms_per_job) * float(dm_line[2])))
        steps = int((dm_end - float(dm_start)) / float(dm_line[2]))
        #last loop to get the <512 steps
        DM_batch = str(bsd_row_num) + '_DM_' + dm_start 
        prepdata_command = 'run "prepsubband" "-ncpus $ncpus -lodm ' + str(dm_start) +\
                            " -dmstep " + str(dm_line[2]) + " -numdms "+str(dms_per_job)+\
                            " -numout " + str(numout) + " -o " + str(obsid) + " " + fits_dir +\
                            '1*.fits" '+str(bsd_row_num)
        commands = []
        commands.append(add_database_function(pbs))
        commands.append("source /group/mwaops/PULSAR/psrBash.profile")
        commands.append("ncpus={0}".format(n_omp_threads))
        commands.append("export OMP_NUM_THREADS={0}".format(n_omp_threads))
        commands.append(prepdata_command)
        
        """ with mask
        batch_line = 'run "prepsubband" " -ncpus $ncpus -lodm '+str(dm_start) +\
                            " -dmstep " + str(dm_line[2]) + " -numdms " + str(steps) + " -numout " +\
                            str(numout) + " -o " + str(obsid) + " -mask " + str(obsid) +\
                            "_rfifind.mask " + fits_dir + str(obsid) + \
                            "/pointings/" + str(pointing) + "/" + str(obsid) + '*.fits" ' +str(bsd_row_num)
        """
        job_id = submit_slurm(DM_batch, commands,
                         batch_dir="{0}{1}/{2}/batch".format(work_dir,pointing,obsid),
                         slurm_kwargs={"time": "3:00:00", "partition": "workq"},#4 hours
                         submit=True, module_list=["presto/master"], export="ALL")
        job_id_list.append(job_id)    
           
    #make a job that simply restarts this program when all prepsubband jobs are complete
    print "Waiting 5 sec to make sure to the dependancy script works"
    sleep(5)
    job_id_str = ""
    for i in job_id_list:
        job_id_str += ":" + str(i)
    
    DM_depend_batch = str(bsd_row_num) + '_dep_prepsubbands'
    commands = []
    commands.append(add_database_function(pbs))
    commands.append("source /group/mwaops/PULSAR/psrBash.profile")
    commands.append("ncpus={0}".format(n_omp_threads))
    commands.append("export OMP_NUM_THREADS={0}".format(n_omp_threads))
    commands.append("cd " + work_dir + sub_dir)
    #commands.append('realfft ' + str(obsid) + '_DM0.00.dat')
    #commands.append("accelsearch -numharm 4 -zmax 0 " +str(obsid) + "_DM0.00.fft")
    commands.append('blindsearch_database.py -c prepsubband -m p -b ' +str(bsd_row_num) )
    commands.append("{0} -m s -r {1}".format(relaunch_script, bsd_row_num))
    
    submit_slurm(DM_depend_batch, commands,
                 batch_dir="{0}{1}/{2}/batch".format(work_dir,pointing,obsid),
                 slurm_kwargs={"time": "20:00", "partition": "workq"},#4 hours
                 submit=True, depend=job_id_str[1:])
    return
                
#-------------------------------------------------------------------------------------------------------------
def sort_fft(obsid, pointing, sub_dir, relaunch_script,
             work_dir='/group/mwaops/nswainston/blindsearch/', pbs=False,
             bsd_row_num=None, pulsar=None,
             fits_dir=None, dm_max=4, n_omp_threads = 20):

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
        
    os.chdir(work_dir + "/" + sub_dir)
    if pulsar==None:
        for i in range(dm_max/2):
            if i < 4:
                if not os.path.exists("DM_00" + str(i*2) + "-00" + str((i+1)*2)):
                    os.mkdir("DM_00" + str(i*2) + "-00" + str((i+1)*2))
            elif i == 4:
                if not os.path.exists("DM_008-010"):
                    os.mkdir("DM_008-010")
            elif ( 4 < i < 49 ):
                if not os.path.exists("DM_0" + str(i*2) + "-0" + str((i+1)*2)):
                    os.mkdir("DM_0" + str(i*2) + "-0" + str((i+1)*2))
            elif i == 49:    
                if not os.path.exists("DM_098-100"):
                    os.mkdir("DM_098-100")
            elif i > 49:
                if not os.path.exists("DM_" + str(i*2) + "-" + str((i+1)*2)):
                    os.mkdir("DM_" + str(i*2) + "-" + str((i+1)*2))

       
    
    DIR=work_dir + sub_dir
    length = len(DIR)
    
    all_files = os.listdir(DIR+ "/")
    
    if pulsar==None:
        dirlist = [f for f in os.listdir(DIR) if os.path.isdir(os.path.join(DIR, f))]
        for d in dirlist:
            if not d.startswith("DM"):
                dirlist.remove(d)
        for d in dirlist:
            #i and j are the start and stop dm to help with sorting
            #for weird bug where batch isn't removed from dirlist
            if d.startswith("DM_"):
                #print d
                i = int(d[3:6])
                if d.endswith("batch"):
                    j = int(d[7:-6])
                elif d.endswith("out"):
                    j = int(d[7:-4])
                else:
                    j = int(d[7:])
                    
                for f in all_files:
                    if ( f.endswith(".dat") or f.endswith(".inf") ) and not f.endswith('rfifind.inf'):
                        temp = f.split("DM")
                        if float(i) <= float(temp[1][:-4]) < float(j):
                            #print temp[2][:-4]
                            os.rename(work_dir + sub_dir + "/" + str(f),work_dir + sub_dir \
                                      + "/" + str(d) + "/" + str(f))
    else:
        dirlist = ['']
                
    dirlist.sort()
    
    #Send off jobs
    
    job_id_list =[]
    srun_commands = []
    for i, d in enumerate(dirlist):
        if d.startswith("DM_"):
            print d
            if pulsar == None:
                os.chdir(work_dir + sub_dir + "/" + d)
                #dir_files = os.listdir(work_dir + sub_dir + "/" + d + "/")
                dir_files = glob.glob(work_dir + sub_dir + "/" + d + "/*dat")
            else:
                #dir_files = os.listdir(work_dir + sub_dir + "/")
                dir_files = glob.glob(work_dir + sub_dir + "/*dat")

            fft_batch = str(bsd_row_num) + "_fft_" + d
            commands = []
            commands.append("source /group/mwaops/PULSAR/psrBash.profile")
            commands.append("ncpus={0}".format(n_omp_threads))
            commands.append("export OMP_NUM_THREADS={0}".format(n_omp_threads))
            commands.append("cd " + work_dir + sub_dir + "/" + d)
            fft_command = 'run "realfft" "'
            for fi, f in enumerate(dir_files):
                if f.endswith(".dat"):     
                    #commands.append('run "realfft" "' + str(f) + '" "'+str(bsd_row_num)+'" "'+str(i)+'"')
                    fft_command += ' ' + str(f).split("/")[-1]
                if fi%16 == 15:
                    #only 16 ffts can be done at once
                    fft_command += '" "'+str(bsd_row_num)+'" "'+str(i)+'"'
                    srun_commands.append(fft_command)
                    fft_command = 'run "realfft" "'
            fft_command += '" "'+str(bsd_row_num)+'" "'+str(i)+'"'
            srun_commands.append(fft_command)
            
            with open(fft_batch+".bash", "w") as srun_bash:
                srun_bash.write(add_temp_database_function(pbs, threads=False))
                for sc in srun_commands:
                    srun_bash.write("{}\n".format(sc))
            
            commands.append("srun -n 1 bash {}.bash".format(fft_batch))
            commands.append('blindsearch_database.py -c realfft -m m -b ' +str(bsd_row_num))
            commands.append('blindsearch_database.py -c realfft -m p -b ' +str(bsd_row_num))
            
            job_id = submit_slurm(fft_batch, commands,
                             batch_dir="{0}{1}/{2}/batch".format(work_dir,pointing,obsid),
                             slurm_kwargs={"time": "2:50:00", "partition": "workq"},#4 hours
                             submit=True, module_list=["presto/master"], export="ALL")
            job_id_list.append(job_id)

    os.chdir(work_dir + "/" + sub_dir)
    
    sleep(1)
    job_id_str = ""
    for i in job_id_list:
        job_id_str += ":" + str(i)
    
    fft_dep_batch = str(bsd_row_num) + "_dep_fft"
    commands = []
    commands.append(add_database_function(pbs))
    commands.append("source /group/mwaops/PULSAR/psrBash.profile")
    commands.append("ncpus={0}".format(n_omp_threads))
    commands.append("export OMP_NUM_THREADS={0}".format(n_omp_threads))
    
    if pulsar == None:
        relaunch_script += " -d 0"
    commands.append("{0} -m a".format(relaunch_script))
    submit_slurm(fft_dep_batch, commands,
                 batch_dir="{0}{1}/{2}/batch".format(work_dir,pointing,obsid),
                 slurm_kwargs={"time": "20:00", "partition": "workq"},#4 hours
                 submit=True, depend=job_id_str[1:])

    print "Sent off fft jobs"
    return
                
                
#-------------------------------------------------------------------------------------------------------------
def accel(obsid, pointing, sub_dir, relaunch_script,
          work_dir='/group/mwaops/nswainston/blindsearch/', pbs=False,
          bsd_row_num=None, pulsar=None, n_omp_threads = 20):
    #blindsearch_pipeline.py -o 1133329792 -p 19:45:14.00_-31:47:36.00 -m a -w /scratch2/mwaops/nswainston/tara_candidates//19:45:14.00_-31:47:36.00/1133329792/DM_058-060
    # sends off the accel search jobs
    #sub_dir = pointing + "/" + obsid + "/"
    
        
    if pulsar == None:
        DM_file_list = glob.glob(work_dir + str(sub_dir) +"/DM*/")
    else:
        DM_file_list = [work_dir + str(sub_dir)]
        dm, p = get_pulsar_dm_p(pulsar)
    for dm_i, DIR in enumerate(DM_file_list):
        os.chdir(DIR)
        dm_file = DIR[:-1].split("/")[-1]
        dir_files = glob.glob("*fft")



        #split ffts into 140 job chunks
        n = 140
        dir_files = list(chunks(dir_files,n))
        
        for dfi, dir_file in enumerate(dir_files):
            #calc processing time
            run_time = len(dir_file)*300
            run_time =  datetime.timedelta(seconds=run_time)
            print  run_time
            print dir_file
            #dir_files = ['1150234552_DM10.92.fft']
            job_id_list =[]
            

            accel_batch = str(bsd_row_num) + "_acl_" + dm_file+ "_"+str(dfi)
            commands = []
            commands.append("source /group/mwaops/PULSAR/psrBash.profile")
            commands.append("ncpus={0}".format(n_omp_threads))
            commands.append("export OMP_NUM_THREADS={0}".format(n_omp_threads))
            commands.append("cd " + DIR )
            
            with open(accel_batch+".bash", "w") as srun_bash:
                srun_bash.write(add_temp_database_function(pbs))
                srun_bash.write("ncpus={0}\n".format(n_omp_threads))
                for f in dir_file:
                    srun_bash.write('run "accelsearch" "-ncpus $ncpus -zmax 0 -flo 0.75 -fhi 500 '+\
                                    '-numharm 8 ' + f + '" "' +str(bsd_row_num) + '" "'+str(dm_i)+'"\n')

            
            commands.append("srun -n 1 -c $ncpus bash {}.bash".format(accel_batch))
            commands.append('blindsearch_database.py -c accelsearch -m m -b ' +str(bsd_row_num))
            commands.append('blindsearch_database.py -c accelsearch -m p -b ' +str(bsd_row_num) +' -d '+str(dm_i))
            
            job_id = submit_slurm(accel_batch, commands,
                                 batch_dir="{0}{1}/{2}/batch".format(work_dir,pointing,obsid),
                                 slurm_kwargs={"time": run_time, "partition": "workq"},#4 hours
                                 submit=True, module_list=["presto/master"], export="ALL")
            job_id_list.append(job_id)
    
    sleep(1)
    print "Dependancy job below"
    job_id_str = ""
    for i in job_id_list:
        job_id_str += ":" + str(i)
    
    accel_dep_batch = str(bsd_row_num) + "_dep_accel"
    commands = []
    commands.append("source /group/mwaops/PULSAR/psrBash.profile")
    commands.append("ncpus={0}".format(n_omp_threads))
    commands.append("export OMP_NUM_THREADS={0}".format(n_omp_threads))
    commands.append("{0} -m f".format(relaunch_script))

    submit_slurm(accel_dep_batch, commands,
                 batch_dir="{0}{1}/{2}/batch".format(work_dir,pointing,obsid),
                 slurm_kwargs={"time": "20:00", "partition": "workq"},#4 hours
                 submit=True, depend=job_id_str[1:])
    
    print "Sent off accel jobs"
    return
       
#-------------------------------------------------------------------------------------------------------------
def fold(obsid, pointing, sub_dir, relaunch_script,
         work_dir='/group/mwaops/nswainston/blindsearch/', pbs=False,
         bsd_row_num=None, pulsar=None, fits_dir=None, n_omp_threads = 20):
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
            
       
    #calcs sn_min for candidates
    numout = numout_calc(fits_dir)
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
    if len(cand_list) > 0:
        #sort by DM
        from operator import itemgetter
        cand_list.sort(key=itemgetter(3))
        #print cand_list
        
               
        print "Number of cands in this file: " + str(len(cand_list))
        job_id_list =[]
        
        
        #number of folds to do per sbatch job (decrease is using .fits files)        
        cands_per_batch = 30 
         

        for i,c in enumerate(cand_list):
            accel_file_name, cand_num, SN, cand_DM, period = c
            #through some stuffing around sort the fold into 100 folds per job
            #the fold option using .dat files which is quicker but inaccurate
            #fold_command = 'run "prepfold" "-ncpus $ncpus -n 128 -nsub 128 '+\
            #           "-accelcand "+c[1]+" -accelfile "+c[0]+".cand  -o " +\
            #           c[0][:-8] + " " + c[0][:-8] + '.dat" "'+str(bsd_row_num)+'" "'+str(dm_i)+'"' 
           
            #the fold options that uses .fits files which is slower but more accurate
            fold_command = 'run "prepfold" "-n 128 -nsub 128 -noclip -o {0}_{1}_{2} -p {3} -dm {4} -nosearch {5}*.fits" "{6}" "0"'.format(accel_file_name, cand_num, pointing, float(period)/1000.,cand_DM, fits_dir, bsd_row_num)
            #TODO remove this for db
              
            if (i == 0) or ((i % cands_per_batch) == 0):
                #move to correct file in batch file
                SUBDIR = work_dir + str(sub_dir)[:-1] + "/presto_profiles"
                
                fold_batch = str(bsd_row_num) + "_fold_c_{0}_{1}_{2}".format(i,i+cands_per_batch,pointing)
                commands = []
                commands.append(add_database_function(pbs))
                commands.append("source /group/mwaops/PULSAR/psrBash.profile")
                commands.append("ncpus={0}".format(n_omp_threads))
                commands.append("export OMP_NUM_THREADS={0}".format(n_omp_threads))
                commands.append("cd " + SUBDIR)
                commands.append(fold_command)
                
            else:
                commands.append(fold_command)
            if ((i+1) % cands_per_batch) == 0:
                job_id = submit_slurm(fold_batch, commands,
                             batch_dir="{0}{1}/{2}/batch".format(work_dir,pointing,obsid),
                             slurm_kwargs={"time": "4:50:00", "partition": "workq"},#4 hours
                             submit=True, module_list=["presto/master"], export="ALL")
                job_id_list.append(job_id)     
                
        if not ((len(cand_list)+1) % cands_per_batch) == 0: 
            job_id = submit_slurm(fold_batch, commands,
                             batch_dir="{0}{1}/{2}/batch".format(work_dir,pointing,obsid),
                             slurm_kwargs={"time": "4:50:00", "partition": "workq"},#4 hours
                             submit=True, module_list=["presto/master"], export="ALL")
            job_id_list.append(job_id)
        
        sleep(1)
        job_id_str = ""
        for i in job_id_list:
            job_id_str += ":" + str(i)
        
    
        
        fold_dep_batch = str(bsd_row_num) + "_dep_fold" 
        commands = []
        commands.append(add_database_function(pbs))
        commands.append("source /group/mwaops/PULSAR/psrBash.profile")
        commands.append("ncpus={0}".format(n_omp_threads))
        commands.append("export OMP_NUM_THREADS={0}".format(n_omp_threads))
        commands.append("cd " + SUBDIR )
        commands.append('echo "Searching for a profile with a sigma greater than 3"')
        """
        if pulsar == None and dm_i > 0:
            i = dm_i * 2 - 1
            commands.append('cd '+ SUBDIRpast + '\n' +\
                         'count=0\n' +\
                         'total=`ls *DM'+str(i)+'.9*.ps | wc -l`\n' +\
                         'for i in $(ls *DM'+str(i)+'.9*.ps); do\n' +\
                         'if (( $count % 100 == 0 )); then\n' +\
                         'echo "$count / $total searched"\n' +\
                         'fi\n' +\
                         'chi=`sed "13q;d" ${i%.ps}.bestprof`\n' +\
                         "if [ ${chi:20:3} -ge 3 ]; then\n" +\
                         'ps_to_png.sh ${i}\n' +\
                         'mv "${i%.ps}".png ../over_3_png/"${i%.ps}".png\n' +\
                         'echo "${i%.ps}.png is over 3"\n' +\
                         "fi\n" +\
                         "count=$(($count+1))\n" +\
                         "done")
        """
        commands.append('count=0\n' +\
                     'total=`ls *.ps | wc -l`\n' +\
                     'for i in $(ls *.ps); do\n' +\
                     'if (( $count % 100 == 0 )); then\n' +\
                     'echo "$count / $total searched"\n' +\
                     'fi\n' +\
                     'chi=`sed "13q;d" ${i%.ps}.bestprof`\n' +\
                     'chi=${chi#*=}\n'+\
                     'if [ ${chi%.*} -ge 2 ]; then\n' +\
                     'mv "${i%.ps}".png ../over_3_png/"${i%.ps}".png\n' +\
                     'echo "${i%.ps}.png is over 3"\n' +\
                     "fi\n" +\
                     "count=$(($count+1))\n" +\
                     "done")
        commands.append('blindsearch_database.py -c prepfold -m p -b ' +str(bsd_row_num) +' -d 0')
        submit_slurm(fold_dep_batch, commands,
                         batch_dir="{0}{1}/{2}/batch".format(work_dir,pointing,obsid),
                         slurm_kwargs={"time": "2:50:00", "partition": "workq"},#4 hours
                         submit=True, depend=job_id_str[1:])
    return


parser = argparse.ArgumentParser(description="""
Does a blind search for a pulsar in MWA data using the galaxy supercomputer.
""")
parser.add_argument('-o','--observation',type=str,help='The observation ID of the fits file to be searched')
parser.add_argument('-p','--pointing',type=str,help='The pointing of the fits file to be searched')
parser.add_argument('-m','--mode',type=str,help='There are three modes or steps to complete the pipeline. ["b","p","s","a","f"]. The inital mode is to beamform "b" everything using process_vcs.py. The first mode is to prepdata "p" which dedisperses the the fits files into .dat files. The second mode sort and search "s" which sorts the files it folders and performs a fft. The third mode is accel search "a" runs an accel search on them. The final mode is fold "f" which folds all possible candidates so they can be visaully inspected.', default="b")
parser.add_argument('-i','--incoh', action="store_true", help='Uses the incoh fits file location')
parser.add_argument('-w','--work_dir',type=str,help='Work directory. Default: /group/mwaops/nswainston/blindsearch/')
parser.add_argument('-s','--sub_dir',type=str,help='Used by the program to keep track of the sub directory its using')
parser.add_argument('-r','--bsd_row_num',type=int,help='Database row reference number for keeping track of the scripts.')
parser.add_argument('-d','--dm_file_int',type=int,help='Used by the program to keep track DM file being used to stagger the jobs and not send off over 9000 jobs.')
parser.add_argument('--dm_max',type=int, default = 4,help='DM max searched. Default 4')
parser.add_argument('--pulsar',type=str,help="Used to search for a known pulsar by inputing it's Jname. The code then looks within 1 DM and 15%% of the pulsar's period.")
parser.add_argument('--pbs',action="store_true",help="PBS queue mode.")
group_beamform = parser.add_argument_group('group_beamform','Beamforming Options')
group_beamform.add_argument("--DI_dir", default=None, help="Directory containing either Direction Independent Jones Matrices (as created by the RTS) or calibration_solution.bin as created by Andre Offringa's tools.[no default]")
group_beamform.add_argument('--cal_obs', '-O', type=int, help="Observation ID of calibrator you want to process.", default=None)
group_beamform.add_argument("--pulsar_file", default=None, help="Location of a file containting the pointings to be processed. Made using grid.py.")
group_beamform.add_argument("-b", "--begin", type=int, help="First GPS time to process [no default]")
group_beamform.add_argument("-e", "--end", type=int, help="Last GPS time to process [no default]")
group_beamform.add_argument("-c", "--check", type=int, help="Number of times the beamformer has attempted to redo the pointings. Stops when it gets to 5. Default 0.", default=0)
group_beamform.add_argument("-a", "--all", action="store_true",  help="Perform on entire observation span. Use instead of -b & -e.")
group_beamform.add_argument("--search", action="store_true",  help="Continue with the blindsearch pipeline after a successful beamforming check. Default False")
group_beamform.add_argument("--relaunch", action="store_true",  help="Relaunch check that doesn't send off pipeline again.")

args=parser.parse_args()

obsid = args.observation
if args.incoh:
    pointing = 'incoh'
    pointing_list = ['incoh']
else:
    pointing = args.pointing

#Default parsing
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

#check begining end times
if args.all and (args.begin or args.end):
    print "Please specify EITHER (-b,-e) OR -a"
    quit()
elif args.all:
    args.begin, args.end = meta.obs_max_min(args.observation)
    

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


#work out start and stop times for beamforming
if args.mode == "b" or args.mode == None:
    if args.pulsar_file:
        with open(args.pulsar_file) as f:
            pointing_list = f.readlines()
    elif args.pointing:
        pointing_list = [args.pointing.replace("_"," ")]
    
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
             search=args.search, bsd_row_num=args.bsd_row_num, incoh=args.incoh,
             pbs=args.pbs, pulsar=args.pulsar, check=args.check, args=args)

    
elif args.mode == "r":
    rfifind(obsid, pointing, sub_dir, relaunch_script,
            work_dir=work_dir, pbs=args.pbs, bsd_row_num=args.bsd_row_num, pulsar=args.pulsar,
            n_omp_threads=n_omp_threads)
elif args.mode == "p":
    prepdata(obsid, pointing, relaunch_script,
             work_dir=work_dir, pbs=args.pbs, bsd_row_num=args.bsd_row_num, pulsar=args.pulsar,
             n_omp_threads=n_omp_threads, fits_dir=fits_dir, dm_max=args.dm_max)
elif args.mode == "s":
    sort_fft(obsid, pointing, sub_dir, relaunch_script,
             work_dir=work_dir, pbs=args.pbs, bsd_row_num=args.bsd_row_num, pulsar=args.pulsar,
             n_omp_threads=n_omp_threads, dm_max=args.dm_max)
elif args.mode == "a":
    accel(obsid, pointing, sub_dir, relaunch_script,
          work_dir=work_dir, pbs=args.pbs, bsd_row_num=args.bsd_row_num, pulsar=args.pulsar,
          n_omp_threads=n_omp_threads)
elif args.mode == "f":
    fold(obsid, pointing, sub_dir, relaunch_script,
         work_dir=work_dir, pbs=args.pbs, bsd_row_num=args.bsd_row_num, pulsar=args.pulsar,
         n_omp_threads=n_omp_threads, fits_dir=fits_dir)

    
        
    
#blindsearch_pipeline.py -o 1166459712 -p 06:30:00.0_-28:34:00.0
