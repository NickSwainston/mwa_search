#! /usr/bin/env python

import subprocess
import os
import argparse
import urllib
import urllib2
import json
from time import sleep
import blindsearch_database

#python /group/mwaops/nswainston/bin/blindsearch_pipeline.py -o 1133329792 -p 19:45:14.00_-31:47:36.00
#python /group/mwaops/nswainston/bin/blindsearch_pipeline.py -o 1150234552 -p 00:34:08.8703_-07:21:53.409 --pulsar J0034-0721
#python /group/mwaops/nswainston/bin/blindsearch_pipeline.py -o 1099414416 -p 05:34:32_+22:00:53 --pulsar J0534+2200

#1163853320 47 tuck data


def getmeta(service='obs', params=None):
    """
    Given a JSON web service ('obs', find, or 'con') and a set of parameters as
    a Python dictionary, return the RA and Dec in degrees from the Python dictionary.
    
    getmeta(service='obs', params=None)
    """
    BASEURL = 'http://mwa-metadata01.pawsey.org.au/metadata/'
    if params:
        data = urllib.urlencode(params)  # Turn the dictionary into a string with encoded 'name=value' pairs
    else:
        data = ''
    #Validate the service name
    if service.strip().lower() in ['obs', 'find', 'con']:
        service = service.strip().lower()
    else:
        print "invalid service name: %s" % service
        return
    #Get the data
    try:
        result = json.load(urllib2.urlopen(BASEURL + service + '?' + data))
    except urllib2.HTTPError as error:
        print "HTTP error from server: code=%d, response:\n %s" % (error.code, error.read())
        return
    except urllib2.URLError as error:
        print "URL or network error: %s" % error.reason
        return
    #Return the result dictionary
    return result
    
    
def add_database_function(pbs):
    if pbs:
        batch_line ='#PBS -q gstar\n'+\
                    '#PBS -l nodes=1:ppn=6:gpus=1\n' +\
                    '#PBS -A p125_astro\n' +\
                    'function run\n' +\
                    '{\n' +\
                    '    # run command and add relevant data to the job database\n' +\
                    '    # 1st parameter is command to be run (e.g. wsclean)\n' +\
                    '    # 2nd parameter is parameters to that command (e.g. "-j $ncpus")\n' +\
                    '    # 3rd parameter is bs_id\n' +\
                    '    # 4th parameter is DM file int [optional]\n' +\
                    '    if [ -z "$4" ]; then\n' +\
                    '        rownum=`blindsearch_database.py -m s -c $1 -a "$2" -b $3 -n 1` \n' +\
                    '    else\n' +\
                    '        rownum=`blindsearch_database.py -m s -c $1 -a "$2" -b $3 -n 1 -d $4`\n' +\
                    '    fi\n' +\
                    '    $1 $2\n' +\
                    '    echo $1 $2\n' +\
                    '    errcode=$?\n' +\
                    '    blindsearch_database.py -m e -c $1 -r $rownum --errorcode $errcode\n' +\
                    '    echo "blindsearch_database.py -m e -c $1 -r $rownum --errorcode $errcode"\n' +\
                    '    if [ "$errcode" != "0" ]; then\n' +\
                    '        exit $errcode\n' +\
                    '    fi\n' +\
                    '}\n'
    else:
        batch_line ="#SBATCH --export=NONE\n" +\
                    "#SBATCH --gid=mwaops\n" +\
                    "#SBATCH --account=mwaops\n" +\
                    "#SBATCH --nodes=1\n" +\
                    'aprun="aprun -b -n 1 -d $ncpus -q "\n' +\
                    'function run\n' +\
                    '{\n' +\
                    '    # run command and add relevant data to the job database\n' +\
                    '    # 1st parameter is command to be run (e.g. wsclean)\n' +\
                    '    # 2nd parameter is parameters to that command (e.g. "-j $ncpus")\n' +\
                    '    # 3rd parameter is bs_id\n' +\
                    '    # 4th parameter is DM file int [optional]\n' +\
                    '    if [ -z "$4" ]; then\n' +\
                    '        rownum=`blindsearch_database.py -m s -c $1 -a "$2" -b $3 -n $ncpus` \n' +\
                    '    else\n' +\
                    '        rownum=`blindsearch_database.py -m s -c $1 -a "$2" -b $3 -n $ncpus -d $4`\n' +\
                    '    fi\n' +\
                    '    $aprun $1 $2\n' +\
                    '    echo $aprun $1 $2\n' +\
                    '    errcode=$?\n' +\
                    '    blindsearch_database.py -m e -c $1 -r $rownum --errorcode $errcode\n' +\
                    '    echo "blindsearch_database.py -m e -c $1 -r $rownum --errorcode $errcode"\n' +\
                    '    if [ "$errcode" != "0" ]; then\n' +\
                    '        exit $errcode\n' +\
                    '    fi\n' +\
                    '}\n'
    return batch_line
    
    
def add_temp_database_function(pbs):
    if pbs:
        batch_line ='#PBS -q gstar\n'+\
                    '#PBS -l nodes=1:ppn=6:gpus=1\n' +\
                    '#PBS -A p125_astro\n' +\
                    'function run\n' +\
                    '{\n' +\
                    '    # run command and add relevant data to the job database\n' +\
                    '    # 1st parameter is command to be run (e.g. wsclean)\n' +\
                    '    # 2nd parameter is parameters to that command (e.g. "-j $ncpus")\n' +\
                    '    # 3rd parameter is bs_id\n' +\
                    '    # 4th parameter is DM file int [optional]\n' +\
                    '    if [ -z "$4" ]; then\n' +\
                    '        echo `date +%Y-%m-%d" "%H:%M:%S`",$1,$2,$3,1" >> ${1}_temp_database_file.csv\n' +\
                    '    else\n' +\
                    '        echo `date +%Y-%m-%d" "%H:%M:%S`",$1,$2,$3,1,$4" >> ${1}_temp_database_file.csv\n' +\
                    '    fi\n' +\
                    '    $1 $2\n' +\
                    '    echo $1 $2\n' +\
                    '    errcode=$?\n' +\
                    '    echo `date +%Y-%m-%d" "%H:%M:%S`",$errcode" >> ${1}_temp_database_file.csv\n' +\
                    '    if [ "$errcode" != "0" ]; then\n' +\
                    '        exit $errcode\n' +\
                    '    fi\n' +\
                    '}\n'
    else:
        batch_line ="#SBATCH --export=NONE\n" +\
                    "#SBATCH --gid=mwaops\n" +\
                    "#SBATCH --account=mwaops\n" +\
                    "#SBATCH --nodes=1\n" +\
                    'aprun="aprun -b -n 1 -d $ncpus -q "\n' +\
                    'function run\n' +\
                    '{\n' +\
                    '    # run command and add relevant data to the job database\n' +\
                    '    # 1st parameter is command to be run (e.g. wsclean)\n' +\
                    '    # 2nd parameter is parameters to that command (e.g. "-j $ncpus")\n' +\
                    '    # 3rd parameter is bs_id\n' +\
                    '    # 4th parameter is DM file int [optional]\n' +\
                    '    if [ -z "$4" ]; then\n' +\
                    '        echo `date +%Y-%m-%d" "%H:%M:%S`",$1,$2,$3,$ncpus" >> ${1}_temp_database_file.csv\n' +\
                    '    else\n' +\
                    '        echo `date +%Y-%m-%d" "%H:%M:%S`",$1,$2,$3,$ncpus,$4" >> ${1}_temp_database_file.csv\n' +\
                    '    fi\n' +\
                    '    $aprun $1 $2\n' +\
                    '    echo $aprun $1 $2\n' +\
                    '    errcode=$?\n' +\
                    '    echo `date +%Y-%m-%d" "%H:%M:%S`",$errcode" >> ${1}_temp_database_file.csv\n' +\
                    '    if [ "$errcode" != "0" ]; then\n' +\
                    '        exit $errcode\n' +\
                    '    fi\n' +\
                    '}\n'
    return batch_line
    #'echo "FFT,'+str(bs_id)+',realfft,' + str(f) + ',$ncpus,'+str(i)+',date +%Y-%m-%d" "%H:%M:%S" >> temp_database_file.csv\n' +\
    
    
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
    
#-------------------------------------------------------------------------------------------------------------
def rfifind(obsid, pointing, work_dir, sub_dir,pbs,pulsar=None):
    code_comment = raw_input("Please right a comment describing the purpose of this blindsearch. eg testing: ")
    bs_id = blindsearch_database.database_blindsearch_start(obsid, pointing, code_comment)
    
    #Set up some directories and move to it
    if not os.path.exists(work_dir + pointing):
            os.mkdir(work_dir + pointing)
    if not os.path.exists(work_dir + pointing + "/" + obsid): 
            os.mkdir(work_dir + pointing + "/" + obsid)
    if not pulsar == None:
        if not os.path.exists(work_dir + pointing + "/" + obsid + "/" + pulsar): 
            os.mkdir(work_dir + pointing + "/" + obsid + "/" + pulsar)
        os.chdir(work_dir + pointing + "/" + obsid + "/" + pulsar)
        sub_dir = pointing + "/" + obsid + "/" + pulsar + "/"
    else:
        os.chdir(work_dir + pointing + "/" + obsid)
        sub_dir = pointing + "/" + obsid + "/"
    
    if not os.path.exists(work_dir + sub_dir + "/batch"): 
            os.mkdir(work_dir + sub_dir+ "/batch")
    if not os.path.exists(work_dir + sub_dir + "/out"): 
            os.mkdir(work_dir + sub_dir + "/out")
    
    if pbs:
        fits_dir = '/lustre/projects/p125_astro/DATA/'
        qsub = 'qsub '
        ncpuscom = ''
    else:
        fits_dir = '/group/mwaops/vcs/'
        qsub = 'sbatch '
        ncpuscom = '-ncpus $ncpus '
            
    #Calculates -numout for prepsubbands
    numout = numout_calc(fits_dir + str(obsid) + "/pointings/" + str(pointing) + "/")
            
    #send off rfi job
    with open('batch/rfifind.batch','w') as batch_file:
        if pbs:
            batch_line = "#!/bin/bash -l\n" +\
                         "#PBS -N rfifind\n" +\
                         "#PBS -o out/rfifind.out\n" +\
                         "#PBS -e out/rfifind.error\n" +\
                         "#PBS -l walltime=3:50:00\n" +\
                         "cd " + work_dir + sub_dir + "\n"
        else:
            batch_line = "#!/bin/bash -l\n" +\
                         "#SBATCH --partition=workq\n" +\
                         "#SBATCH --job-name=rfifind\n" +\
                         "#SBATCH --output=out/rfifind.out\n" +\
                         "#SBATCH --time=3:50:00\n" +\
                         "ncpus=20\n" +\
                         "export OMP_NUM_THREADS=$ncpus\n"
        batch_file.write(batch_line)
        batch_file.write(add_database_function(pbs))
        batch_line = 'run "rfifind" "' + ncpuscom + '-noclip -time 12.0 '+\
                        '-o ' + str(obsid) + ' -zapchan 0:19,108:127,128:147,236:255,256:275,364:383,384:403,492:511,512:531,620:639,640:659,748:767,768:787,876:895,896:915,1004:1023,1024:1043,1132:1151,1152:1171,1260:1279,1280:1299,1388:1407,1408:1427,1516:1535,1536:1555,1644:1663,1664:1683,1772:1791,1792:1811,1900:1919,1920:1939,2028:2047,2048:2067,2156:2175,2176:2195,2284:2303,2304:2323,2412:2431,2432:2451,2540:2559,2560:2579,2668:2687,2688:2707,2796:2815,2816:2835,2924:2943,2944:2963,3052:3071 ' + fits_dir + str(obsid) + \
                        '/pointings/' + str(pointing) + '/' + str(obsid) + '*.fits" '+str(bs_id)+"\n"+\
                        'blindsearch_database.py -c rfifind -m p -b ' +str(bs_id) + '\n'+\
                        "blindsearch_pipeline.py -o "\
                          + str(obsid) + " -p " + str(pointing) + " -m p -w " + work_dir +\
                          " -s " +str(sub_dir)+ " -r " +str(bs_id)
        batch_file.write(batch_line)
        if not pulsar == None:
            batch_line = " --pulsar " + str(pulsar)
            batch_file.write(batch_line)
        if pbs:
            batch_file.write(" --pbs ")
        batch_file.write("\n")
        if pbs:
            batch_line ="prepdata " + ncpuscom + " -dm 0 " +\
                        "-numout " + str(numout) + " -o " + str(obsid) + \
                        "_DM0.00 " + fits_dir + str(obsid) + \
                        "/pointings/" + str(pointing) + "/" + str(obsid) + "*.fits"
        else:
            batch_line ="$aprun prepdata " + ncpuscom + " -dm 0 " +\
                        "-numout " + str(numout) + " -o " + str(obsid) + \
                        "_DM0.00 " + fits_dir + str(obsid) + \
                        "/pointings/" + str(pointing) + "/" + str(obsid) + "*.fits"
        batch_file.write(batch_line)
    submit_line = qsub + 'batch/rfifind.batch'
    submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
    for line in submit_cmd.stdout:
        print line,
    return


#-------------------------------------------------------------------------------------------------------------
def prepdata(obsid, pointing, work_dir, sub_dir, bs_id,pbs,pulsar=None):
    if not pulsar == None:
        os.chdir(work_dir + pointing + "/" + obsid + "/" + pulsar)
        sub_dir = pointing + "/" + obsid + "/" + pulsar + "/"
    else:
        os.chdir(work_dir + pointing + "/" + obsid)
        sub_dir = pointing + "/" + obsid + "/"
    
    #Get the centre freq channel and then run DDplan.py to work out the most effective DMs
    print "Obtaining metadata from http://mwa-metadata01.pawsey.org.au/metadata/ for OBS ID: " + str(obsid)
    beam_meta_data = getmeta(service='obs', params={'obs_id':obsid})
    channels = beam_meta_data[u'rfstreams'][u"0"][u'frequencies']
    minfreq = float(min(channels))
    maxfreq = float(max(channels))
    centrefreq = 1.28 * (minfreq + (maxfreq-minfreq)/2) #in MHz
    
    """
    if not pulsar == None:
        dm, p = get_pulsar_dm_p(pulsar)
        output = subprocess.Popen(['DDplan.py','-l',str(float(dm) - 1.),'-d',str(float(dm) + 1.),'-f',str(centrefreq),'-b','30.7200067160534','-t','0.0001','-n','3072'],stdout=subprocess.PIPE).communicate()
    else:
        output = subprocess.Popen(['DDplan.py','-l','0','-d','300','-f',str(centrefreq),'-b','30.7200067160534','-t','0.0001','-n','3072'],stdout=subprocess.PIPE).communicate()
    subprocess.check_call("\n", shell=True)
    dm_list = []
    print output[0]
    lines = output[0].split('\n')
    for l in lines[13:-4]: 
        columns = l.split()
        dm_list.append(columns)
    """
    dm_list = [['1.000','100.000','0.20','1','245','1']]
    if pbs:
        fits_dir = '/lustre/projects/p125_astro/DATA/'
        qsub = 'qsub '
        ncpuscom = ''
    else:
        fits_dir = '/group/mwaops/vcs/'
        qsub = 'sbatch '
        ncpuscom = '-ncpus $ncpus '
           
    #Calculates -numout for prepsubbands
    numout = numout_calc(fits_dir + str(obsid) + "/pointings/" + str(pointing) + "/")
    
    
    #Submit a bunch some prepsubbands to create our .dat files
    job_id_list = []
    for dm_line in dm_list:
        dm_start = dm_line[0]
        dm_end = float(dm_line[2]) * float(dm_line[4]) + float(dm_start)
        while ( (dm_end - float(dm_start)) / float(dm_line[2])) > 512. :
            #dedisperse for only 512 steps
            with open('batch/DM_' + dm_start + '.batch','w') as batch_file:
                if pbs:
                    batch_line = "#!/bin/bash -l\n" +\
                                 "#PBS -N prepsub_" + dm_start + "\n" +\
                                 "#PBS -o out/DM_" + dm_start + ".out\n" +\
                                 "#PBS -e out/DM_" + dm_start + ".error\n" +\
                                 "#PBS -l walltime=3:50:00\n" +\
                                 "cd " + work_dir + sub_dir + "\n"
                    batch_file.write(batch_line)
                else:
                    batch_line = "#!/bin/bash -l\n" +\
                                 "#SBATCH --partition=workq\n" +\
                                 "#SBATCH --job-name=prepsub_" + dm_start + "\n" +\
                                 "#SBATCH --output=out/DM_" + dm_start + ".out\n" +\
                                 "#SBATCH --time=3:50:00\n"+\
                                 "ncpus=20\n" +\
                                 "export OMP_NUM_THREADS=$ncpus\n"
                    batch_file.write(batch_line)
                batch_file.write(add_database_function(pbs))
                batch_line = 'run "prepsubband" "'+ncpuscom + '-lodm ' +\
                                str(dm_start) + " -dmstep " + str(dm_line[2]) + " -numdms 512 -numout " +\
                                str(numout) + " -o " + str(obsid) + " -mask " + str(obsid) +\
                                "_rfifind.mask " + fits_dir + str(obsid) + \
                                "/pointings/" + str(pointing) + "/" + str(obsid) + '*.fits" '+str(bs_id)
                batch_file.write(batch_line)
            submit_line = qsub + 'batch/DM_' + dm_start + '.batch'
            submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
            for line in submit_cmd.stdout:
                print line,
                print 'batch/DM_' + dm_start + '.batch'
                if pbs:
                    if '.pbs.hpc.swin.edu.au' in line:
                        job_id_list.append(line[:7])
                else:
                    if "Submitted" in line:
                        (word1,word2,word3,jobid) = line.split()
                        job_id_list.append(jobid)
            
            dm_start = str(float(dm_start) + (512. * float(dm_line[2])))
        steps = int((dm_end - float(dm_start)) / float(dm_line[2]))
        #last loop to get the <512 steps
        with open('batch/DM_' + dm_start + '.batch','w') as batch_file:
            if pbs:
                batch_line = "#!/bin/bash -l\n" +\
                             "#PBS -N prepsub_" + dm_start + "\n" +\
                             "#PBS -o out/DM_" + dm_start + ".out\n" +\
                             "#PBS -e out/DM_" + dm_start + ".error\n" +\
                             "#PBS -l walltime=3:50:00\n"+\
                             "cd " + work_dir + sub_dir + "\n"
                batch_file.write(batch_line)
            else:
                batch_line = "#!/bin/bash -l\n" +\
                             "#SBATCH --partition=workq\n" +\
                             "#SBATCH --job-name=prepsub_" + dm_start + "\n" +\
                             "#SBATCH --output=out/DM_" + dm_start + ".out\n" +\
                             "#SBATCH --time=3:50:00\n" +\
                             "ncpus=20\n"+\
                             "export OMP_NUM_THREADS=$ncpus\n"
                batch_file.write(batch_line)
            batch_file.write(add_database_function(pbs))
            batch_line = 'run "prepsubband" "'+ncpuscom + '-lodm '+str(dm_start) +\
                                " -dmstep " + str(dm_line[2]) + " -numdms " + str(steps) + " -numout " +\
                                str(numout) + " -o " + str(obsid) + " " + fits_dir + str(obsid) + \
                                "/pointings/" + str(pointing) + "/" + str(obsid) + '*.fits" ' +str(bs_id)
            """ with mask
            batch_line = 'run "prepsubband" "'+ncpuscom + '-lodm '+str(dm_start) +\
                                " -dmstep " + str(dm_line[2]) + " -numdms " + str(steps) + " -numout " +\
                                str(numout) + " -o " + str(obsid) + " -mask " + str(obsid) +\
                                "_rfifind.mask " + fits_dir + str(obsid) + \
                                "/pointings/" + str(pointing) + "/" + str(obsid) + '*.fits" ' +str(bs_id)
            """
            batch_file.write(batch_line)
            
        submit_line = qsub + 'batch/DM_' + dm_start + '.batch'
        submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
        for line in submit_cmd.stdout:
            print line,
            if pbs:
                if '.pbs.hpc.swin.edu.au' in line:
                    job_id_list.append(line[:7])
            else:
                if "Submitted" in line:
                    (word1,word2,word3,jobid) = line.split()
                    job_id_list.append(jobid)
    
    
    #make a job that simply restarts this program when all prepsubband jobs are complete
    print "Waiting 5 sec to make sure to the dependancy script works"
    sleep(5)
    job_id_str = ""
    for i in job_id_list:
        job_id_str += ":" + str(i)
    with open('batch/dependancy_prepsubbands.batch','w') as batch_file:
        if pbs:#TODO fix
            batch_line = "#!/bin/bash -l\n" +\
                         "#PBS -N dependancy\n" +\
                         "#PBS -o out/dependancy_prepsubbands.out\n" +\
                         "#PBS -e out/dependancy_prepsubbands.error\n" +\
                         "#PBS -l walltime=3:50:00\n"+\
                         '#PBS -q gstar\n'+\
                         '#PBS -l nodes=1:ppn=6:gpus=1\n'+\
                         '#PBS -A p125_astro\n' +\
                         "#PBS -W depend=afterany" + job_id_str + "\n" +\
                         "cd " + work_dir + sub_dir + "\n"+\
                         'realfft ' + str(obsid) + '_DM0.00.dat\n'+\
                         "accelsearch -numharm 4 -zmax 0 " +str(obsid) + "_DM0.00.fft\n"+\
                         'blindsearch_database.py -c prepsubband -m p -b ' +str(bs_id) + '\n'+\
                         "blindsearch_pipeline.py -o "\
                              + str(obsid) + " -p " + str(pointing) + " -m s -w " + work_dir +\
                              " -s " +str(sub_dir) + " -r " + str(bs_id)
            batch_file.write(batch_line)
        else:
            batch_line = "#!/bin/bash -l\n" +\
                         "#SBATCH --job-name=dependancy\n" +\
                         "#SBATCH --output=out/dependancy_prepsubbands.out\n" +\
                         "#SBATCH --export=NONE\n" +\
                         "#SBATCH --partition=workq\n" +\
                         "#SBATCH --time=0:05:00\n" +\
                         "#SBATCH --gid=mwaops\n" +\
                         "#SBATCH --account=mwaops\n" +\
                         "#SBATCH --nodes=1\n" +\
                         "#SBATCH --dependency=afterany" + job_id_str + "\n" +\
                         "export OMP_NUM_THREADS=8\n"+\
                         'aprun -b -n 1 -d 8 -q realfft ' + str(obsid) + '_DM0.00.dat\n'+\
                         "accelsearch -numharm 4 -zmax 0 " +str(obsid) + "_DM0.00.fft\n"+\
                         'blindsearch_database.py -c prepsubband -m p -b ' +str(bs_id) + '\n'+\
                         "blindsearch_pipeline.py -o "\
                              + str(obsid) + " -p " + str(pointing) + " -m s -w " + work_dir +\
                              " -s " +str(sub_dir) + " -r " + str(bs_id)
            batch_file.write(batch_line)
        if not pulsar == None:
            batch_line = " --pulsar " + pulsar
            batch_file.write(batch_line)
        
            
    submit_line = qsub + 'batch/dependancy_prepsubbands.batch'
    submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
    for line in submit_cmd.stdout:
            print line,
    print "Sent off prepsubband jobs"
    return
                
#-------------------------------------------------------------------------------------------------------------
def sort_fft(obsid, pointing, work_dir, sub_dir, bs_id,pbs, pulsar=None):
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
        for i in range(4):
            if not os.path.exists("DM_00" + str(i*2) + "-00" + str((i+1)*2)):
                os.mkdir("DM_00" + str(i*2) + "-00" + str((i+1)*2))
        if not os.path.exists("DM_008-010"):
            os.mkdir("DM_008-010")
        for i in range(5,49):
            if not os.path.exists("DM_0" + str(i*2) + "-0" + str((i+1)*2)):
                os.mkdir("DM_0" + str(i*2) + "-0" + str((i+1)*2))
        if not os.path.exists("DM_098-100"):
            os.mkdir("DM_098-100")
        """
        for i in range(50,150):
            if not os.path.exists("DM_" + str(i*2) + "-" + str((i+1)*2)):
                os.mkdir("DM_" + str(i*2) + "-" + str((i+1)*2))
        """

    if pbs:
        fits_dir = '/lustre/projects/p125_astro/DATA/'
        qsub = 'qsub '
        ncpuscom = ''
    else:
        fits_dir = '/group/mwaops/vcs/'
        qsub = 'sbatch '
        ncpuscom = '-ncpus $ncpus '
        
    
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
            if d == 'batch' or d == 'out' or d == 'other_png':
                i = None
            else:
                print d
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
    for i, d in enumerate(dirlist):
        if pulsar == None:
            if not os.path.exists(work_dir + sub_dir + "/" + d + "/batch"):
                os.mkdir(work_dir + sub_dir + "/" + d + "/batch")
            if not os.path.exists(work_dir + sub_dir + "/" + d + "/out"):
                os.mkdir(work_dir + sub_dir + "/" + d + "/out")
            os.chdir(work_dir + sub_dir + "/" + d)
            dir_files = os.listdir(work_dir + sub_dir + "/" + d + "/")
        else:
            dir_files = os.listdir(work_dir + sub_dir + "/")
        with open('batch/fft' + d + '.batch','w') as batch_file:
            if pbs:
                batch_line = "#!/bin/bash -l\n" +\
                             "#PBS -N fft_" + d + "\n" +\
                             "#PBS -o " + work_dir + sub_dir +\
                                      "/" + d + "/out/fft_" + d + ".out\n" +\
                             "#PBS -e " + work_dir + sub_dir +\
                                      "/" + d + "/out/fft_" + d + ".error\n" +\
                             "#PBS -l walltime=4:50:00\n"+\
                             "cd " + work_dir + sub_dir + "/" + d + "\n"
                batch_file.write(batch_line)
            else:
                batch_line = "#!/bin/bash -l\n" +\
                             "#SBATCH --partition=gpuq\n" +\
                             "#SBATCH --job-name=fft_" + d + "\n" +\
                             "#SBATCH --output=" + work_dir + sub_dir +\
                                               "/" + d + "/out/fft_" + d + ".out\n" +\
                             "#SBATCH --time=4:50:00\n" +\
                             "ncpus=8\n" +\
                             "export OMP_NUM_THREADS=$ncpus\n"
                batch_file.write(batch_line)   
            batch_file.write(add_temp_database_function(pbs))         
            for f in dir_files:
                if f.endswith(".dat"):     
                    batch_line = 'run "realfft" "' + str(f) + '" "'+str(bs_id)+'" "'+str(i)+'"\n'
                    #batch_line = 'run realfft "' + str(f) + '" ' + work_dir + ' blindsearch ' + obsid +'\n'
                    batch_file.write(batch_line)
            
            batch_line = 'blindsearch_database.py -c realfft -m m -b ' +str(bs_id) + '\n'
            batch_file.write(batch_line)
        submit_line = qsub + 'batch/fft' + d + '.batch'
        submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
        for line in submit_cmd.stdout:
                print line,
                if pbs:
                    if '.pbs.hpc.swin.edu.au' in line:
                        job_id_list.append(line[:7])
                else:
                    if "Submitted" in line:
                        (word1,word2,word3,jobid) = line.split()
                        job_id_list.append(jobid)
                #print submit_cmd.communicate()[0]
    
    os.chdir(work_dir + "/" + sub_dir)
    
    sleep(1)
    job_id_str = ""
    for i in job_id_list:
        job_id_str += ":" + str(i)
    with open('batch/dependancy_fft' + d + '.batch','w') as batch_file:
        if pbs:#TODO fix
            batch_line = "#!/bin/bash -l\n" +\
                         "#PBS -N dependancy_fft\n" +\
                         "#PBS -o out/dependancy_fft.out\n" +\
                         "#PBS -e out/dependancy_fft.error\n" +\
                         "#PBS -l walltime=0:05:00\n"+\
                         '#PBS -q gstar\n'+\
                         '#PBS -l nodes=1:ppn=6:gpus=1\n'+\
                         '#PBS -A p125_astro\n' +\
                         "#PBS -W depend=afterany" + job_id_str + "\n" +\
                         "cd " + work_dir + sub_dir + "\n"
        else:
            batch_line = "#!/bin/bash -l\n" +\
                         "#SBATCH --job-name=dependancy_fft\n" +\
                         "#SBATCH --output=out/dependancy_fft.out\n" +\
                         "#SBATCH --export=NONE\n" +\
                         "#SBATCH --partition=workq\n" +\
                         "#SBATCH --time=0:05:00\n" +\
                         "#SBATCH --gid=mwaops\n" +\
                         "#SBATCH --account=mwaops\n" +\
                         "#SBATCH --nodes=1\n" +\
                         "#SBATCH --dependency=afterany" + job_id_str + "\n" 
        batch_file.write(batch_line)
        batch_line = 'blindsearch_database.py -c realfft -m p -b ' +str(bs_id) + '\n'+\
                     "blindsearch_pipeline.py -o " +\
                                 str(obsid) + " -p " + str(pointing) + " -m a -w " + work_dir +\
                                 " -s " + str(sub_dir) + ' -r ' + str(bs_id)
        batch_file.write(batch_line) 
        if not pulsar == None:
            batch_line = " --pulsar " + pulsar
            batch_file.write(batch_line)
        else:
            batch_line = " -d 0"
            batch_file.write(batch_line)
        if pbs:
            batch_file.write(" --pbs ")
    submit_line = qsub + 'batch/dependancy_fft' + d + '.batch'
    submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
    #print submit_cmd.stdout
    print "Sent off fft jobs"
    return
                
                
#-------------------------------------------------------------------------------------------------------------
def accel(obsid, pointing, work_dir, sub_dir, dm_i, bs_id,pbs, pulsar=None):
    #blindsearch_pipeline.py -o 1133329792 -p 19:45:14.00_-31:47:36.00 -m a -w /scratch2/mwaops/nswainston/tara_candidates//19:45:14.00_-31:47:36.00/1133329792/DM_058-060
    # sends off the accel search jobs
    #sub_dir = pointing + "/" + obsid + "/"
    
    if pbs:
        fits_dir = '/lustre/projects/p125_astro/DATA/'
        qsub = 'qsub '
        ncpuscom = ''
    else:
        fits_dir = '/group/mwaops/vcs/'
        qsub = 'sbatch '
        ncpuscom = '-ncpus $ncpus '
    
    print dm_i
    dm_file = dm_i_to_file(dm_i)
    
    if pulsar == None:
        DIR = work_dir + str(sub_dir) + dm_file
    else:
        DIR = work_dir + str(sub_dir)
        dm, p = get_pulsar_dm_p(pulsar)
    os.chdir(DIR)
    dir_files = os.listdir(DIR)
    #dir_files = ['1150234552_DM10.92.fft']
    job_id_list =[]
    with open('batch/accel_' + dm_file + '.batch','w') as batch_file:
        if pbs:
            batch_line = "#!/bin/bash -l\n" +\
                         "#PBS -N accel_" + dm_file + "\n" +\
                         "#PBS -o " + DIR + "/out/accel_" + dm_file + ".out\n" +\
                         "#PBS -e " + DIR + "/out/accel_" + dm_file + ".error\n" +\
                         "#PBS -l walltime=2:50:00\n"+\
                         "cd " + DIR + "\n"
            batch_file.write(batch_line)
        else:
            batch_line = "#!/bin/bash -l\n" +\
                         "#SBATCH --partition=workq\n" +\
                         "#SBATCH --job-name=accel_" + dm_file + "\n" +\
                         "#SBATCH --output=" + DIR + "/out/accel_" + dm_file + ".out\n" +\
                         "#SBATCH --time=1:00:00\n" +\
                         "ncpus=20\n"+\
                         "export OMP_NUM_THREADS=$ncpus\n"
            batch_file.write(batch_line)
        batch_file.write(add_temp_database_function(pbs))
        
        for f in dir_files:
            if f.endswith(".fft"):
                #if pulsar == None:
                batch_line = 'run "accelsearch" "'  + ncpuscom + ' -flo 1 -fhi 500 '+\
                                '-numharm 8 -zmax 0 ' + f + '" "' +str(bs_id) + '" "'+str(dm_i)+'"\n'
                                # + ' ' + DIR[:len(work_dir)] + ' ' + DIR[len(work_dir):] + ' ' + obsid 
                batch_file.write(batch_line)
                """
                else:
                
                batch_line = 'run "accelsearch" "-ncpus $ncpus -flo ' +\
                                     str(1./(float(p)*1.15)) + ' -fhi ' +\
                                     str(1./(float(p)*0.85)) + ' -numharm 8 -zmax 0 ' + f +\
                                     '" "' +str(bs_id) + '" "'+str(dm_i)+'"\n'
                                     #+ '" ' +\
                                     #DIR[:len(work_dir)] + ' ' + DIR[len(work_dir):] + ' ' + obsid 
                batch_file.write(batch_line)
                """
            
    submit_line = qsub + 'batch/accel_' + dm_file + '.batch'
    submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
    for line in submit_cmd.stdout:
        print line,
        if pbs:
            if '.pbs.hpc.swin.edu.au' in line:
                job_id_list.append(line[:7])
        else:
            if "Submitted" in line:
                (word1,word2,word3,jobid) = line.split()
                job_id_list.append(jobid)
                #print submit_cmd.communicate()[0]
    sleep(1)
    print "Dependancy job below"
    job_id_str = ""
    for i in job_id_list:
        job_id_str += ":" + str(i)
    with open('batch/dependancy_accel.batch','w') as batch_file:
        if pbs:#TODO fix
            batch_line = "#!/bin/bash -l\n" +\
                         "#PBS -N dependancy_accel\n" +\
                         "#PBS -o out/dependancy_accel.out\n" +\
                         "#PBS -e out/dependancy_accel.error\n" +\
                         "#PBS -l walltime=0:05:00\n"+\
                         '#PBS -q gstar\n'+\
                         '#PBS -l nodes=1:ppn=6:gpus=1\n'+\
                         '#PBS -A p125_astro\n' +\
                         "#PBS -W depend=afterany" + job_id_str + "\n" +\
                         "cd " + DIR + "\n"
        else:                 
            batch_line = "#!/bin/bash -l\n" +\
                         "#SBATCH --job-name=dependancy_accel\n" +\
                         "#SBATCH --output=out/dependancy_accel.out\n" +\
                         "#SBATCH --export=NONE\n" +\
                         "#SBATCH --partition=workq\n" +\
                         "#SBATCH --time=0:05:00\n" +\
                         "#SBATCH --gid=mwaops\n" +\
                         "#SBATCH --account=mwaops\n" +\
                         "#SBATCH --nodes=1\n" +\
                         "#SBATCH --dependency=afterany" + job_id_str + "\n" 
        batch_file.write(batch_line) 
        batch_line = 'blindsearch_database.py -c accelsearch -m m -b ' +str(bs_id) + '\n' +\
                     'blindsearch_database.py -c accelsearch -m p -b ' +str(bs_id) +' -d '+str(dm_i)+'\n'+\
                     "blindsearch_pipeline.py -o " + str(obsid) +\
                             " -p " + str(pointing) + " -m f -w " + work_dir + " -s " + str(sub_dir) +\
                             " -r " + str(bs_id)
        batch_file.write(batch_line)  
        if pulsar == None:
             batch_line = " -d " + str(dm_i)
             batch_file.write(batch_line)
        else:
            batch_line = " --pulsar " + str(pulsar)
            batch_file.write(batch_line)
        if pbs:
            batch_file.write(" --pbs ")
    submit_line = qsub + 'batch/dependancy_accel.batch'
    submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
    for line in submit_cmd.stdout:
                print line,
    print "Sent off accel jobs"
    return
       
#-------------------------------------------------------------------------------------------------------------
def fold(obsid, pointing, work_dir, sub_dir, dm_i, bs_id,pbs, pulsar = None):
    from math import floor
    dm_file_orig = dm_i_to_file(dm_i)
    
    DIR=work_dir + str(sub_dir)
    os.chdir(DIR)
    
    #run accel_sift.py to find candidates with DM structure
    if pulsar == None: 
        if pbs:
            submit_line = 'python ~/My-Scripts/ACCEL_sift.py ' + dm_file_orig
            file_loc = 'cand_files/cands_'+dm_file_orig+'.txt'
        else:
            submit_line = 'python /group/mwaops/nswainston/bin/ACCEL_sift.py ' + dm_file_orig
            file_loc = 'cand_files/cands_'+dm_file_orig+'.txt'
    else:
        if pbs:
            submit_line = 'python ~/My-Scripts/ACCEL_sift.py .'
            file_loc = 'ACCEL_sift_cands.txt'
        else:
            submit_line = 'python /group/mwaops/nswainston/bin/ACCEL_sift.py .'
            file_loc = 'ACCEL_sift_cands.txt'
            
    if pbs:
        fits_dir = '/lustre/projects/p125_astro/DATA/'
        qsub = 'qsub '
        ncpuscom = ''
    else:
        fits_dir = '/group/mwaops/vcs/'
        qsub = 'sbatch '
        ncpuscom = '-ncpus $ncpus '
        
    #calcs sn_min for candidates
    numout = numout_calc(fits_dir + str(obsid) + "/pointings/" + str(pointing) + "/")
    from math import sqrt,log
    sn_min = ( sqrt(log(numout)) - 0.88 ) / 0.47
    
    cand_list = [] 
    #TODO temporary fix because ACCEL_Sift doesn't work on gstar    
    if pbs:
        #don't do accelsift
        all_files = os.listdir(DIR+ "/" + dm_file_orig + "/")
        for f in all_files:
            if f.endswith("ACCEL_0"):
                with open(DIR+ "/" + dm_file_orig + "/" +f,'rb') as accel_cand_file:
                    lines = accel_cand_file.readlines()
                    for l in lines[3:]:
                        l = l.split()
                        if len(l) == 0:
                            break
                        if float(l[1]) > sn_min:
                            cand_list.append([f,l[0],l[1],f[13:-8],l[5]])
    else:
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
                            cand_list.append([cand_line[0].split(':')[0],cand_line[0].split(':')[1],cand_line[2],cand_line[1],cand_line[7]])
                        #print cand_line[2]
        #print len(cand_list)
    
    
    
    
    
    if len(cand_list) > 0:
        #sort by DM
        from operator import itemgetter
        cand_list.sort(key=itemgetter(3))
        #print cand_list
        
        #checks you're in the right file (Will change to SUBDIR in the batch file only)
        if pulsar == None:
            dm_i_temp = floor( floor(float(cand_list[0][3]))/2 )
            dm_file = dm_i_to_file(dm_i_temp)    
            SUBDIRpast = work_dir + str(sub_dir) + dm_file
            SUBDIRorig = work_dir + str(sub_dir) + dm_file_orig #current
            #if DM_000-002 current = past
        else:
            SUBDIR = work_dir + str(sub_dir)[:-1]
        
        
        """ #simple case if you want to send off one fold a job  
        for c in cand_list:
            #TODO add a check that you're not folding the same thing again
            with open('batch/fold_' + dm_file_orig + '_' + str(fold_num) + '.batch','w') as batch_file:
                batch_line = "#!/bin/bash -l\n" +\
                             "#SBATCH --partition=gpuq\n" +\
                             "#SBATCH --job-name=fold_" + c[0] + '_' + c[2] + "\n" +\
                             "#SBATCH --output=" + DIR + "/out/fold_" + c[0] + '_' + c[2] +".out\n" +\
                             "#SBATCH --time=12:00:00\n" +\
                             "ncpus=8\n"+\
                             "export OMP_NUM_THREADS=$ncpus\n" +\
                             "time aprun -b -n 1 -d $ncpus -q prepfold -ncpus $ncpus -n 128 -nsub 128 "+\
                               "-accelcand "+c[1]+" -accelfile "+c[0]+".cand  -o " +\
                               c[0][:-8] + "  " + c[0][:-8] + ".dat\n"
                batch_file.write(batch_line)
            submit_line = 'sbatch batch/fold_' + c[0] + '_' + c[1] +  '.batch'
            submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
        """
                
        print "Number of cands in this file: " + str(len(cand_list))
        fold_num = 1
        job_id_list =[]
        
        last_dm_file = dm_file
         
        #number of folds to do per sbatch job (decrease is using .fits files)        
        cands_per_batch = 100      
        #through some stuffing around sort the fold into 100 folds per job
        for i,c in enumerate(cand_list):
            if (i == 0) or ((i % cands_per_batch) == 0):
                #move to correct file in batch file
                if pulsar == None:
                    dm_i_temp = int(floor( floor(float(c[3]))/2 ))
                    dm_file = dm_i_to_file(dm_i_temp)
                    SUBDIR = work_dir + str(sub_dir) + dm_file
                else:
                    SUBDIR = work_dir + str(sub_dir)[:-1]
                batch_file = open('batch/fold_' + dm_file_orig + '_' + str(fold_num) + '.batch','w')
                #TODO add a check that you're not folding the same thing again
                if pbs:
                    batch_line = "#!/bin/bash -l\n" +\
                                 "#PBS -N fold_" + dm_file_orig + '_' + str(fold_num) + "\n" +\
                                 "#PBS -o " + DIR + "/out/fold_" + dm_file_orig + '_' +\
                                                             str(fold_num) +".out\n" +\
                                 "#PBS -e " + DIR + "/out/fold_" + dm_file_orig + '_' +\
                                                             str(fold_num) +".error\n" +\
                                 "#PBS -l walltime=1:50:00\n"+\
                                 "cd " + SUBDIR + "\n"+\
                                 "module unload pgplot\n"+\
                                 "module load pgplot/x86_64/gnu/5.2-gcc-4.7.1\n"
                    batch_file.write(batch_line)
                else:
                    batch_line = "#!/bin/bash -l\n" +\
                                 "#SBATCH --partition=gpuq\n" +\
                                 "#SBATCH --job-name=fold_" + dm_file_orig + '_' + str(fold_num) + "\n" +\
                                 "#SBATCH --output=" + DIR + "/out/fold_" + dm_file_orig + '_' +\
                                                             str(fold_num) +".out\n" +\
                                 "#SBATCH --time=2:00:00\n" +\
                                 "ncpus=8\n"+\
                                 "export OMP_NUM_THREADS=$ncpus\n" +\
                                 "cd " + SUBDIR + "\n"
                    batch_file.write(batch_line)
                batch_file.write(add_database_function(pbs))
                #fits
                """
                batch_line = "time aprun -b -n 1 -d $ncpus -q prepfold -ncpus $ncpus -n 128 -nsub 128 "+\
                               "-noclip -mask " + str(obsid) + "_rfifind.mask -o " + c[0] + '_' + c[1] +\
                               " -p " + c[3] + " -dm " + c[2] + " -nosearch" +\
                               " /group/mwaops/vcs/" + str(obsid) + "/pointings/" + str(pointing) +\
                               "/" + str(obsid) + "*.fits\n" 
                """
                #batch
                batch_line = 'run "prepfold" "' + ncpuscom + ' -n 128 -nsub 128 '+\
                               "-accelcand "+c[1]+" -accelfile "+c[0]+".cand  -o " +\
                               c[0][:-8] + " " + c[0][:-8] + '.dat" "'+str(bs_id)+'" "'+str(dm_i)+'"\n'
                batch_file.write(batch_line)
                batch_file.close()
            else:
                batch_file = open('batch/fold_' + dm_file_orig + '_' + str(fold_num) + '.batch','a')
                #moves to other directory once it moves past DM%2
                if not dm_file == last_dm_file and pulsar == None:
                    batch_line = "cd " + SUBDIRorig + "\n"
                    batch_file.write(batch_line)
                #fits
                """
                batch_line = "time aprun -b -n 1 -d $ncpus -q prepfold -ncpus $ncpus -n 128 -nsub 128 "+\
                               "-noclip -mask " + str(obsid) + "_rfifind.mask -o " + c[0] + '_' + c[1] +\
                               " -p " + c[3] + " -dm " + c[2] + " -nosearch" +\
                               " /group/mwaops/vcs/" + str(obsid) + "/pointings/" + str(pointing) +\
                               "/" + str(obsid) + "*.fits\n" 
                """
                #batch
                batch_line = 'run "prepfold" "' + ncpuscom + ' -n 128 -nsub 128 '+\
                               "-accelcand "+c[1]+" -accelfile "+c[0]+".cand  -o " +\
                               c[0][:-8] + " " + c[0][:-8] + '.dat" "'+str(bs_id)+'" "'+str(dm_i)+'"\n'
                batch_file.write(batch_line)
                batch_file.close()
            if ((i+1) % cands_per_batch) == 0:
                 
                submit_line = qsub + 'batch/fold_' + dm_file_orig + '_' + str(fold_num) + '.batch'
                submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
                for line in submit_cmd.stdout:
                    print line,
                    if pbs:
                        if '.pbs.hpc.swin.edu.au' in line:
                            job_id_list.append(line[:7])
                    else:
                        if "Submitted" in line:
                            (word1,word2,word3,jobid) = line.split()
                            job_id_list.append(jobid)
                fold_num += 1
                
            last_dm_file = dm_file
        if not ((len(cand_list)+1) % cands_per_batch) == 0: 
            submit_line = qsub + 'batch/fold_' + dm_file_orig + '_' + str(fold_num) + '.batch'
            submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
            for line in submit_cmd.stdout:
                    print line,
                    if pbs:
                        if '.pbs.hpc.swin.edu.au' in line:
                            job_id_list.append(line[:7])
                    else:
                        if "Submitted" in line:
                            (word1,word2,word3,jobid) = line.split()
                            job_id_list.append(jobid)
        
        

        sleep(1)
        job_id_str = ""
        for i in job_id_list:
            job_id_str += ":" + str(i)
        
        if pulsar != None:
            SUBDIRorig = SUBDIR
    
        #TODO check this fires off in the right file
        with open('batch/dependancy_fold_' + dm_file_orig + '.batch','w') as batch_file:
            if pbs:#TODO fix
                batch_line = "#!/bin/bash -l\n" +\
                             "#PBS -N dependancy_fold\n" +\
                             "#PBS -o out/dependancy_fold_" + dm_file_orig + ".out\n" +\
                             "#PBS -e out/dependancy_fold_" + dm_file_orig + ".error\n" +\
                             "#PBS -l walltime=2:05:00\n"+\
                             '#PBS -q gstar\n'+\
                             '#PBS -l nodes=1:ppn=6:gpus=1\n'+\
                             '#PBS -A p125_astro\n' +\
                             "#PBS -W depend=afterany" + job_id_str + "\n" +\
                             "cd " + SUBDIR + "\n"
            else:
                batch_line = "#!/bin/bash -l\n" +\
                             "#SBATCH --job-name=dependancy_fold\n" +\
                             "#SBATCH --output=out/dependancy_fold_" + dm_file_orig + ".out\n" +\
                             "#SBATCH --export=NONE\n" +\
                             "#SBATCH --partition=workq\n" +\
                             "#SBATCH --time=2:05:00\n" +\
                             "#SBATCH --gid=mwaops\n" +\
                             "#SBATCH --account=mwaops\n" +\
                             "#SBATCH --nodes=1\n" +\
                             "#SBATCH --dependency=afterany" + job_id_str + "\n" +\
                             'echo "Searching for a profile with a sigma greater than 3"\n'
            batch_file.write(batch_line)
            if pulsar == None and dm_i > 0:
                i = dm_i * 2 - 1
                batch_line = 'cd '+ SUBDIRpast + '\n' +\
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
                             "done\n"
                batch_file.write(batch_line)
            
            batch_line = 'cd '+ SUBDIRorig + '\n' +\
                         'count=0\n' +\
                         'total=`ls *.ps | wc -l`\n' +\
                         'for i in $(ls *.ps); do\n' +\
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
                         "done\n" +\
                         'blindsearch_database.py -c prepfold -m p -b ' +str(bs_id) + ' -d '+str(dm_i)+'\n'+\
                         "blindsearch_pipeline.py -o " +\
                                     str(obsid) + " -p " + str(pointing) + " -m a -w " + work_dir +\
                                     " -s " + str(sub_dir) + " -r " + str(bs_id)
            batch_file.write(batch_line) 
            if pulsar == None:  
                batch_line = " -d " + str(dm_i + 1) 
                batch_file.write(batch_line)  
            if pbs:
                batch_file.write(" --pbs ")
        submit_line = qsub + 'batch/dependancy_fold_' + dm_file_orig + '.batch'
        submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
        print submit_cmd.communicate()[0],
    else:
        #if there is no cand file assumed there are no cands
        submit_line = "blindsearch_pipeline.py -o " +\
                                     str(obsid) + " -p " + str(pointing) + " -m a -w " + work_dir +\
                                     " -s " + str(sub_dir) + " -d "+str(dm_i + 1)
        
        if pbs:
            submit_line = submit_line + " --pbs "
        submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
        print submit_cmd.communicate()[0], 
    
    return


parser = argparse.ArgumentParser(description="""
Does a blind search for a pulsar in MWA data using the galaxy supercomputer.
""")
parser.add_argument('-o','--observation',type=str,help='The observation ID of the fits file to be searched')
parser.add_argument('-p','--pointing',type=str,help='The pointing of the fits file to be searched')
parser.add_argument('-m','--mode',type=str,help='There are three modes or steps to complete the pipeline. The first mode is to prepdata "p" which dedisperses the the fits files into .dat files. The second mode sort and search "s" which sorts the files it folders and performs a fft. The third mode is accel search "a" runs an accel search on them. The final mode is fold "f" which folds all possible candidates so they can be visaully inspected.')
parser.add_argument('-w','--work_dir',type=str,help='Work directory. Default: /group/mwaops/nswainston/blindsearch/')
parser.add_argument('-s','--sub_dir',type=str,help='Used by the program to keep track of the sub directory its using')
parser.add_argument('-r','--row_num',type=int,help='Database row reference number for keeping track of the scripts.')
parser.add_argument('-d','--dm_file_int',type=int,help='Used by the program to keep track DM file being used to stagger the jobs and not send off over 9000 jobs.')
parser.add_argument('--pulsar',type=str,help="Used to search for a known pulsar by inputing it's Jname. The code then looks within 1 DM and 15%% of the pulsar's period.")
parser.add_argument('--pbs',action="store_true",help="PBS queue mode.")
args=parser.parse_args()


if args.work_dir:
    w_d = args.work_dir
elif args.pbs:
    w_d = '/home/nswainst/blindsearch/'
else:
    w_d = '/group/mwaops/nswainston/blindsearch/'

obs = args.observation
#obsid =  1133329792
point = args.pointing
#19:45:14.00_-31:47:36.00
s_d = args.sub_dir

if args.mode == "r" or args.mode == None:
    print "check"
    rfifind(obs, point, w_d, s_d,args.pbs,args.pulsar)
if args.mode == "p":
    prepdata(obs, point, w_d, s_d,args.row_num,args.pbs,args.pulsar)
elif args.mode == "s":
    sort_fft(obs, point, w_d, s_d,args.row_num,args.pbs,args.pulsar)
elif args.mode == "a":
    accel(obs, point, w_d, s_d,args.dm_file_int,args.row_num,args.pbs,args.pulsar)
elif args.mode == "f":
    fold(obs, point, w_d, s_d,args.dm_file_int,args.row_num,args.pbs,args.pulsar)
    
        
    
#blindsearch_pipeline.py -o 1166459712 -p 06:30:00.0_-28:34:00.0
