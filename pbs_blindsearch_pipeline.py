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
    
    
def add_database_function():
    batch_line ='#PBS -q gstar\n'+\
                '#PBS -l nodes=1:ppn=1,mem=2g\n' +\
                '#PBS -A p125_astro\n' +\
                'cd "$PBS_O_WORKDIR"\n' +\
                'function run\n' +\
                '{\n' +\
                '    # run command and add relevant data to the job database\n' +\
                '    # 1st parameter is command to be run (e.g. wsclean)\n' +\
                '    # 2nd parameter is parameters to that command (e.g. "-j $ncpus")\n' +\
                '    # 3rd parameter is bs_id\n' +\
                '    # 4th parameter is DM file int [optional]\n' +\
                '    if [ "$1" == "rfifind" ] || [ "$1" == "prepsubband" ]; then\n' +\
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
                '    else\n'+\
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
parser.add_argument('--sn_min',type=float,help="Sn_min cut off for accel candidates")
parser.add_argument('--pbs',action="store_true",help="PBS queue mode.")
args=parser.parse_args()


work_dir = '/lustre/projects/p125_astro/blindsearch/'

obsid = args.observation
#obsid =  1133329792
pointing = args.pointing
#19:45:14.00_-31:47:36.00
sub_dir = pointing + "/" + obsid + "/"

pbs = args.pbs
pulsar = args.pulsar

fits_dir = '/lustre/projects/p125_astro/DATA/'
bs_id = args.row_num

if args.mode == "c":
    numout = numout_calc(fits_dir + str(obsid) + "/pointings/" + str(pointing) + "/")
    from math import sqrt,log
    if args.sn_min:
        sn_min = args.sn_min
    else:
        #added a temp 0.8 fudge factor to detect fainter candidates
        #sn_min = ( sqrt(log(numout)) - 0.88 ) / 0.47
        sn_min = 3.0
    print "SN_min: ",sn_min
    
    os.chdir(work_dir + sub_dir)
    cand_list = [] 
    sn_min3_count =0
    submit_cmd = subprocess.Popen('python ~/blindsearch_scripts/ACCEL_sift.py .',shell=True,stdout=subprocess.PIPE)
    print submit_cmd.communicate()[0],
    if os.path.exists('ACCEL_sift_cands.txt'):
        #creat a list of s/n>7 cands
        #[[accel_file,cand,s/n,DM,period(ms)]]
        with open('ACCEL_sift_cands.txt',"rb") as sifted_cands:
            lines = sifted_cands.readlines()
            for l in lines:
                if l.startswith(obsid):
                    cand_line = l.split()
                    if float(cand_line[2]) > sn_min:
                        cand_list.append([cand_line[0].split(':')[0],cand_line[0].split(':')[1],cand_line[2],cand_line[1],cand_line[7]])
                    if float(cand_line[2]) > 3.0:
                        sn_min3_count += 1
        with open('/lustre/projects/p125_astro/blindsearch/ACCEL_sift_sn_list.txt',"a") as sn_list_file:
            sn_list_file.write(str(pointing) + " " +str(sn_min) + " " + str(len(cand_list)) + " " + str(sn_min3_count) +"\n")
    #print cand_list
    output = []
    if len(cand_list) == 0:
        print "No cands over " + str(sn_min)
    print "Candidate number: " + str(len(cand_list))
    output.append(add_database_function())
    for c in cand_list:
        
        batch_line = '~/blindsearch_scripts/run_function.sh "prepfold" " -n 128 -nsub 128 '+\
                               "-noclip -o " + c[0][:-8] + " -p " + str(float(c[4])/1000.) + " -dm " + \
                               c[3] + " /lustre/projects/p125_astro/DATA/" + str(obsid) + "/pointings/" + \
                               str(pointing) + "/" + str(obsid) + '*.fits" "'+str(bs_id)+'" "0"'
        """
        batch_line = 'prepfold -ncpus 12 -n 128 -nsub 128 '+\
                               "-noclip -o " + c[0][:-8] + " -p " + str(float(c[4])/1000.) + " -dm " + \
                               c[3] + " /lustre/projects/p125_astro/DATA/" + str(obsid) + "/pointings/" + \
                               str(pointing) + "/" + str(obsid) + '*.fits'
        
        batch_line = '~/blindsearch_scripts/run_function.sh "prepfold" " -ncpus 12 -n 128 -nsub 128 '+\
                           "-accelcand "+c[1]+" -accelfile "+c[0]+".cand  -o " +\
                           c[0][:-8] + " " + c[0][:-8] + '.dat" "'+str(bs_id)+'" "0"'
        
        
        batch_line = 'run "prepfold" " -n 128 -nsub 128 '+\
                           "-accelcand "+c[1]+" -accelfile "+c[0]+".cand  -o " +\
                           c[0][:-8] + " " + c[0][:-8] + '.dat" "'+str(bs_id)+'" "0"'
        """
        print batch_line
        output.append(batch_line)
    
    for o in output:
        submit_cmd = subprocess.Popen(o,shell=True,stdout=subprocess.PIPE)
        for line in submit_cmd.stdout:
            print line,
    #print output[-1]
    
                        
else:
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


    
    qsub = 'qsub '
    ncpuscom = ''
    #ncpuscom = ' -ncpus'
    """
    #Get the centre freq channel and then run DDplan.py to work out the most effective DMs
    print "Obtaining metadata from http://mwa-metadata01.pawsey.org.au/metadata/ for OBS ID: " + str(obsid)
    beam_meta_data = getmeta(service='obs', params={'obs_id':obsid})
    channels = beam_meta_data[u'rfstreams'][u"0"][u'frequencies']
    minfreq = float(min(channels))
    maxfreq = float(max(channels))
    centrefreq = 1.28 * (minfreq + (maxfreq-minfreq)/2) #in MHz


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
    dm_line = [['1.00','100.00','0.20','1','495','1']]
            
    #Calculates -numout for prepsubbands
    numout = numout_calc(fits_dir + str(obsid) + "/pointings/" + str(pointing) + "/")
    from math import sqrt,log

     
    #super script
    with open('batch/super_script.batch','w') as batch_file:
        batch_line = "#!/bin/bash\n" +\
                     "#PBS -N " + pointing + "_blinsearch_pipeline\n" +\
                     "#PBS -o batch/superscript.out\n" +\
                     "#PBS -e batch/superscript.error\n" +\
                     "#PBS -l walltime=1:00:00\n" 
        batch_file.write(batch_line)
        batch_file.write(add_database_function())
        batch_line = 'run "prepsubband" "'+ncpuscom + '-lodm ' + str(dm_line[0][0]) + " -dmstep " +\
                        str(dm_line[0][2]) + " -numdms " +str(dm_line[0][4]) + " -numout " +\
                        str(numout) + " -o " + str(obsid) +\
                         " -mask /lustre/projects/p125_astro/DATA/1166459712/ics/1166459712_rfifind.mask "+\
                         fits_dir + str(obsid) + \
                        "/pointings/" + str(pointing) + "/" + str(obsid) + '*.fits" '+str(bs_id) + "\n" +\
                     'blindsearch_database.py -c prepsubband -m p -b ' +str(bs_id) + '\n'+\
                     'for i in $(seq '+ str(dm_line[0][0]) + ' ' + str(dm_line[0][2]) +' '+\
                                str(dm_line[0][1])+'); do \n'+\
                     'run "realfft" "' + str(obsid) + '_DM${i}.dat" "'+str(bs_id)+'" "0"\n'+\
                     'run "accelsearch" "'  + ncpuscom + ' -flo 1 -fhi 500 '+\
                                '-numharm 32 -zmax 0 ' + str(obsid) + '_DM${i}.fft" "' +str(bs_id) +\
                                '" "0";done\n'+\
                     'blindsearch_database.py -c realfft -m m -b ' +str(bs_id) + '\n' +\
                     'blindsearch_database.py -c realfft -m p -b ' +str(bs_id) + '\n'+\
                     'blindsearch_database.py -c accelsearch -m m -b ' +str(bs_id) + '\n' +\
                     'blindsearch_database.py -c accelsearch -m p -b ' +str(bs_id) +' -d 0\n'        
        batch_file.write(batch_line)
    submit_line = qsub + 'batch/super_script.batch'
    submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
    for line in submit_cmd.stdout:
        print line,
        if '.pbs.hpc.swin.edu.au' in line:
            job_id = line[:7]


#export OMP_NUM_THREADS=8

    #sort by DM
    from operator import itemgetter
    #print cand_list
    with open('batch/fold_script.batch','w') as batch_file:
        batch_line = "#!/bin/bash -l\n" +\
                     "#PBS -N fold_"+ str(pointing) + "\n" +\
                     "#PBS -o batch/fold_.out\n" +\
                     "#PBS -e batch/fold_.error\n" +\
                     "#PBS -l walltime=6:00:00\n"+\
                     "#PBS -W depend=afterany:" + str(job_id) + "\n" 
        batch_file.write(batch_line)
        batch_file.write(add_database_function())
        batch_line = "module unload pgplot\n"+\
                     "module load pgplot/x86_64/gnu/5.2-gcc-4.7.1\n" +\
                     'pbs_blindsearch_pipeline.py -o ' + str(obsid) + ' -p ' + str(pointing) + \
                            ' -m c -r ' + str(bs_id) + '\n' +\
                     'count=0\n' +\
                     'total=`ls *.ps | wc -l`\n' +\
                     'for i in $(ls *.ps); do\n' +\
                     'if (( $count % 100 == 0 )); then\n' +\
                     'echo "$count / $total searched"\n' +\
                     'fi\n' +\
                     'chi=`sed "13q;d" ${i%.ps}.bestprof`\n' +\
                     "if [ ${chi:20:3} -ge 3 ]; then\n" +\
                     'ps_to_png.sh ${i}\n' +\
                     'mv "${i%.ps}".png ../../over_2_png/'+str(pointing)+'_"${i%.ps}".png\n' +\
                     'echo "${i%.ps}.png is over 3"\n' +\
                     "fi\n" +\
                     "count=$(($count+1))\n" +\
                     "done\n" +\
                     'blindsearch_database.py -c prepfold -m m -b ' +str(bs_id) + ' -d 0\n'+\
                     'blindsearch_database.py -c prepfold -m p -b ' +str(bs_id) + ' -d 0\n'+\
                     'single_pulse_search.py *dat'+\
                     'ps_to_png.sh '+str(obsid)+'_singlepulse.ps\n' +\
                     'mv '+str(obsid)+'_singlepulse.png ../../over_2_png/'+str(pointing)+'_'+str(obsid)+'_singlepulse.png\n' 
        batch_file.write(batch_line)
    submit_line = 'qsub batch/fold_script.batch'
    submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
    print submit_cmd.communicate()[0],
        
