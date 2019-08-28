#! /usr/bin/env python3

import subprocess
import os
import argparse
import textwrap
import pipes
import glob
from time import sleep
import datetime
import socket
import numpy as np

#vcstools imports
import search_database
import mwa_metadb_utils as meta
from mwa_metadb_utils import get_channels
import process_vcs as pvcs
from job_submit import submit_slurm
import config

import logging
logger = logging.getLogger(__name__)

if 'SEARCH_WORK_DIR' in os.environ:
    DEFAULT_WORK_DIR = os.environ['SEARCH_WORK_DIR']
else:
    DEFAULT_WORK_DIR = '/fred/oz125/nswainst/pulsar_search'
if 'SEARCH_DB' in os.environ:
    DB_FILE_LOC = os.environ['SEARCH_DB']
else:
    print("#WARNING: No environment variable DB_FILE_LOC. Search tracking won't work")

DB_FILE_LOC = os.environ['SEARCH_DB']
DEFAULT_WORK_DIR = os.environ['SEARCH_WORK_DIR']


class search_options_class:
    # Initializer / Instance Attributes
    def __init__(self, obsid,
                 pointing=None, cal_id=None, begin=None, end=None,
                 channels=None, incoh=False, args=None,
                 work_dir=None, sub_dir=None, fits_dir_base=None,
                 dm_min=0, dm_max=250,
                 relaunch_script=None, search=False, bsd_row_num=None, cold_storage_check=False,
                 table='Prepdata', attempt=0,
                 nice=None, search_ver='master',
                 DI_dir=None, pointing_dir=None, n_omp_threads=None):
        hostname = socket.gethostname()
        comp_config = config.load_config_file()

        #obs/beamforming options
        self.obsid    = obsid
        self._pointing = pointing
        self.cal_id   = cal_id
        self.begin    = begin
        self.end      = end
        self.channels = get_channels(obsid, channels=channels)
        self.incoh    = incoh
        self.args     = args

        #directories
        if DI_dir is None:
            self.DI_dir = "{0}{1}/cal/{2}/rts/".format(comp_config['base_product_dir'],
                                                       args.observation, args.cal_obs)
        else:
            self.DI_dir    = DI_dir
        if work_dir is None:
            self.work_dir  = DEFAULT_WORK_DIR
        else:
            self.work_dir  = work_dir
        self.sub_dir       = sub_dir
        if fits_dir_base is None:
            self.fits_dir_base = '{0}{1}/'.format(comp_config['base_product_dir'], obsid)
        else:
            self.fits_dir_base = fits_dir_base
        if pointing_dir is None:
            if incoh:
                self._pointing_dir = '{0}incoh/'.format(self.fits_dir_base)
            elif self.pointing is not None:
                self._pointing_dir = '{0}pointings/{1}/'.format(self.fits_dir_base,
                                                               self.pointing)
            else:
                self._pointing_dir = None
        else:
            self._pointing_dir = pointing_dir

        #search parameters
        self.dm_min           = dm_min
        self.dm_max           = dm_max

        #search database
        self._relaunch_script    = relaunch_script
        self.search             = search
        self._bsd_row_num        = bsd_row_num
        self.cold_storage_check = cold_storage_check
        self._table              = table
        self.attempt            = attempt

        #job options
        if nice is None:
             if hostname.startswith('john') or hostname.startswith('farnarkle'):
                 self.nice = 0
             else:
                 self.nice = 100
        else:
            self.nice = nice
        self.search_ver = search_ver
        if n_omp_threads is None:
            if hostname.startswith('john') or hostname.startswith('farnarkle'):
                #on ozstar use single core jobs
                self._n_omp_threads = 1
            else:
                #on galaxy use all 8 cores because can't split them up
                self._n_omp_threads = 8
        else:
            self._n_omp_threads = n_omp_threads

    #Set up variables that need to be editted
    def getTable(self):
        return self._table
    def setTable(self, value):
        self._table = value
    table = property(getTable, setTable)

    def getBRN(self):
        return self._bsd_row_num
    def setBRN(self, value):
        self._bsd_row_num = value
    bsd_row_num = property(getBRN, setBRN)

    def getPdir(self):
        return self._pointing_dir
    def setPdir(self, value):
        self._pointing_dir = value
    pointing_dir = property(getPdir, setPdir)

    def getPoint(self):
        return self._pointing
    def setPoint(self, value):
        self._pointing = value
    pointing = property(getPoint, setPoint)

    def getNOT(self):
        return self._n_omp_threads
    def setNOT(self, value):
        self._n_omp_threads = value
    n_omp_threads = property(getNOT, setNOT)

    def getRLS(self):
        return self._relaunch_script
    def setRLS(self, value):
        self._relaunch_script = value
    relaunch_script = property(getRLS, setRLS)



def send_cmd_shell(cmd):
    output = subprocess.Popen(cmd, stdin=subprocess.PIPE,
                              stdout=subprocess.PIPE, shell=True,
                              stderr=subprocess.STDOUT).communicate()[0].decode()
    return output

def send_cmd(cmd):
    output = subprocess.Popen(cmd.split(' '), stdin=subprocess.PIPE,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT).communicate()[0].decode()
    return output


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


def your_slurm_queue_check(max_queue = 200, queue=None, grep=None):
    """
    Checks if you have over 100 jobs on the queue, if so waits until your queue clears
    """
    submit_line = 'squeue -u $USER'
    if queue is not None:
        submit_line += ' --partition={0}'.format(queue)
    else:
        queue = ""
    if grep is not None:
        submit_line += ' | grep "{0}"'.format(grep)
    submit_line += ' | wc -l'
    submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
    q_num = ""
    for line in submit_cmd.stdout:
        q_num += line.decode()
    #remove header line
    q_num = int(q_num) - 1
    while (q_num > max_queue ):
        print("{}/{} jobs on the {} queue. Waiting 2 mins for queue to clear".format(q_num,
              max_queue, queue))
        sleep(120)
        submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
        q_num = ""
        for line in submit_cmd.stdout:
            q_num += line.decode()
        q_num = int(q_num) - 1
    print("{}/{} jobs on the {} queue. Continuing".format(q_num, max_queue, queue))
    return


def add_database_function():
    batch_line ='function run\n' +\
                '{\n' +\
                '    # run command and add relevant data to the job database\n' +\
                '    # 1st parameter is command\n' +\
                '    # 2nd parameter is argument\n' +\
                '    # 3rd parameter is bsd_row_num\n' +\
                '    # 4th parameter is script_row_num\n' +\
                '    # 5th parameter is attempt_num\n' +\
                '    echo search_database.py -m s -c $1 -b $3 -r $4 -a $5 \n' +\
                '    search_database.py -m s -c $1 -b $3 -r $4 -a $5 \n' +\
                '    echo $1 $2\n' +\
                '    srun --export=ALL -n 1 -c $ncpus $1 $2\n' +\
                '    errcode=$?\n' +\
                '    search_database.py -m e -c $1 -b $3 -r $4 -a $5 --errorcode $errcode\n' +\
                '    if [ "$errcode" != "0" ]; then\n' +\
                '        exit $errcode\n' +\
                '    fi\n' +\
                '}\n' +\
                '\n'
    return batch_line


def add_temp_database_function(job_num, attempt_num, n_omp_threads=1):
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
                '    echo $1 $2\n' +\
                '    $1 $2\n' +\
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
    dirlist = glob.glob("{0}/{1}_*.fits".format(fits_dir, obsid))
    if not dirlist:
        print("ERROR. No files in {0}/{1}_*.fits".format(fits_dir, obsid))
        exit()
    numout = 0
    for d in dirlist:
        submit_line = 'readfile ' + d
        submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
        for line in submit_cmd.stdout:
            if b'Time per file (sec) =' in line:
                subint = int(line.split(b'=')[-1])
        numout += int(subint * 1e4)

    if numout%2:
        numout += 1
    return numout

def get_pulsar_dm_p(pulsar):
    #Gets the ra and dec from the output of PSRCAT
    cmd = 'psrcat -c dm {}'.format(pulsar)
    output = send_cmd(cmd)
    lines = output.split('\n')
    for l in lines[4:-1]:
        columns = l.split()

        if len(columns) > 1:
            dm = columns[1]
    cmd = 'psrcat -c p0 {}'.format(pulsar)
    output = send_cmd(cmd)
    lines = output.split('\n')
    for l in lines[4:-1]:
        columns = l.split()
        if len(columns) > 1:
            p = columns[1]
    return [dm, p]

def process_vcs_wrapper(search_opts, pointings,
                        pulsar_list_list=None,
                        vdif=False, summed=False,
                        code_comment=None, pointing_id=None,
                        channels=None):
    """
    Does some basic checks and formating before
    if args.pulsar_file:
        code_comment += using beamforming from process_vcs.py
    """
    comp_config = config.load_config_file()
    #check queue
    your_slurm_queue_check(queue=comp_config['gpuq_partition'], max_queue=70)
    your_slurm_queue_check(max_queue=500)

    #check for search_opts.incoh file which is used to predict if you have used rfifind
    search_opts.incoh_check = False
    bf_formats = " -p"
    if not os.path.exists('{0}{1}/incoh'.format(comp_config['base_product_dir'], search_opts.obsid)):
        bf_formats += " -i"
        os.mkdir('{0}{1}/incoh'.format(comp_config['base_product_dir'], search_opts.obsid))
        search_opts.incoh_check = True
    if vdif:
        bf_formats += " -u"
    if summed:
        bf_formats += " -s"

    #set up and launch beamfroming
    comp_config = config.load_config_file()
    data_dir = '{0}{1}'.format(comp_config['base_data_dir'], search_opts.obsid)
    metafits_dir = "{0}/{1}_metafits_ppds.fits".format(data_dir, search_opts.obsid)
    rts_flag_file = '{0}/cal/{1}/rts/flagged_tiles.txt'.\
                     format(search_opts.fits_dir_base, search_opts.cal_id)
    if not os.path.exists(rts_flag_file):
        print('RTS flag file is not in the default location of: '
              '{} please move it there. Exiting'.format(rts_flag_file))
        exit()
    pvcs.ensure_metafits(data_dir, search_opts.obsid, metafits_dir)
    channels = get_channels(search_opts.obsid, channels=channels)
    job_id_list_list = pvcs.coherent_beam(search_opts.obsid, search_opts.begin, search_opts.end,
                      data_dir, search_opts.fits_dir_base,
                      "{0}/batch".format(search_opts.fits_dir_base),
                      metafits_dir, 128, pointings, search_opts.args,
                      rts_flag_file=rts_flag_file, bf_formats=bf_formats,
                      DI_dir=search_opts.DI_dir,
                      calibration_type="rts", nice=search_opts.nice,
                      vcstools_version="multi-pixel_beamform",
                      channels=channels)

    code_comment_in = code_comment
    dep_job_id_list = []
    for job_id_list in job_id_list_list:
        for pn, pointing in enumerate(pointings):
            search_opts.setPoint(pointing)
            if pulsar_list_list is None:
                pulsar_list = None
            else:
                pulsar_list = pulsar_list_list[pn]
            if code_comment_in is not None:
                code_comment = "{0} pn {1}".format(code_comment_in, pointing_id[pn])
            search_opts.setPdir('{0}/pointings/{1}'.format(search_opts.fits_dir_base,
                                                                  search_opts.pointing))

            search_opts.setBRN(search_database.database_search_start(search_opts.obsid,
                                          search_opts.pointing, "{0}".format(code_comment)))
            dep_job_id_list.append(dependant_splice_batch(search_opts, job_id_list=job_id_list,
                                                          pulsar_list=pulsar_list))
    return dep_job_id_list

def multibeam_binfind(search_opts, pointing_dir_list, job_id_list, pulsar, loglvl="INFO"):
    """
    Takes many pointings and launches data_processing_pipeline which folds on all of the pointings and finds the best one. This will by default continue running the processing pipeline
    """
    pointing_str = " ".join(pointing_dir_list)
    logger.info("pointing string: {0}".format(pointing_str))
    commands = []
    commands.append("echo 'Folding on multiple pointings'")
    commands.append("data_process_pipeline.py -m m -d {0} -o {1} -O {2} -p {3} -L {4} "
                    "--mwa_search {5}".format(pointing_str, search_opts.obsid,
                    search_opts.cal_id, pulsar, loglvl, search_opts.search_ver))

    name="multibeam_fold_{0}_{1}".format(pulsar, search_opts.obsid)
    batch_dir = "{0}/batch/".format(search_opts.fits_dir_base)
    submit_slurm(name, commands,\
                batch_dir=batch_dir,\
                slurm_kwargs={"time": "00:05:00"},\
                module_list=['mwa_search/{0}'.format(search_opts.search_ver),\
                              'presto/no-python'],\
                submit=True, vcstools_version='multi-pixel_beamform',\
                depend=job_id_list)


def dependant_splice_batch(search_opts, job_id_list=None, pulsar_list=None):
    """
    Launches a script that splices the beamformed files and, where approriate,
    launches the search pipeline or folds on known pulsars.
    """
    comp_config = config.load_config_file()

    #create a split wrapper dependancy
    splice_wrapper_batch = 'splice_wrapper_{0}_{1}'.format(search_opts.obsid, search_opts.pointing)
    if os.path.exists("{0}/batch/".format(search_opts.fits_dir_base)):
        batch_dir = "{0}/batch/".format(search_opts.fits_dir_base)
    else:
        batch_dir = search_opts.fits_dir_base

    commands = []
    if search_opts.bsd_row_num is not None:
        #record beamforming processing time
        commands.append('search_database.py -m b -b {0} -f {1}mb_{2}'.\
                        format(search_opts.bsd_row_num, batch_dir, search_opts.pointing))

    #grabbing chan ids for splicing
    search_opts.channels = get_channels(search_opts.obsid, channels=search_opts.channels)

    #add splice command
    splice_command = 'splice_wrapper.py -o {0} -d -c {1}'.format(search_opts.obsid,
                     ' '.join(map(str, search_opts.channels)))
    commands.append('{0} -w {1}'.format(splice_command, search_opts.pointing_dir))
    if search_opts.incoh:
        commands.append('{0} -i -w {1}{2}/incoh/'.format(splice_command,
                         comp_config['base_product_dir'], search_opts.obsid))

    """
    if pulsar_list is not None:
        #check_known_pulsars.py uses this to check if it was detected and if so upload it
        commands.append('cd {0}'.format(search_opts.pointing_dir))
        for pulsar in pulsar_list:
            #Assign number of bins based on the period
            period = get_pulsar_dm_p(pulsar)[1]
            nbin = 2048*float(period)**0.75
            nbin = int(32 * round(nbin/32))
            if nbin<32:
                nbin=32

            #load presto module here because it uses python 2
            commands.append('echo "Folding on known pulsar {0}"'.format(pulsar))
            commands.append('psrcat -e {0} > {0}.eph'.format(pulsar))
            commands.append("sed -i '/UNITS           TCB/d' {0}.eph".format(pulsar))
            commands.append("prepfold -o {0} -noxwin -runavg -noclip -timing {1}.eph -nsub 256\
                            {2}/1*fits -n {3}".format(search_opts.obsid, pulsar,
                                                      search_opts.pointing_dir, nbin))
            commands.append('errorcode=$?')
            commands.append('pulsar={}'.format(pulsar[1:]))
            pulsar_bash_string = '${pulsar}'
            #Somre old ephems don't have the correct ra and dec formating and
            #causes an error with -timing but not -psr
            commands.append('if [ "$errorcode" != "0" ]; then')
            commands.append('   echo "Folding using the -psr option"')
            commands.append('   prepfold -o {0} -noxwin -runavg -noclip -psr {1} -nsub 256\
                            {2}/1*fits -n {3}'.format(search_opts.obsid, pulsar,
                                                      search_opts.pointing_dir, nbin))
            commands.append('   pulsar={}'.format(pulsar))
            commands.append('fi')
            commands.append('rm {0}.eph'.format(pulsar))

            commands.append('echo "Checking profile that it has a chi of over 4"')
            commands.append('chi=`sed "13q;d" {0}_PSR_{1}.pfd.bestprof`'.format(search_opts.obsid,
                                                                                pulsar_bash_string))
            commands.append('chi=${chi#*=}')
            commands.append('echo "Chi value of ${chi%.*}"')
            commands.append('if [ ${chi%.*} -ge 4 ]; then')
            commands.append('   echo "Strong detection. Uploading to MWA Pulsar Database"')
            commands.append('   submit_to_database.py -o {0} --cal_id {1} -p {2} '
                                '--bestprof {0}_PSR_{3}.pfd.bestprof --ppps '
                                '{0}_PSR_{3}.pfd.ps'.format(search_opts.obsid,
                                search_opts.cal_id, pulsar, pulsar_bash_string))
            commands.append('   echo "Searching for pulsar using the pipeline to test the '
                                      'pipelines effectivness"')
            commands.append('   mwa_search_pipeline.py -o {0} -b {1} -e {2} --search '
                                '--pulsar {3} -O {4} --code_comment "Known pulsar auto '
                                'test {3}"'.format(search_opts.obsid, search_opts.begin,
                                            search_opts.end, pulsar, search_opts.cal_id))
            commands.append("fi")
    """

    #add search_opts.relaunch script
    if search_opts.relaunch_script is not None:
        relaunch_script = "{0} -p {1} -s {1}/{2}".format(search_opts.relaunch_script,
                                            search_opts.pointing, search_opts.obsid)
        if search_opts.bsd_row_num is not None:
            relaunch_script += ' -r {0}'.format(search_opts.bsd_row_num)
        if search_opts.incoh:
            #run rfi job
            relaunch_script += ' -m r --incoh'
        else:
            relaunch_script += ' -m b'
        commands.append(relaunch_script)

    hostname = socket.gethostname()
    if hostname.startswith('john') or hostname.startswith('farnarkle'):
        mem = 2048
        temp_mem = 48
    else:
        mem = 1024
        temp_mem = None

    job_id = submit_slurm(splice_wrapper_batch, commands,
                 batch_dir=batch_dir, nice=search_opts.nice,
                 slurm_kwargs={"time": "8:00:00"},
                 module_list=['mwa_search/{}'.format(search_opts.search_ver),
                              'presto/no-python'],
                 submit=True, depend=job_id_list, depend_type='afterany',
                 mem=mem, temp_mem=temp_mem,
                 vcstools_version="multi-pixel_beamform")
    return job_id


def beamform(search_opts, pointing_list, code_comment=None,
             pulsar_list_list=None, vdif=False, summed=False,
             relaunch=False):
    search_opts.channels = get_channels(search_opts.obsid, channels=search_opts.channels)
    bsd_row_num_input = search_opts.bsd_row_num
    pointing_dir_input = search_opts.pointing_dir

    #create a dictionary for recording pointdirs and pulsars to make a depdant fold job
    pulsar_fold_dict = {}
    if pulsar_list_list is None:
        pulsar_fold_dict[None] = []
    else:
        for pulsar_list in pulsar_list_list:
            pulsar_fold_dict[" ".join(pulsar_list)] = []

    #work out maximum number of pointings
    hostname = socket.gethostname()
    if ( hostname.startswith('john') or hostname.startswith('farnarkle') ) and summed:
        #temp_mem = int(5. * (float(stop) - float(start) + 1.) * \
        #           float(len(pointing_list)) / 1000.) + 1
        temp_mem_max = 300. #GB
        # Use the maximum number of pointings that the SSD memory can handle
        max_pointing = int(( temp_mem_max - 1 ) * 1000. / \
                           (5. * (float(search_opts.end) - float(search_opts.begin) + 1.)))
        if max_pointing > 29:
            # More than 30 won't fit on the GPU mem
            max_pointing = 29
    else:
        max_pointing = 15

    print("Maximum of pointings per beamforming job: {}".format(max_pointing))

    pointings_to_beamform = []
    pulsar_list_list_to_beamform = []
    search_opts.relaunch_script_in = search_opts.relaunch_script
    for n, line in enumerate(pointing_list):
        if line.startswith("#"):
            continue
        print("Checking pointing {0} out of {1}".format(n+1, len(pointing_list)))
        if search_opts.incoh:
            search_opts.setPoint("incoh")
        elif ':' not in line:
            #search_opts.pointing = search_opts.fits_dir_base.split("/")[-1]
            search_opts.setPoint(line)
        else:
            ra, dec = line.split(" ")
            if dec.endswith("\n"):
                dec = dec[:-1]
            search_opts.setPoint(ra + "_" + dec)

        #pulsar check parsing
        if pulsar_list_list is None:
            pulsar_list = None
            pulsar = None
        else:
            pulsar_list = pulsar_list_list[n]
            pulsar = pulsar_list_list[n][0]
        logger.debug("pulsar_list: {}".format(pulsar_list))

        #set up search_opts.relaunch scripts
        if pulsar is None:
            search_opts.sub_dir = '{0}/{1}'.format(search_opts.pointing, search_opts.obsid)
        else:
            search_opts.sub_dir = '{0}/{1}'.format(pulsar, search_opts.obsid)
        search_opts.relaunch_script = "{0} -p {1} -s {2}".format(search_opts.relaunch_script_in,
                                                      search_opts.pointing, search_opts.sub_dir)

        #fits dir parsing
        comp_config = config.load_config_file()
        if pointing_dir_input is None:
            if search_opts.incoh:
                search_opts.setPdir('{0}incoh/'.format(search_opts.fits_dir_base))
            else:
                search_opts.setPdir('{0}pointings/{1}/'.format(search_opts.fits_dir_base,
                                                                   search_opts.pointing))
        else:
            search_opts.setPdir(pointing_dir_input)

        if search_opts.cold_storage_check:
            #Check if search_opts.pointing in cold storage
            try :
                exists_remote_check = exists_remote("hpc-hsm.pawsey.org.au",
                        "/project/mwaops/nswainston/yogesh_low_DM_candiate/{0}_pointing.tar.gz".\
                        format(search_opts.pointing))
                if exists_remote_check and len(pointing_list) > 1:
                    print("The pointing is in cold storage so assumed it is analysised "
                          "so not reprocessing")
                    continue
            except Exception('SSH failed'):
                print("Connection to cold storage failed. Will only check for local files")

        #Go through some file checks (true is something is missing)
        path_check = False
        missing_file_check = False
        unspliced_check = False
        searched_check = False
        missing_chan_list = []
        searched_check = search_database.database_search_done_check(search_opts.obsid,
                                                                    search_opts.pointing)
        if os.path.exists(search_opts.pointing_dir):
            #first check is there's already spliced files
            #does check if they have the same start time
            expected_file_num = int( (search_opts.end-search_opts.begin)/200 ) + 2
            for fnc in range(1,expected_file_num):
                logger.debug(search_opts.pointing_dir+search_opts.obsid+"_*"+str(fnc)+".fits")
                if not glob.glob(search_opts.pointing_dir+search_opts.obsid+"_*"+str(fnc)+".fits"):
                    missing_file_check = True
            #if search_opts.fits_dir_base is not None:
            #    #assumes that maybe they didn't use the whole obs and that's ok
            #    if len(glob.glob(search_opts.pointing_dir+search_opts.obsid+"_*.fits")) > 0:
            #        missing_file_check = False

            if missing_file_check:
                #check if we have any unspliced files
                #there are some so going to resubmit jobs
                search_opts.channels = get_channels(search_opts.obsid, channels=search_opts.channels)
                job_id_list =[]
                unspliced_check = False
                for ch in search_opts.channels:
                    for ne in range(1,expected_file_num):
                        logger.debug("{0}*_{1}_{2}_ch*{3}_00*{4}.fits".format(
                                     search_opts.pointing_dir, search_opts.pointing,
                                     search_opts.obsid, ch, ne))
                        if not glob.glob("{0}*_{1}_{2}_ch*{3}_00*{4}.fits".format(
                                         search_opts.pointing_dir, search_opts.obsid,
                                         search_opts.pointing, ch, ne)):
                            unspliced_check = True
                            if ch not in missing_chan_list:
                                logger.debug("missing chan {0}".format(ch))
                                missing_chan_list.append(ch)
            else:
                # Replacing the above with a check of the number of samples
                expected_nsamples = (search_opts.end-search_opts.begin) * 10000
                nsamples = numout_calc(search_opts.pointing_dir, search_opts.obsid)
                if expected_nsamples > nsamples:
                    print("Not enough fits files so deleteing fits files and beamforming")
                    path_check = True
                    for rm_file in glob.glob("{}/{}*fits".format(search_opts.pointing_dir,
                                                                 search_opts.obsid)):
                        os.remove(rm_file)

        else:
            path_check = True


        if ((missing_file_check and not unspliced_check and search_opts.search) or \
            (search_opts.search and ((not searched_check or relaunch) \
                or len(pointing_list) == 1) ) )\
            and bsd_row_num_input is None:
                search_opts.setBRN(search_database.database_search_start(search_opts.obsid,
                                   search_opts.pointing, "{0} pn {1}".format(code_comment,n)))
                search_opts.setRLS("{0} -r {1}".format(search_opts.relaunch_script,
                                                       search_opts.bsd_row_num))
        else:
            search_opts.setBRN(bsd_row_num_input)


        #work out what needs to be done
        if search_opts.search and searched_check and not relaunch:
            print("Already searched so not searching again")
        elif missing_file_check and not unspliced_check and not path_check:
            #splice files
            print("Splicing the files in {0}".format(search_opts.pointing))
            dep_job_id = dependant_splice_batch(search_opts, pulsar_list=pulsar_list)
            if pulsar_list is None:
                pulsar_fold_dict[pulsar_list].append([search_opts.pointing_dir, dep_job_id])
            else:
                pulsar_fold_dict[" ".join(pulsar_list)].append([search_opts.pointing_dir,
                                                                dep_job_id])
        elif path_check or len(missing_chan_list) == 24:
            logger.debug(missing_file_check, unspliced_check, path_check, len(missing_chan_list))
            # do beamforming
            print("No pointing directory or files for {0}, will beamform shortly".\
                    format(search_opts.pointing))
            pointings_to_beamform.append(search_opts.pointing)
            #pulsars that will be checked after beamforming
            if pulsar_list is None:
                pulsar_list_list_to_beamform = None
            else:
                pulsar_list_list_to_beamform.append(pulsar_list)
        elif unspliced_check:
            #resubmit any search_opts.channels that are incomplete
            print("Some channels missing, beamforming on {0} for {1}".format(missing_chan_list,
                  search_opts.pointing))
            if len(pointing_list) > 1:
                your_slurm_queue_check(queue=comp_config['gpuq_partition'], max_queue=70)
            
            temp_pointing_id = pointing_list.index(search_opts.pointing.replace("_", " ")) + 1
            dep_job_id = process_vcs_wrapper(search_opts, [search_opts.pointing],
                                pulsar_list_list=[pulsar_list],
                                vdif=vdif, summed=summed,
                                code_comment=code_comment,
                                pointing_id=temp_pointing_id,
                                channels=missing_chan_list)
            
            logger.debug(pulsar_list, search_opts.pointing, dep_job_id)
            pointing_dir_temp = '{0}/pointings/{1}'.format(search_opts.fits_dir_base,
                                                           search_opts.pointing)
            if pulsar_list is None:
                pulsar_fold_dict[pulsar_list].append([pointing_dir_temp, dep_job_id])
            else:
                pulsar_fold_dict[" ".join(pulsar_list)].append([pointing_dir_temp,
                                                                dep_job_id])


        else:
            #All files there so the check has succeded and going to start the pipeline
            if search_opts.search and ((not searched_check or relaunch)\
                   or len(pointing_list) == 1):
                print("Fits files available, begining pipeline for {0}".format(search_opts.pointing))
                if len(pointing_list) > 1:
                    your_slurm_queue_check(max_queue=500)
                prepdata(search_opts)
            else:
                print("Fits files available, not beamforming or searching")
                if pulsar_list is None:
                    pulsar_fold_dict[pulsar_list].append([search_opts.pointing_dir, None])
                else:
                    pulsar_fold_dict[" ".join(pulsar_list)].append([search_opts.pointing_dir, None])

            #remove any extra unspliced files
            for fr in glob.glob(search_opts.pointing_dir+"*_"+search_opts.obsid+"_*.fits"):
                os.remove(fr)

            #Sending off the splice wrapper just for the folding
            #splice_wrapper should fail
            #if pulsar_list is not None:
            #    dependant_splice_batch(search_opts, pulsar_list=pulsar_list)

        if ( n + 1 == len(pointing_list) and len(pointings_to_beamform) != 0 )\
            or len(pointings_to_beamform) == max_pointing:
            #Send of beamforming job at the end or the loop or when you have 15 or 30 pointings
            # list of search_opts.pointing ids, always 1 for single search_opts.pointings and
            # will be relivant to the line number of search_opts.pointing files
            pointing_id = []
            for point in pointings_to_beamform:
                pointing_id.append(pointing_list.index(point.replace("_", " ")) + 1)

            print("Sending of {0} pointings for beamforming".format(len(pointings_to_beamform)))
            dep_job_id_list = process_vcs_wrapper(search_opts, pointings_to_beamform,
                                pulsar_list_list=pulsar_list_list_to_beamform,
                                vdif=vdif, summed=summed,
                                code_comment=code_comment, pointing_id=pointing_id)
            logger.debug("pulsar_list_list_to_beamform: {}".format(pulsar_list_list_to_beamform))
            if pulsar_list_list is not None:
                for pulsar_list, pointing, dep_job_id in zip(pulsar_list_list_to_beamform,
                                                             pointings_to_beamform,
                                                             dep_job_id_list):
                    logger.debug(pulsar_list, pointing, dep_job_id)
                    pointing_dir_temp = '{0}/pointings/{1}'.format(search_opts.fits_dir_base,
                                                                   pointing)
                    if pulsar_list is None:
                        pulsar_fold_dict[pulsar_list].append([pointing_dir_temp, dep_job_id])
                    else:
                        pulsar_fold_dict[" ".join(pulsar_list)].append([pointing_dir_temp,
                                                                        dep_job_id])
            pointings_to_beamform = []
            pulsar_list_list_to_beamform = []
    #send off pulsar fold jobs
    if pulsar_list_list is not None:
        for pulsar_list in pulsar_fold_dict:
            pulsar_list = pulsar_list.split()
            logger.debug("pulsar_list: {}".format(pulsar_list))
            pointing_dir_list = []
            dep_job_id_list = []
            dict_list = pulsar_fold_dict[" ".join(pulsar_list)]
            logger.debug("dict_list: {}".format(dict_list))
            for pointing_dir, dep_job_id in dict_list:
                pointing_dir_list.append(pointing_dir)
                if dep_job_id is not None:
                    dep_job_id_list.append(dep_job_id)

            if len(dep_job_id_list) == 0:
                dep_job_id_list = None
            for pulsar in pulsar_list:
                print("Sending off data processing for {0}".format(pulsar))
                logger.debug("{0} Pointing list: {1}".format(pulsar, pointing_dir_list))
                #send of a fold job for each pulsar
                multibeam_binfind(search_opts, pointing_dir_list, dep_job_id_list, pulsar)
    return


#-------------------------------------------------------------------------------------------------------------
def rfifind(search_opts):
    comp_config = config.load_config_file()

    #Calculates -numout for prepsubbands
    numout = numout_calc(search_opts.pointing_dir, search_opts.obsid)

    #Set up some directories and move to it
    if not os.path.exists("{0}{1}".format(search_opts.work_dir, search_opts.pointing)):
        os.mkdir("{0}{1}".format(search_opts.work_dir, search_opts.pointing))
    if not os.path.exists("{0}{1}/{2}".format(search_opts.work_dir, search_opts.pointing,
                          search_opts.obsid)):
        os.mkdir("{0}{1}/{2}".format(search_opts.work_dir, search_opts.pointing,
                 search_opts.obsid))
    if not os.path.exists("{0}{1}/{2}/batch".format(search_opts.work_dir,
                          search_opts.pointing, search_opts.obsid)):
        os.mkdir("{0}{1}/{2}/batch".format(search_opts.work_dir,
                 search_opts.pointing, search_opts.obsid))

    hostname = socket.gethostname()
    if hostname.startswith('john') or hostname.startswith('farnarkle'):
        #fft needs more memory, only need to change on ozstar
        mem=2048
    else:
        mem=1024


    rfi_batch = str(search_opts.bsd_row_num) + '_rfi_{0}'.format(search_opts.obsid)
    commands = []
    commands.append("ncpus={0}".format(search_opts.n_omp_threads))
    commands.append("export OMP_NUM_THREADS={0}".format(search_opts.n_omp_threads))
    if not os.path.exists("{0}/rfi_masks/{1}_rfifind.mask".format(search_opts.work_dir,
                                                                  search_opts.obsid)):
        #if there is not already a rfi mask make one
        commands.append('cd {0}{1}/{2}'.format(search_opts.work_dir, search_opts.pointing,
                                               search_opts.obsid))
        commands.append('rfifind -ncpus $ncpus -noclip -time 12.0 -o {0} -zapchan 0:19,108:127,'
                '128:147,236:255,256:275,364:383,384:403,492:511,512:531,620:639,640:659,748:767,'
                '768:787,876:895,896:915,1004:1023,1024:1043,1132:1151,1152:1171,1260:1279,1280:'
                '1299,1388:1407,1408:1427,1516:1535,1536:1555,1644:1663,1664:1683,1772:1791,1792:'
                '1811,1900:1919,1920:1939,2028:2047,2048:2067,2156:2175,2176:2195,2284:2303,2304:'
                '2323,2412:2431,2432:2451,2540:2559,2560:2579,2668:2687,2688:2707,2796:2815,2816:'
                '2835,2924:2943,2944:2963,3052:3071 {1}*.fits'.format(search_opts.obsid,
                                                                      search_opts.pointing_dir))
        commands.append('mv {0}_rfifind.mask {1}/rfi_masks/'.format(search_opts.obsid,
                                                                    search_opts.work_dir))
        commands.append('search_database.py -c rfifind -m p -b ' +str(search_opts.bsd_row_num))
        commands.append("prepdata -ncpus $ncpus -dm 0 " "-numout {0} -o {1}_DM0.00 "
                        "{2}*incoh*.fits".format(numout, search_opts.obsid,
                                                 search_opts.pointing_dir))
        commands.append('mv {0}_rfifind.* {1}/rfi_masks/'.format(search_opts.obsid,
                                                                 search_opts.work_dir))
        commands.append('mv {0}_DM0.00.dat {1}/rfi_masks/'.format(search_opts.obsid,
                                                                  search_opts.work_dir))
        commands.append('mv {0}_DM0.00.inf {1}/rfi_masks/'.format(search_opts.obsid,
                                                                  search_opts.work_dir))
    #commands.append("{0} -m p -r {1} -s {2}/{3}".format(search_opts.relaunch_script,
                     #search_opts.bsd_row_num, search_opts.pointing, search_opts.obsid))

    submit_slurm(rfi_batch, commands,
                 batch_dir="{0}{1}/{2}/batch".format(search_opts.work_dir,
                                search_opts.pointing, search_opts.obsid),
                 slurm_kwargs={"time": "2:00:00"}, mem=mem,
                 nice=search_opts.nice,
                 submit=True, module_list=[comp_config['presto_module']])

    return


#-------------------------------------------------------------------------------------------------------------
def prepdata(search_opts):
    comp_config = config.load_config_file()

    #Set up some directories and move to it
    if not os.path.exists("{0}/rfi_masks".format(search_opts.work_dir)):
        os.mkdir("{0}/rfi_masks".format(search_opts.work_dir))
    #make subdir
    logger.debug("{0}{1}".format(search_opts.work_dir, search_opts.sub_dir.split("/")[0]))
    if not os.path.exists("{0}{1}".format(search_opts.work_dir, search_opts.sub_dir.split("/")[0])):
        os.mkdir("{0}{1}".format(search_opts.work_dir, search_opts.sub_dir.split("/")[0]))
    if not os.path.exists("{0}{1}".format(search_opts.work_dir, search_opts.sub_dir)):
        os.mkdir("{0}{1}".format(search_opts.work_dir, search_opts.sub_dir))

    os.chdir(search_opts.work_dir + search_opts.sub_dir)

    if not os.path.exists("{0}{1}/batch".format(search_opts.work_dir, search_opts.sub_dir)):
        os.mkdir("{0}{1}/batch".format(search_opts.work_dir, search_opts.sub_dir))


    #Get the centre freq channel and then run DDplan.py to work out the most effective DMs
    search_opts.channels = get_channels(search_opts.obsid, channels=search_opts.channels)
    minfreq = float(min(search_opts.channels))
    maxfreq = float(max(search_opts.channels))
    centrefreq = 1.28 * (minfreq + (maxfreq-minfreq)/2) #in MHz

    #Old DDplan method
    """
    output = subprocess.Popen(['DDplan.py','-l',str(search_opts.dm_min),
                               '-d',str(search_opts.dm_max),
                               '-f',str(centrefreq),
                               '-b','30.7200067160534',
                               '-t','0.0001',
                               '-n','3072',
                               '-o','dm_temp'],stdout=subprocess.PIPE).communicate()
    subprocess.check_call("\n", shell=True)
    os.remove('dm_temp.eps')
    dm_list = []
    print(output[0])
    lines = output[0].split('\n')
    for l in lines[13:-4]:
        columns = l.split()
        dm_list.append(columns)

    #dm_list = [['1.500','3.500','0.01','4','200','1']]
    """

    #new lfDDplan method #TODO make this more robust

    from lfDDplan import dd_plan
    dm_list = dd_plan(centrefreq, 30.72, 3072, 0.1, search_opts.dm_min, search_opts.dm_max)
    #dm_list = [[low_dm, high_dm, DM_step, number_of_steps, time_res]]


    #Calculates -numout for prepsubbands
    numout = numout_calc(search_opts.pointing_dir, search_opts.obsid)

    #Calculates the expected procesing time (conservative) in seconds
    #the number of DMs doesn't appear to greatly affect the processing time so not included in the calculation
    expe_proc_time = 3.*10**11 / numout
    dm_list.sort()
    print(dm_list)
    dm_list_list = [] # a range of DMs for each prepsubband command

    hostname = socket.gethostname()
    if hostname.startswith('john') or hostname.startswith('farnarkle'):
        #If on ozstar use their SSD to improve I/O
        SSD_file_dir = '$JOBFS/'
    else:
        SSD_file_dir = ''
    #Create a list of all the commands needed
    if os.path.exists('{0}rfi_masks/{1}_rfifind.mask'.format(search_opts.work_dir,
                                                             search_opts.obsid)):
        mask_command = ' -mask {0}rfi_masks/{1}_rfifind.mask'.format(search_opts.work_dir,
                                                                     search_opts.obsid)
    else:
        mask_command = ''
    #dms_per_job = 1024
    dms_per_job = 256
    nsub = 24
    commands_list = []
    for dm_line in dm_list:
        dm_start = dm_line[0]
        dm_end = dm_line[1] #float(dm_line[2]) * float(dm_line[4]) + float(dm_start)
        while ( (dm_end - float(dm_start)) / float(dm_line[2])) > float(dms_per_job) :
            #dedisperse for only 1024 steps
            commands_list.append('-ncpus $ncpus -lodm {0} {1} -nsub {2} -dmstep {3} '
                '-numdms {4} -numout {5} -zerodm -o {6}{7} {8}{7}_*.fits'.format(dm_start,
                mask_command, nsub, dm_line[2], dms_per_job+1, numout, SSD_file_dir,
                search_opts.obsid, search_opts.pointing_dir))
            dm_list_list.append(np.around(np.arange(float(dm_start),
                                          float(dm_start) + float(dms_per_job) * float(dm_line[2]),
                                          float(dm_line[2])), decimals=2))
            dm_start = str(float(dm_start) + (float(dms_per_job) * float(dm_line[2])))
        steps = int((dm_end - float(dm_start)) / float(dm_line[2]))
        #last loop to get the <1024 steps
        commands_list.append('-ncpus $ncpus -lodm {0} {1} -nsub {2} -dmstep {3} '
                '-numdms {4} -numout {5} -zerodm -o {6}{7} {8}{7}_*.fits'.format(dm_start,
                mask_command, nsub, dm_line[2], steps+1, numout, SSD_file_dir,
                search_opts.obsid, search_opts.pointing_dir))
        dm_list_list.append(np.around(np.arange(float(dm_start),
                                      float(dm_start) + float(steps) * float(dm_line[2]),
                                      float(dm_line[2])), decimals=2))

    if "-r" not in search_opts.relaunch_script:
        search_opts.setRLS("{0} -r {1}".format(search_opts.relaunch_script, search_opts.bsd_row_num))


    search_opts.setTable('Prepdata')
    hostname = socket.gethostname()
    if hostname.startswith('john') or hostname.startswith('farnarkle'):
        #It is more efficient on ozstar to use their SSDs for the intermediate files
        #such as .dat and fft files so the PRESTO commands must be run in series in a
        #single job
        sort_fft(search_opts, dm_list_list=dm_list_list, prepsub_commands=commands_list)
    else:
        #Puts all the expected jobs on the databse
        #search_database_script_id_list
        search_database.database_script_list(search_opts.bsd_row_num, 'prepsubband', commands_list,
                             search_opts.n_omp_threads, expe_proc_time)


        error_check(search_opts, bash_job=True, total_job_time=7200)
    return

#-------------------------------------------------------------------------------------------------------------
def sort_fft(search_opts, dm_list_list=None, prepsub_commands=None):

    #Makes 90 files to make this all a bit more managable and sorts the files.
    os.chdir(search_opts.work_dir + "/" + search_opts.sub_dir)
    if not os.path.exists("over_3_png"):
        os.mkdir("over_3_png")
    if not os.path.exists("{}/cand_files".format(search_opts.work_dir)):
        os.mkdir("{}/cand_files".format(search_opts.work_dir))

    DIR=search_opts.work_dir + search_opts.sub_dir
    os.chdir(DIR)

    hostname = socket.gethostname()
    if hostname.startswith('john') or hostname.startswith('farnarkle'):
        #If on ozstar use their SSD to improve I/O
        SSD_file_dir = '$JOBFS/'
    else:
        SSD_file_dir = ''
    #Set up jobs on database
    expe_proc_time = 180. #giving it a generous 60 seconds as these jobs don't take long at all
    commands_list = []
    if dm_list_list is None:
        dat_files = glob.glob(search_opts.work_dir + search_opts.sub_dir + "/*dat")
        dat_files.sort()
        for dat in dat_files:
            commands_list.append(dat.split("/")[-1])
    else:
        #no files yet in single job methog
        dm_list = []
        for dml in dm_list_list:
            dm_list.extend(dml)
        for dm in dm_list:
            commands_list.append("{0}{1}_DM{2:.2f}.dat".format(SSD_file_dir, search_opts.obsid, dm))

    #Send off jobs
    search_opts.setTable('FFT')
    hostname = socket.gethostname()
    if hostname.startswith('john') or hostname.startswith('farnarkle'):
        accel(search_opts, dm_list_list=dm_list_list,
              prepsub_commands=prepsub_commands, fft_commands=commands_list)
    else:
        search_database.database_script_list(search_opts.bsd_row_num, 'realfft', commands_list,
                                             search_opts.n_omp_threads, expe_proc_time)

        error_check(search_opts, bash_job=True, total_job_time=3600)

    print("Sent off fft jobs")
    return


#-------------------------------------------------------------------------------------------------------------
def accel(search_opts, dm_list_list=None, prepsub_commands=None, fft_commands=None):

    os.chdir(search_opts.work_dir + search_opts.sub_dir)

    nharm = 16. # number of harmonics to search
    min_period = 0.001 # min period to search for in sec (ANTF min = 0.0013)
    max_period = 30. # max period to search for in sec  (ANTF max = 23.5)
    #convert to freq
    min_freq = 1. / max_period
    max_freq = 1. / min_period
    #adjust the freq to include the harmonics
    min_f_harm = min_freq
    max_f_harm = max_freq * nharm

    #For initial search we will save processing by not doing an acceleration search
    max_search_accel = 0

    hostname = socket.gethostname()
    if hostname.startswith('john') or hostname.startswith('farnarkle'):
        #If on ozstar use their SSD to improve I/O
        SSD_file_dir = '$JOBFS/'
    else:
        SSD_file_dir = ''

    commands_list = []
    if dm_list_list is None:
        dir_files = glob.glob("*fft")
    else:
        dm_list  = []
        dir_files = []
        for dml in dm_list_list:
            dm_list.extend(dml)
        for dm in dm_list:
            dir_files.append("{0}{1}_DM{2:.2f}.fft".format(SSD_file_dir, search_opts.obsid, dm))

    for fft_file in dir_files:
        #calc processing time
        expe_proc_time = 300.
        commands_list.append('-ncpus $ncpus -zmax {0:d} -flo {1} -fhi {2} -numharm {3:d} {4}'.\
                             format(int(max_search_accel), min_f_harm, max_f_harm,
                                    int(nharm), fft_file))


    #Send off jobs
    if dm_list_list is None:
        search_opts.setTable('Accel')
        error_check(search_opts, bash_job=True)
        search_database.database_script_list(search_opts.bsd_row_num, 'accelsearch',
                                             commands_list, search_opts.n_omp_threads,
                                             expe_proc_time)


    else:
        presto_single_job(search_opts, dm_list_list,
                          prepsub_commands=prepsub_commands,
                          fft_commands=fft_commands,
                          accel_commands=commands_list)

    print("Sent off accel jobs")
    return

#-------------------------------------------------------------------------------------------------------------
def fold(search_opts):

    DIR=search_opts.work_dir + str(search_opts.sub_dir)
    os.chdir(DIR)
    if not os.path.exists("presto_profiles"):
        os.mkdir("presto_profiles")

    #run accel_sift.py to find candidates with DM structure
    file_loc = '{0}cand_files/cands_{1}.txt'.format(search_opts.work_dir,
                                        search_opts.sub_dir.replace("/","_"))

    #calcs sn_min for candidates
    numout = numout_calc(search_opts.pointing_dir, search_opts.obsid)
    from math import sqrt,log
    #sn_min = ( sqrt(log(numout)) - 0.88 ) / 0.47
    sn_min = 5.

    cand_list = []
    import shutil
    if shutil.which('ACCEL_sift.py') is None:
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
        os.chdir(search_opts.work_dir)
        #This is now done in it's own batch script
        #print("Removing Candidates without DM structure using ACCEL_sift.py")
        #print(submit_line)
        #submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
        #print(submit_cmd.communicate()[0],)
        if os.path.exists(file_loc):
            #creat a list of s/n>7 cands
            #[[accel_file,cand,s/n,DM,period(ms)]]
            with open(file_loc,"r") as sifted_cands:
                lines = sifted_cands.readlines()
                for l in lines:
                    if l.startswith(search_opts.obsid):
                        cand_line = l.split()
                        if float(cand_line[2]) > sn_min:
                            cand_list.append([cand_line[0].split(':')[0],cand_line[0].split(':')[1],
                                              cand_line[2],cand_line[1],cand_line[7]])
        else:
            print("Can't find ACCEL cand file: {}. ACCEL_sift.py likely failed".format(file_loc))
            print("Exiting here")
            exit()
        os.chdir(DIR + "/presto_profiles")

    #Uses the mask if it's available
    if os.path.exists('{0}rfi_masks/{1}_rfifind.mask'.format(search_opts.work_dir, search_opts.obsid)):
        mask_command = ' -mask {0}rfi_masks/{1}_rfifind.mask'.format(search_opts.work_dir,
                                                                     search_opts.obsid)
    else:
        mask_command = ''

    #cand_list = [accel_file_name, cand_num, SN, DM, period(ms)]
    print("Ssearch_opts.ending off jobs with fft sn greater than {}".format(sn_min))
    expe_proc_time = 5400.
    commands_list = []
    if len(cand_list) > 0:
        #sort by DM
        from operator import itemgetter
        cand_list.sort(key=itemgetter(3))
        print("Number of cands in this file: " + str(len(cand_list)))

        for c in cand_list:
            accel_file_name, cand_num, _, cand_DM, period = c
            #through some stuffing around sort the fold into 100 folds per job
            #the fold option using .dat files which is quicker but inaccurate
            #fold_command = 'run "prepfold" "-ncpus $ncpus -n 128 -nsub 128 '+\
            #           "-accelcand "+c[1]+" -accelfile "+c[0]+".cand  -o " +\
            #           c[0][:-8] + " " + c[0][:-8] + '.dat" "'+str(search_opts.bsd_row_num)+'" "'+str(dm_i)+'"'

            #the fold options that uses .fits files which is slower but more accurate
            if float(cand_DM) > 1.:
                commands_list.append('-n 128 -noxwin -noclip -o {0}_{1}_{2} -p {3} -dm {4} '
                    '-nosearch {5} {6}{7}_*.fits'.format(accel_file_name, cand_num,
                    search_opts.pointing, float(period)/1000., cand_DM, mask_command,
                    search_opts.pointing_dir, search_opts.obsid))
            else:
                print("Skipping cand with DM {0} and period {1} ms".format(cand_DM, period))
    search_database.database_script_list(search_opts.bsd_row_num, 'prepfold', commands_list,
                                                  search_opts.n_omp_threads, expe_proc_time)
    if len(cand_list) > 0:
        #Send off jobs
        search_opts.setTable('Fold')
        error_check(search_opts)
    else:
        wrap_up(search_opts)



    return



def wrap_up(search_opts):

    #put all significant candidates into the over directory


    #update final database logs
    comp_config = config.load_config_file()
    wrap_batch = str(search_opts.bsd_row_num) + "_wrap_up"
    commands = []
    commands.append(add_database_function())
    commands.append('cd {0}{1}/presto_profiles'.format(search_opts.work_dir, search_opts.sub_dir))
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
    commands.append('       if [ -f "${i%.ps}".png ]; then')
    commands.append('           mv "${i%.ps}".png ../over_3_png/"${i%.ps}".png')
    commands.append('           echo "${i%.ps}.png is over 3"')
    commands.append('       else')
    commands.append('           mv ${i} ../over_3_png/${i}')
    commands.append('           echo "${i} is over 3"')
    commands.append('       fi')
    commands.append('       mv "${i%.ps}" ../over_3_png/"${i%.ps}"')
    commands.append('   fi')
    commands.append('   count=$(($count+1))')
    commands.append('done')
    commands.append('if ls ../over_3_png/*.png 1> /dev/null 2>&1; then')
    commands.append('   over_sn=`ls ../over_3_png/*.png | wc -l`')
    commands.append('else')
    commands.append('   over_sn=`ls ../over_3_png/*.ps | wc -l`')
    commands.append('fi')
    commands.append('if [ $over_sn -eq 0 ]; then')
    commands.append('   echo "No candidates so deleting pointing, dat and fft files"')
    commands.append('   rm {0}{1}/*dat'.format(search_opts.work_dir, search_opts.sub_dir))
    commands.append('   rm {0}{1}/*fft'.format(search_opts.work_dir, search_opts.sub_dir))
    commands.append('   rm {0}{1}/*ACCEL_0*'.format(search_opts.work_dir, search_opts.sub_dir))
    commands.append('   rm -rf {0}{1}/{2}'.format(comp_config['base_product_dir'],
                                                  search_opts.obsid, search_opts.pointing))
    commands.append('fi')
    commands.append('search_database.py -m w -b {0} --cand_val "$total $over_sn 0"'.\
                    format(search_opts.bsd_row_num))
    submit_slurm(wrap_batch, commands,
                 batch_dir="{0}{1}/batch".format(search_opts.work_dir, search_opts.sub_dir),
                 slurm_kwargs={"time": "2:50:00"}, nice=search_opts.nice,
                 submit=True, module_list=['mwa_search/{}'.format(search_opts.search_ver)])
    return

def presto_single_job(search_opts, dm_list_list, prepsub_commands=None,
                      fft_commands=None, accel_commands=None):
    """
    A simpler version of error_check() that sends off prepsubband, fft and accelsearch
    commands one after the other to take advantage of Ozstars SSDs
    """
    job_id_list = []
    # Get prepsubband commands
    if prepsub_commands is None:
        prepsub_commands = search_database.database_script_check('Prepdata',
                                                search_opts.bsd_row_num, 1)
    # Get fft commands
    if fft_commands is None:
        fft_commands = search_database.database_script_check('FFT', search_opts.bsd_row_num, 1)

    # Get accel commands
    if accel_commands is None:
        accel_commands = search_database.database_script_check('Accel', search_opts.bsd_row_num, 1)

    temp_mem = 100 #GB

    dat_start = 0 #id of the first file to use dat/fft
    for dmi, command_data in enumerate(prepsub_commands):
        processing_time = 0.0
        check_batch = "{0}_presto_a{1}_{2}".format(search_opts.bsd_row_num,
                                        search_opts.attempt+1, dmi)
        commands = []
        commands.append(add_database_function())
        commands.append('cp -r $TEMPO2 $JOBFS/tempo2')
        commands.append('cp -r $TEMPO2_CLOCK_DIR $JOBFS/tempo2_clock_dir')
        commands.append('export TEMPO2=$JOBFS/tempo2')
        commands.append('export TEMPO2_CLOCK_DIR=$JOBFS/tempo2_clock_dir')
        commands.append('cd {0}{1}/'.format(search_opts.work_dir,
                                               search_opts.sub_dir))
        commands.append('')

        #add prepsubband command
        if isinstance(command_data, list):
            commands.append('run "{0}" "{1}" "{2}" "{3}" "{4}"'.format(command_data[0],
                                 command_data[1], search_opts.bsd_row_num, dmi,
                                 search_opts.attempt+1))
            processing_time += float(command_data[2])
        else:
            commands.append("prepsubband {}".format(command_data))

        #work out the DMs of the commands to add
        dat_num = len(dm_list_list[dmi])
        dat_range = range(dat_start, dat_start + dat_num)

        #make the fft bash file
        with open('{0}{1}/{2}_fft_a{3}_{4}.bash'.format(search_opts.work_dir,
                  search_opts.sub_dir, search_opts.bsd_row_num,
                  search_opts.attempt+1, dmi), "w") as srun_bash:
            srun_bash.write(add_temp_database_function(dmi,
                            search_opts.attempt + 1,
                            n_omp_threads=1))
            for ffti in dat_range:
                command_fft = fft_commands[ffti]
                if isinstance(command_fft, list):
                    srun_bash.write('run "{0}" "{1}" "{2}" "{3}" "{4}"\n'.format(
                                     command_fft[0], command_fft[1], search_opts.bsd_row_num,
                                     ffti, search_opts.attempt+1))
                    processing_time += command_fft[2]
                else:
                    srun_bash.write("realfft {}\n".format(command_fft))
        commands.append("srun --export=ALL -n 1 -c 1 bash {0}_fft_a{1}_{2}.bash".\
                        format(search_opts.bsd_row_num, search_opts.attempt+1, dmi))

        #make the accel bash file
        with open('{0}{1}/{2}_accel_a{3}_{4}.bash'.format(search_opts.work_dir,
                  search_opts.sub_dir, search_opts.bsd_row_num,
                  search_opts.attempt+1, dmi), "w") as srun_bash:
            srun_bash.write(add_temp_database_function(dmi,
                            search_opts.attempt + 1,
                            n_omp_threads=search_opts.n_omp_threads))
            for ai in dat_range:
                command_accel = accel_commands[ai]
                if isinstance(command_accel, list):
                    srun_bash.write('run "{0}" "{1}" "{2}" "{3}" "{4}"\n'.format(
                                     command_accel[0], command_accel[1], search_opts.bsd_row_num,
                                     ai, search_opts.attempt+1))
                    processing_time += command_accel[2]
                else:
                    srun_bash.write("accelsearch {}\n".format(command_accel))


        commands.append("srun --export=ALL -n 1 -c {0} bash {1}_accel_a{2}_{3}.bash".\
                        format(search_opts.n_omp_threads, search_opts.bsd_row_num,
                               search_opts.attempt+1, dmi))

        # load python modules required to do the single pulse python script
        commands.append('ml numpy/1.14.1-python-2.7.14')
        commands.append('ml scipy/1.0.0-python-2.7.14')
        commands.append('single_pulse_search.py $JOBFS/*.dat')
        commands.append('cp $JOBFS/{0}_singlepulse.ps {1}{2}'.format(search_opts.obsid,
                        search_opts.work_dir, search_opts.sub_dir))

        #Remove accel files off ssd
        commands.append('cp $JOBFS/*ACCEL* {0}{1}'.format(search_opts.work_dir,
                                                          search_opts.sub_dir))
        commands.append('cp $JOBFS/*inf {0}{1}'.format(search_opts.work_dir,
                                                          search_opts.sub_dir))

        #load python modules and run database scripts
        #TODO Temporarily removed to avoid database timeout issues
        '''
        commands.append('module load mwa_search/{0}'.format(search_opts.search_ver))
        commands.append('search_database.py -c realfft -m m --file_location '
                        'realfft_temp_database_file_{0}_{1}.csv\n'.format(
                             search_opts.attempt + 1, dmi))
        commands.append('search_database.py -c accelsearch -m m --file_location '
                        'accelsearch_temp_database_file_{0}_{1}.csv\n'.format(
                             search_opts.attempt + 1, dmi))
        '''
        if processing_time == 0:
            processing_time = 10800
        if processing_time > 10800:
            processing_time = 10800 #max at 3 hours because that should be plenty of time
        total_job_time_str = datetime.timedelta(seconds=int(processing_time))
        job_id = submit_slurm(check_batch, commands,
                 batch_dir="{0}{1}/batch".format(search_opts.work_dir,search_opts.sub_dir),
                 slurm_kwargs={"time": total_job_time_str},
                 module_list=['presto/min_path', 'psrcat/1.49',
                              'tempo2/278e129', 'tempo2_clock/2019-03-21'],
                 submit=True, cpu_threads=search_opts.n_omp_threads,
                 mem=2048, temp_mem=temp_mem,
                 vcstools_version=None)
        job_id_list.append(job_id)

        dat_start += dat_num


    #Dependancy job
    print("Waiting 1 sec to make sure to the dependancy script works")
    sleep(1)

    check_depend_batch = '{0}_dep_presto_a{1}'.format(search_opts.bsd_row_num,
                                                      search_opts.attempt +1)
    commands = []
    #TODO temporarily moving right to fold to avoid database issues
    #commands.append("{0} --attempt {1} -m c --table Accel".format(search_opts.relaunch_script,
    #                                                              search_opts.attempt + 1))
    #TODO ONLY NEED FOR database removal
    module_list = []
    # on ozstar so use their modules
    module_list.append('mwa_search/{}'.format(search_opts.search_ver))
    module_list.append('module use /apps/users/pulsar/skylake/modulefiles\n'
                       'module load presto/d6265c2')
    module_list.append('matplotlib/2.2.2-python-2.7.14')
    #find ACCEL_sift path
    import shutil
    accel_sift = shutil.which("ACCEL_sift.py")
    commands = []
    commands.append(add_database_function())
    commands.append('cd {0}'.format(search_opts.work_dir))
    commands.append('srun --export=ALL -n 1 -c 1 {0} {1}/'.format(accel_sift,
                                                              search_opts.sub_dir))
    commands.append('module use {}'.format(comp_config['module_dir']))
    commands.append('module load mwa_search/{}'.format(search_opts.search_ver))
    #TODO end temp sec

    commands.append("{0} -m f -r {1}".format(search_opts.relaunch_script, search_opts.bsd_row_num))

    submit_slurm(check_depend_batch, commands,
                 batch_dir="{0}{1}/batch".format(search_opts.work_dir, search_opts.sub_dir),
                 slurm_kwargs={"time": "20:00"}, nice=search_opts.nice,
                 submit=True, depend=job_id_list, depend_type="afterany",
                 #module_list=['mwa_search/{}'.format(search_opts.search_ver)],
                 module_list=module_list,
                 cpu_threads=search_opts.n_omp_threads)

    return

def error_check(search_opts, bash_job=False,
                mem=1024, total_job_time=18000.):
    """
    Checkes the database for any jobs that didn't complete (or didn't even start)
    and reruns any errors before continuing to the next step
    """
    comp_config = config.load_config_file()
    hostname = socket.gethostname()

    sub_sub_dir = ''
    bash_job = False
    temp_mem = None
    if search_opts.table == 'Prepdata':
        total_job_time = 24000
        next_mode = 't'
        cur_mode = 'p'
        if hostname.startswith('john') or hostname.startswith('farnarkle'):
            mem = 2048
    elif search_opts.table == 'FFT':
        search_opts.setNOT(1) #fft is not parrelised
        if hostname.startswith('john') or hostname.startswith('farnarkle'):
            #fft needs more memory, only need to change on ozstar
            mem = 2048
        total_job_time = 6000
        bash_job = True
        next_mode = 'a'
        cur_mode = 't'
    elif search_opts.table == 'Accel':
        bash_job = True
        next_mode = 'f'
        cur_mode = 'a'
    elif search_opts.table == 'Fold':
        sub_sub_dir = 'presto_profiles'
        next_mode = 'w'
        cur_mode = 'f'

    logger.debug("Table: {0}".format(search_opts.table))
    total_job_time_str = datetime.timedelta(seconds=total_job_time)
    if search_opts.attempt == 0:
        print("Launching {} scripts for the first time.".format(search_opts.table))
        error_data = search_database.database_script_check(search_opts.table,
                                                           search_opts.bsd_row_num, 1)
    else:
        print("Checking for any {} scripts that failed.".format(search_opts.table))
        error_data = search_database.database_script_check(search_opts.table,
                                search_opts.bsd_row_num, search_opts.attempt)
    if len(error_data) == 0:
        print("No incomplete jobs, moving on to next part of the pipeline")
        if cur_mode == 'a':
            # changing to python 2 to use ACCELsift
            print('Sending off job to run ACCELsift')

            module_list = []
            if hostname.startswith('john') or hostname.startswith('farnarkle'):
                # on ozstar so use their modules
                module_list.append('mwa_search/{}'.format(search_opts.search_ver))
                module_list.append('module use /apps/users/pulsar/skylake/modulefiles\n'
                                   'module load presto/d6265c2')
                module_list.append('matplotlib/2.2.2-python-2.7.14')
            else:
                # use galaxy modules
                module_list.append('module unload vcstools')
                module_list.append('matplotlib')
                module_list.append('python/2.7.14')
                module_list.append('numpy')
                module_list.append('presto')

            #find ACCEL_sift path
            import shutil
            accel_sift = shutil.which("ACCEL_sift.py")
            commands = []
            commands.append(add_database_function())
            commands.append('cd {0}'.format(search_opts.work_dir))
            commands.append('srun --export=ALL -n 1 -c 1 {0} {1}/'.format(accel_sift,
                                                                      search_opts.sub_dir))
            commands.append('module use {}'.format(comp_config['module_dir']))
            if not (hostname.startswith('john') or hostname.startswith('farnarkle')):
                # If on galaxy it sometimes needs the python version explictedly stated
                commands.append('module unload matplotlib')
                commands.append('module unload numpy')
                commands.append('module unload presto')
                commands.append('module unload python/2.7.14')
                commands.append('PYTHON_VERSION=3.6.3')
                commands.append('MAALI_PYTHON_LIB_VERSION=3.6')
                commands.append('module load python/3.6.3')
                commands.append('module load scipy')
            commands.append('module load mwa_search/{}'.format(search_opts.search_ver))
            commands.append("{0} -m {1}".format(search_opts.relaunch_script, next_mode))
            submit_slurm("{0}_ACCEL_sift".format(search_opts.bsd_row_num), commands,
                         batch_dir="{0}{1}/batch".format(search_opts.work_dir, search_opts.sub_dir),
                         slurm_kwargs={"time": '59:00'}, nice=search_opts.nice,
                         module_list=module_list,
                         cpu_threads=1, mem=mem, submit=True,
                         vcstools_version="Error_on_purpose")
        else:
            print("{0} -m {1}".format(search_opts.relaunch_script, next_mode))
            print(send_cmd_shell("{0} -m {1}".format(search_opts.relaunch_script, next_mode)))
    elif search_opts.attempt > 10:
        print("Still failing after 10 attempts. Exiting Here.")
    else:
        print("Submiting {} commands".format(len(error_data)))
        presto_command = error_data[0][0]
        if search_opts.attempt != 0 :
            command_list = []
            for er in error_data:
               command_list.append(er[1])
            search_database.database_script_list(search_opts.bsd_row_num, presto_command,
                      command_list, search_opts.n_omp_threads, error_data[0][2],
                      attempt=(search_opts.attempt+1))

        #submit jobs
        check_job_num = 1
        job_id_list = []

        processing_time = 0.0
        check_batch = "{0}_{1}_a{2}_{3}".format(search_opts.bsd_row_num, presto_command,
                            search_opts.attempt+1, check_job_num)
        commands = []
        bash_commands = []

        if not bash_job:
            commands.append(add_database_function())
        commands.append('cd {0}{1}/{2}'.format(search_opts.work_dir,
                                         search_opts.sub_dir, sub_sub_dir))

        if hostname.startswith('john') or hostname.startswith('farnarkle'):
            #If on Ozstar move tempo 2 files onto the SSD to make their frequent
            #reads and write faster
            temp_mem = 1 #GB
            commands.append('cp -r $TEMPO2 $JOBFS/tempo2')
            commands.append('cp -r $TEMPO2_CLOCK_DIR $JOBFS/tempo2_clock_dir')
            commands.append('export TEMPO2=$JOBFS/tempo2')
            commands.append('export TEMPO2_CLOCK_DIR=$JOBFS/tempo2_clock_dir')

        for ei, er in enumerate(error_data):
            bash_commands.append('run "{0}" "{1}" "{2}" "{3}" "{4}"'.format(presto_command,
                                er[1], search_opts.bsd_row_num, ei, (search_opts.attempt+1)))
            processing_time += float(er[2])
            if search_opts.table == 'Prepdata' and \
               (hostname.startswith('john') or hostname.startswith('farnarkle')):
                #prepdata is very I/O heavy so on ozstar write to the SSDs then move it off
                bash_commands.append('echo "Moving data off SSD onto {0}{1}/{2}"'.format(
                                     search_opts.work_dir, search_opts.sub_dir, sub_sub_dir))
                bash_commands.append('time cp $JOBFS/* {0}{1}/{2}'.format(search_opts.work_dir,
                                                           search_opts.sub_dir, sub_sub_dir))
                bash_commands.append('echo "Finished moving data off SSD"')
                temp_mem = 51 #GB

            # if next step will go over the expected time limit then create a new script
            if (processing_time + float(er[2])) > total_job_time:
                if bash_job:
                    #if it's a bash job make a file for it to run
                    with open('{0}{1}/{2}.bash'.format(search_opts.work_dir,
                              search_opts.sub_dir, check_batch), "w") as srun_bash:
                        srun_bash.write(add_temp_database_function(check_job_num,
                                        search_opts.attempt + 1,
                                        n_omp_threads=search_opts.n_omp_threads))
                        for sc in bash_commands:
                            srun_bash.write("{}\n".format(sc))

                    commands.append("srun --export=ALL -n 1 -c {0} bash {1}.bash".\
                                    format(search_opts.n_omp_threads, check_batch))
                    commands.append('search_database.py -c {0} -m m --file_location '
                                    '{0}_temp_database_file_{1}_{2}.csv\n'.format(
                                        presto_command, search_opts.attempt + 1, check_job_num))
                else:
                    commands = commands + bash_commands

                job_id = submit_slurm(check_batch, commands,
                         batch_dir="{0}{1}/batch".format(search_opts.work_dir,search_opts.sub_dir),
                         slurm_kwargs={"time": total_job_time_str},
                         nice=search_opts.nice,
                         module_list=['mwa_search/{}'.format(search_opts.search_ver)],
                         submit=True, cpu_threads=search_opts.n_omp_threads,
                         mem=mem, temp_mem=temp_mem)
                job_id_list.append(job_id)

                check_job_num += 1
                processing_time = 0.0
                check_batch = "{0}_{1}_a{2}_{3}".format(search_opts.bsd_row_num, presto_command,
                        search_opts.attempt + 1, check_job_num)
                commands = []
                bash_commands = []

                if not bash_job:
                    commands.append(add_database_function())
                commands.append('cd {0}{1}/{2}'.format(search_opts.work_dir, search_opts.sub_dir,
                                                       sub_sub_dir))

        #Check there is a command to run
        if len(bash_commands) > 0:
            if bash_job:
                #if it's a bash job make a file for it to run
                with open('{0}{1}/{2}.bash'.format(search_opts.work_dir,
                          search_opts.sub_dir, check_batch), "w") as srun_bash:
                    srun_bash.write(add_temp_database_function(check_job_num,
                                    search_opts.attempt + 1,
                                    n_omp_threads=search_opts.n_omp_threads))
                    for sc in bash_commands:
                        srun_bash.write("{}\n".format(sc))

                commands.append("srun --export=ALL -n 1 -c {0} bash {1}.bash".\
                                format(search_opts.n_omp_threads, check_batch))
                commands.append('search_database.py -c {0} -m m --file_location '
                                '{0}_temp_database_file_{1}_{2}.csv\n'.format(presto_command,
                                    search_opts.attempt + 1, check_job_num))
            else:
                commands = commands + bash_commands

            job_id = submit_slurm(check_batch, commands,
                     batch_dir="{0}{1}/batch".format(search_opts.work_dir, search_opts.sub_dir),
                     slurm_kwargs={"time": total_job_time_str},
                     nice=search_opts.nice,
                     module_list=['mwa_search/{}'.format(search_opts.search_ver)],
                     cpu_threads=search_opts.n_omp_threads, submit=True,
                     mem=mem, temp_mem=temp_mem)
            job_id_list.append(job_id)

        #Depsearch_opts.endancy job
        print("Waiting 5 sec to make sure to the dependancy script works")
        sleep(5)

        check_depend_batch = '{0}_dep_{1}_a{2}'.format(search_opts.bsd_row_num,
                                presto_command, search_opts.attempt +1)
        commands = []
        commands.append("{0} --attempt {1} -m c --table {2}".format(search_opts.relaunch_script,
                                              search_opts.attempt + 1, search_opts.table))

        submit_slurm(check_depend_batch, commands,
                     batch_dir="{0}{1}/batch".format(search_opts.work_dir, search_opts.sub_dir),
                     slurm_kwargs={"time": "20:00"},
                     nice=search_opts.nice,
                     submit=True, depend=job_id_list, depend_type="afterany",
                     module_list=['mwa_search/{}'.format(search_opts.search_ver)],
                     cpu_threads=search_opts.n_omp_threads)
    return




if __name__ == "__main__":
    mode_options = ["b", "r", "p", "t", "a", "f", "c", "w"]
    # Dictionary for choosing log-levels
    loglevels = dict(DEBUG=logging.DEBUG,
                     INFO=logging.INFO,
                     WARNING=logging.WARNING)
    parser = argparse.ArgumentParser(description="""
    Used to automate mass beamforming of MWA data and pulsar searches using the
    galaxy supercomputer (Ozstar coming soon).
    """, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-o', '--observation', type=str,
        help='The observation ID of the MWA observation')
    parser.add_argument('-c', '--channels', type=int, nargs=24, default=None,
        help='A list of the observations channel IDs for example "-c 109 110 111 '
             '112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 '
             '129 130 131 132". If this option is not used a metadata call will be '
             'used to find the channel IDs.')
    parser.add_argument('-v', '--mwa_search_version', type=str, default='master',
        help="The module version of mwa_search to use. Default: master")
    parser.add_argument("-L", "--loglvl", type=str, choices=loglevels.keys(), default="INFO",
        help="Logger verbosity level. Default: INFO")

    pointing_options = parser.add_argument_group('Pointing Options')
    pointing_options.add_argument('-p', '--pointing', type=str,
        help='The search_opts.pointing of the fits file to be searched. The RA and Dec '
             'seperated by _ in the format HH:MM:SS_+DD:MM:SS.')
    pointing_options.add_argument("--pulsar_file", default=None,
        help='Location of a file containting the pointings to be processed. '
             'Each line contains a pointing with an RA and Dec seperated by a '
             'space in the format "HH:MM:SS +DD:MM:SS". Can be made using grid.py.')
    pointing_options.add_argument('-i', '--incoh', action="store_true",
        help='Uses the incoh fits file location')

    beamform_options = parser.add_argument_group('Beamforming Options')
    beamform_options.add_argument("--DI_dir", default=None,
        help="Directory containing either Direction Independent Jones Matrices "
             "(as created by the RTS) or calibration_solution.bin as created by "
             "Andre Offringa's tools.[no default]")
    beamform_options.add_argument('-O', '--cal_obs', type=int,
        help="Observation ID of calibrator you want to process. Used to work out "
             "default DI_dir directory.")
    beamform_options.add_argument("-b", "--begin", type=int,
        help="First GPS time to process [no default]")
    beamform_options.add_argument("-e", "--end", type=int,
        help="Last GPS time to process [no default]")
    beamform_options.add_argument("-a", "--all", action="store_true",
        help="Perform on entire observation span. Use instead of -b & -e.")
    beamform_options.add_argument("--fits_dir", default=None,
        help="The directory containing the fits files. Only required if the fits "
             "files aren't in the standard vcs location.")

    search_options = parser.add_argument_group('Pulsar Search Options')
    search_options.add_argument("--search", action="store_true",
        help="Continue with the search pipeline after a successful beamforming "
             "check. Default False")
    search_options.add_argument("--relaunch", action="store_true",
        help="Will rerun the search pipeline even when no beamforming is required. "
             "Otherwise assumes that if the beamforming is complete then a search "
             "has already been performed.")
    search_options.add_argument('--code_comment', type=str,
        help='A comment describing the purpose of this search. eg testing')
    search_options.add_argument("--attempt", type=int, default=0,
        help="The number of attempts that a script has made. Default is 1.")
    search_options.add_argument('-w', '--work_dir', type=str, default=None,
        help='Work directory. Default: {}'.format(DEFAULT_WORK_DIR))
    search_options.add_argument('-s', '--sub_dir', type=str,
        help='Used by the program to keep track of the sub directory its using')
    search_options.add_argument('-r', '--bsd_row_num', type=int,
        help='Database row reference number for keeping track of the scripts.')
    search_options.add_argument('--dm_min', type=float, default = 0.0,
        help='DM max searched. Default 1')
    search_options.add_argument('--dm_max', type=float, default = 250.0,
        help='DM max searched. Default 250')
    search_options.add_argument('--pulsar', type=str,
        help="Used to search for a known pulsar by inputing it's Jname. The code "
             "then looks within 1 DM and 15%% of the pulsar's period.")
    search_options.add_argument('-m', '--mode', type=str, default="b",
        help=textwrap.dedent('''Modes used by the pipeline do indicate which step in the process it is up to. The default is beamform mode.
"b" Checks the fits file directory to decide if all files are there or if splicing or beamforming is required.
"r" Uses PRESTO's rfifind to create a RFI mask.
"p" Uses PRESTO's prepdata to dedisperses the the fits files into .dat files.
"t" Uses PRESTO's realfft to do a fast fourier transform on the .dat files.
"a" Uses PRESTO's accelsearch to search for periodic signals.
"f" Uses PRESTO's prepfold to fold on any significant candidates.
"w" Wraps up the pipeline by moving significant pulse profiles to it's own file for inspection.
"c" Check script that runs after each step to make sure no jobs have failed.'''))
    search_options.add_argument("--table", type=str,
        help="The table name for the search database to check, similar to the commands used.")
    search_options.add_argument("--csc", action="store_true",
        help="If option used will check if the pointing is being stored in Galaxy's "
             "cold storage HSM.")
    args=parser.parse_args()
    comp_config = config.load_config_file()

    # set up the logger for stand-alone execution
    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    ch.setLevel(loglevels[args.loglvl])
    formatter = logging.Formatter('%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)


    #argument parsing
    if args.mode not in mode_options:
        print("Unrecognised mode [{0}]. Please use a mode from {1}. Exiting".format(
               args.mode, mode_options))
        quit()

    if not args.observation:
        print("No observation id given. Please supply one using -o. Exiting")
        quit()

    if args.pulsar and not args.pointing:
        #if no pointing given grab it from psrcat
        import find_pulsar_in_obs as fpio
        temp = fpio.get_psrcat_ra_dec(pulsar_list=[args.pulsar])
        temp = fpio.format_ra_dec(temp, ra_col = 1, dec_col = 2)
        jname, raj, decj = temp[0]
        args.pointing = '{0}_{1}'.format(raj, decj)
    elif not (args.pointing or args.pulsar_file or args.incoh):
        print("No pointing option supplied. Please either use -p or --pulsar_file. Exiting")
        quit()

    if args.pulsar:
        #only search aroung the pulsar DM
        dm, p = get_pulsar_dm_p(args.pulsar)
        args.dm_min = float(dm) - 2.0
        if args.dm_min < 1.0:
            args.dm_min = 1.0
        args.dm_max = float(dm) + 2.0
        print("Searching DMs from {0} to {1}".format(args.dm_min, args.dm_max))

    if args.mode == "b":
        if args.incoh:
            args.cal_obs = None
        elif not args.DI_dir:
            if args.cal_obs:
                args.DI_dir = "{0}{1}/cal/{2}/rts/".format(comp_config['base_product_dir'],
                                                          args.observation,
                                                          args.cal_obs)
                print("No DI_dir given so assuming {0} is the directory".format(args.DI_dir))
            else:
                print("Please use either --DI_dir to explicitly or use --cal_obs if "
                      "the DI Jones matrices are in the standard directory. Exiting.")
                quit()
        #check begining end times
        if args.all and (args.begin or args.end):
            print("Please specify EITHER (-b,-e) OR -a")
            quit()
        elif args.all:
            args.begin, args.end = meta.obs_max_min(args.observation)
        elif not (args.all or args.begin or args.end):
            print("Please specify the beginning and end times using (-b,-e) OR -a. Exiting")
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
    else:
        work_dir = DEFAULT_WORK_DIR
    #work_dir = os.environ['BLINDSEARCH_WORK_DIR']
    if not work_dir.endswith('/'):
        work_dir = "{0}/".format(work_dir)

    if args.sub_dir:
        sub_dir = args.sub_dir
    elif not args.mode == 'b':
        sub_dir = pointing + "/" + obsid + "/"
    else:
        sub_dir = None

    fits_dir_base = '{0}{1}/'.format(comp_config['base_product_dir'], obsid)

    if args.fits_dir:
        pointing_dir = args.fits_dir
    elif args.pointing:
        if args.incoh:
            pointing_dir = '{0}/incoh/'.format(fits_dir_base)
        else:
            pointing_dir = '{0}/pointings/{1}/'.format(fits_dir_base, args.pointing)
    else:
        pointing_dir = None
    #get channels
    args.channels = get_channels(obsid, channels=args.channels)

    #Base launch of this code (everything except mode and dmfile int)
    relaunch_script = "mwa_search_pipeline.py -o " + str(obsid) + " -w " + work_dir

    if not args.mode =='b':
        relaunch_script +=  " -s " + str(sub_dir)
    if args.DI_dir:
        relaunch_script += " --DI_dir="+args.DI_dir
    if args.bsd_row_num:
        relaunch_script += " -r " + str(args.bsd_row_num)
    if not args.pulsar == None:
        relaunch_script += " --pulsar " + args.pulsar
    if args.pointing:
        relaunch_script += " -p " + str(pointing)
    if args.begin and args.end:
        relaunch_script += " -b " + str(args.begin) + " -e " + str(args.end)
    if args.search and args.mode == 'b':
        relaunch_script += " --search"
    if args.dm_max:
        relaunch_script += " --dm_max " + str(args.dm_max)
    if args.dm_min:
        relaunch_script += " --dm_min " + str(args.dm_min)
    if args.incoh:
        relaunch_script +=  " --incoh "
    if args.mwa_search_version:
        relaunch_script +=  " -v " + str(args.mwa_search_version)
    if args.fits_dir:
        relaunch_script +=  " --fits_dir " + str(args.fits_dir)
    if args.channels:
        relaunch_script +=  " --channels"
        for ch in args.channels:
            relaunch_script +=  " {0}".format(ch)
    if args.cal_obs:
        relaunch_script +=  " -O " + str(args.cal_obs)
    if args.relaunch:
        relaunch_script +=  " --relaunch"


    search_opts = search_options_class(obsid, pointing=pointing, cal_id=args.cal_obs,
                 begin=args.begin, end=args.end, channels=args.channels, incoh=args.incoh,
                 args=args, work_dir=work_dir, sub_dir=sub_dir, fits_dir_base=fits_dir_base,
                 pointing_dir=pointing_dir, dm_min=args.dm_min, dm_max=args.dm_max,
                 relaunch_script=relaunch_script, search=args.search,
                 bsd_row_num=args.bsd_row_num, cold_storage_check=args.csc,
                 table=args.table, attempt=args.attempt, search_ver=args.mwa_search_version,
                 DI_dir=args.DI_dir)

    #work out start and stop times for beamforming
    if args.mode == "b" or args.mode == None:
        if args.pulsar_file:
            with open(args.pulsar_file) as f:
                pointing_list_raw = f.readlines()
            #remove \n
            pointing_list = []
            for point in pointing_list_raw:
                pointing_list.append(point.strip())
        elif args.pointing:
            pointing_list = [args.pointing.replace("_"," ")]
        """
        elif args.fits_dir:
            #gives it a pointing name based on the last dir in the fits
            #directory as that's normaly a pulsar id or pointing
            if args.fits_dir.endswith("/"):
                args.fits_dir = args.fits_dir[:-1]
            pointing_list = [args.fits_dir.split("/")[-1]]
        """

        #If in search mode start up the database entry
        if args.code_comment:
            code_comment = args.code_comment
        elif args.search and not args.bsd_row_num:
            code_comment = input("Please write a comment describing the purpose "
                                 "of this search. eg testing: ")
        else:
            code_comment = None
        if args.pulsar_file and code_comment is not None:
            code_comment += " (using: {0}) ".format(args.pulsar_file)

        if args.search:
            #if searching use summed polarisation mode to save on psrfits data size
            summed=True
        else:
            summed=False

        #pulsar check parsing
        if args.pulsar is None:
            pulsar_list_list=None
        else:
            pulsar_list_list = [[args.pulsar]]
        beamform(search_opts, pointing_list, code_comment=code_comment,
                 pulsar_list_list=pulsar_list_list, summed=summed,
                 relaunch=args.relaunch)
    elif args.mode == "c":
        error_check(search_opts)

    elif args.mode == "r":
        rfifind(search_opts)
    elif args.mode == "p":
        prepdata(search_opts)
    elif args.mode == "t":
        sort_fft(search_opts)
    elif args.mode == "a":
        accel(search_opts)
    elif args.mode == "f":
        fold(search_opts)
    elif args.mode == "w":
        wrap_up(search_opts)
