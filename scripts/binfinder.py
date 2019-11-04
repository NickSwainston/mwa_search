#!/usr/bin/env python3

import os
import glob
import logging
import argparse
import sys
import config

import data_process_pipeline
from job_submit import submit_slurm
import plotting_toolkit
import find_pulsar_in_obs as fpio
import check_known_pulsars as checks
import sn_flux_est as snfe
import file_maxmin
import psrqpy
logger = logging.getLogger(__name__)

#get ATNF db location
try:
    ATNF_LOC = os.environ['PSRCAT_FILE']
except KeyError:
    logger.warn("ATNF database could not be loaded on disk. This may lead to a connection failure")
    ATNF_LOC = None

#----------------------------------------------------------------------
def copy_to_product_dir(pulsar, pointing_dir, obsid):
    """
    Copies the files fitting the input specifications from a prepfold operation to a data product directory
    Parameters:
    -----------
    bins: int
        The number of bins from the prepfold operation
    pulsar: string
        The J name of the pulsar
    pointing_dir: string
        The pointing directory holding the data products
    obsid: int
        The obsid of the observation
    """
    comp_config = config.load_config_file()
    base_dir = comp_config['base_product_dir']
    product_dir = os.path.join(base_dir, obsid, "data_products", pulsar)
    all_bins = find_bins_in_dir(pointing_dir)
    bin_limit = bin_sampling_limit(pulsar)
    copy_bins=[]
    if bin_limit<50:
        copy_bins.append(all_bins[0])
    else:
        copy_bins.append( min(all_bins[1:]))
    if bin_limit<100:
        copy_bins.append(50)
    else:
        copy_bins.append(100)

    data_products=[]
    for bins in copy_bins:
        files = glob.glob(os.path.join(pointing_dir, "*{0}*bins*{1}*".format(bins, pulsar[1:])))
        for afile in files:
            data_products.append(afile)

    for product in data_products:
        data_process_pipeline.copy_data(product, product_dir)

#----------------------------------------------------------------------
def add_prepfold_to_commands(commands, pointing, pulsar, obsid, beg, end, nbins,\
                            use_mask=True, ntimechunk=120, dmstep=1, period_search_n=1):
    """
    Adds prepfold commands to a list

    Parameters:
    -----------
    commands: list
        A list of commands. Can be empty, this list will be appended to by the function
    pointing: string
        The directory to work in. Typically the pointing directory.
    puslar: string 
        The J name of the pulsar
    obisd: int
        The ID of the observation
    use_mask: boolean
        Whether or not to use a mask if it exists
    start: float
        OPTIONAL - The normalized time that the pulsar enters the beam. Probably from snfe.pulsar_beam_coverage()Default=None
    end: float
        OPTOINAL - The normalized time that the pulsar exits the beam. Probably from snfe.pulsar_beam_coverage(). Default=None
    nbins: int
        OPTIONAL - The number of bins to fold over. Default=100
    ntimechunk: int
        OPTIONAL - The ntimechunk option for prepfold. Default=120
    dmstep: float
        OPTIONAL - The dmstep option for prepfold. Default=1
    period_search_n: float
        OPTIONAL - The period_search_n option for perpfold. Default=1

    Returns: 
    --------
    commands: list
        The commands list that was input with the prepfold commands appended
    """
    #find the beginning and end of the pulsar's beam coverage for this obs
    enter, exit = snfe.pulsar_beam_coverage(obsid, pulsar, beg=beg, end=end)
    logger.info("start and end of pulsar beam coverage for beginng time: {0} and end time: {1}:"\
                ":{2}, {3}".format(beg, end, enter, exit))
    if enter is None or exit is None:
        logger.error("pulsar is not in beam for any of the on-disk files. Ending...")
        sys.exit(1)

    comp_config = config.load_config_file()
    #Figure out whether or not to input a mask
    if use_mask == True:
        check_mask = glob.glob(os.path.join(comp_config['base_product_dir'], obsid, "incoh", "*mask"))
        if check_mask:
            mask = "-mask " + check_mask[0]
        else:
            mask = ""
    else:
        mask=""

    #make the prepfold command
    constants = "-pstep 1 -pdstep 2 -ndmfact 1 -noxwin -nosearch -runavg -noclip -nsub 256 1*fits "
    variables = "-o {0}_{1}_bins ".format(obsid, nbins)
    variables += mask
    variables += "-n {0} ".format(nbins)
    variables += "-start {0} -end {1} ".format(enter, exit)
    variables += "-dmstep {0} ".format(dmstep)
    variables += "-npart {0} ".format(ntimechunk)
    variables += "-npfact {0} ".format(period_search_n)

    #load presto module here because it uses python 2
    commands.append('cd {0}'.format(pointing))
    commands.append('echo "Folding on known pulsar {0}"'.format(pulsar))
    commands.append('psrcat -e {0} > {0}.eph'.format(pulsar))
    commands.append("sed -i '/UNITS           TCB/d' {0}.eph".format(pulsar))
    commands.append("prepfold -timing {0}.eph {1} {2}"\
                    .format(pulsar, variables, constants))
    commands.append('errorcode=$?')
    commands.append('pulsar={}'.format(pulsar[1:]))

    #Some old ephems don't have the correct ra and dec formating and
    #causes an error with -timing but not -psr
    commands.append('if [ "$errorcode" != "0" ]; then')
    commands.append('   echo "Folding using the -psr option"')
    commands.append('   prepfold -psr {0} {1} {2}'\
                    .format(pulsar, variables, constants))
    commands.append('   pulsar={}'.format(pulsar))
    commands.append('fi')
    commands.append('rm {0}.eph'.format(pulsar))

    return commands

#----------------------------------------------------------------------
def bestprof_info(filename):
    """
    Finds various information on a .bestprof file

    Parameters:
    filename: string
        The path of the bestprof file
    Returns:
    info_dict: dictionary
        A dictionary consisting of the following:
        obsid: int
            The ID of the observation
        puslar: string
            The J name of the pulsar
        nbins: int
            The number of bins used to fold this profile
        chi: float
            The reduced Chi squared value of the fold
        sn: float
            The signal to noise ratio of the fold
        dm: float
            The pulsar's dispersion measure
        period: float
            The pulsar's period
        period_error: float
            The error in the pulsar's period measurement
    """
    #open the file and read the info into a dictionary
    info_dict = {}
    f = open(filename, "r")
    lines = f.read()
    lines = lines.split("\n")
    #info:
    info_dict["obsid"] = int(lines[0].split()[4].split("_")[0])
    info_dict["pulsar"] = lines[1].split()[3].split("_")[1]
    info_dict["nbins"] = int(lines[9].split()[4])
    info_dict["chi"] = float(lines[12].split()[4])
    info_dict["sn"] = float(lines[13].split()[4][2:])
    info_dict["dm"] = float(lines[14].split()[4])
    info_dict["period"] = float(lines[15].split()[4]) #in ms
    info_dict["period_error"] = float(lines[15].split()[6])
    f.close()
    return info_dict

#----------------------------------------------------------------------
def bin_sampling_limit(pulsar, sampling_rate=1e-4):
    """
    Finds the sampling limit of the input pulsar in units of number of bins

    Parameters:
    -----------
    puslar: string
        The J name of the pulsar to check the sampling limit for
    sampling_rate: float
        OPTIONAL - the sampling rate of the instrument. Default=1e-4 (MWA VCS)

    Returns:
    bin_lim: int
        The highest number of bins that can be logically folded on
    """
    query = psrqpy.QueryATNF(params=["P0"], psrs=[pulsar], loadfromdb=ATNF_LOC).pandas
    period = query["P0"][0]
    bin_lim = int(period/sampling_rate + 1) #the +1 is to round the limit up every time
    logger.debug("Bin limit: {0}".format(bin_lim))
    return bin_lim


#----------------------------------------------------------------------
def submit_to_db_and_continue(run_params, best_bins):
    """
    Submits the best fold profile to the pulsar database. Will also submit .ppps, .png and .pfd

    Parameters:
    -----------
    run_params: object
        The run_params object defined in data_process_pipeline
    best_bins: int
    """
    #Add path to filenames for submit script
    cwd = os.getcwd()

    ppps = glob.glob("*{0}_bins*{1}*.pfd.ps".format(best_bins, run_params.pulsar[1:]))[0]
    ppps = os.path.join(cwd, ppps)
    bestprof = glob.glob("*{0}_bins*{1}*.pfd.bestprof".format(best_bins, run_params.pulsar[1:]))[0]
    bestprof = os.path.join(cwd, bestprof)
    png = glob.glob("*{0}_bins*{1}*.png".format(best_bins, run_params.pulsar[1:]))[0]
    png = os.path.join(cwd, png)
    pfd = glob.glob("*{0}_bins*{1}*.pfd".format(best_bins, run_params.pulsar[1:]))[0]
    pfd = os.path.join(cwd, pfd)
    
    commands = []
    commands.append("echo 'Submitting profile to database with {} bins'".format(best_bins))
    commands.append('submit_to_database.py -o {0} --cal_id {1} -p {2} --bestprof {3} --ppps {4}'\
    .format(run_params.obsid, run_params.cal_id, run_params.pulsar, bestprof, ppps))
    commands.append('echo "submitted profile to database: {0}"'.format(bestprof))


    #Make a nice plot
    plotting_toolkit.plot_bestprof(os.path.join(run_params.pointing_dir, bestprof),\
                                    out_dir=run_params.pointing_dir)
    #copy all data to product directoy for easy viewing
    copy_to_product_dir(run_params.pulsar, run_params.pointing_dir, run_params.obsid) 

    bin_lim = bin_sampling_limit(run_params.pulsar)
    if bin_lim>100:
        b_standard=100
    else:
        b_standard=50
   
    if best_bins != b_standard: 
        #do the same for 100/50 bin profiles depending on whether this is an msp or not
        ppps = glob.glob("*{0}_bins*{1}*.pfd.ps".format(b_standard, run_params.pulsar[1:]))[0]
        ppps = os.path.join(cwd, ppps)
        bestprof = glob.glob("*{0}_bins*{1}*.pfd.bestprof".format(b_standard, run_params.pulsar[1:]))[0]
        bestprof = os.path.join(cwd, bestprof)
        png = glob.glob("*{0}_bins*{1}*.png".format(b_standard, run_params.pulsar[1:]))[0]
        png = os.path.join(cwd, png)
        pfd = glob.glob("*{0}_bins*{1}*.pfd".format(b_standard, run_params.pulsar[1:]))[0]
        pfd = os.path.join(cwd, pfd)
 
        commands = []
        commands.append("echo 'Submitting profile to database with {} bins'".format(b_standard))
        commands.append('submit_to_database.py -o {0} --cal_id {1} -p {2} --bestprof {3} --ppps {4}'\
        .format(run_params.obsid, run_params.cal_id, run_params.pulsar, bestprof, ppps))
        commands.append('echo "submitted profile to database: {0}"'.format(bestprof))

    #submit job
    name = "Submit_db_{0}_{1}".format(run_params.pulsar, run_params.obsid)
    comp_config = config.load_config_file()
    batch_dir = os.path.join(comp_config['base_product_dir'], run_params.obsid, "batch")
    logger.info("Submitting submission script for profile: {0}".format(bestprof))
    logger.info("Job name: {}".format(name))

    submit_slurm(name, commands,\
                 batch_dir=batch_dir,\
                 slurm_kwargs={"time": "00:05:00"},\
                 module_list=['mwa_search/{0}'.format(run_params.mwa_search)],\
                 submit=True, vcstools_version="{0}".format(run_params.vcs_tools))

    #Run stokes fold
    commands = []
    commands.append("data_process_pipeline.py -d {0} -O {1} -p {2} -o {3} -b {4} -L {5}\
                    --mwa_search {6} --vcs_tools {7} -m s"\
                    .format(run_params.pointing_dir, run_params.cal_id, run_params.pulsar,\
                    run_params.obsid, best_bins, run_params.loglvl, run_params.mwa_search,\
                    run_params.vcs_tools))

    name = "dpp_stokes_{0}_{1}".format(run_params.pulsar, run_params.obsid)
    comp_config = config.load_config_file()
    batch_dir = os.path.join(comp_config['base_product_dir'], run_params.obsid, "batch")
    logger.info("Submitting Stokes Fold script")
    logger.info("Job Name: {}".format(name))
    
    submit_slurm(name, commands,\
                 batch_dir=batch_dir,\
                 slurm_kwargs={"time": "00:05:00"},\
                 module_list=['mwa_search/{0}'.format(run_params.mwa_search)],\
                 submit=True, vcstools_version="{0}".format(run_params.vcs_tools))

#----------------------------------------------------------------------
def get_best_profile_in_dir(pointing_dir, pulsar):
    """
    Finds the profile with the highest bins with acceptable SN and Chi values in a directory
    
    Parameters:
    -----------
    potining_dir: string
        The poitning directory for the function to look through
    pulsar: string
        The J name of the pulsar
    
    Returns:
    --------
    prof_name: string
        The name of the best profile found in the directory
    """
    #find all of the relevant bestprof profiles in the pointing directory
    bestprof_names = glob.glob("*bins*{0}*.bestprof".format(pulsar[1:]))
    if len(bestprof_names)==0:
        logger.error("No bestprofs found in directory! Exiting")
        sys.exit(1)

    #throw all of the information from each bestprof into an array
    bin_order = []
    sn_order = []
    chi_order = []
    for prof in bestprof_names:
        prof_info = bestprof_info(filename=prof)
        bin_order.append(prof_info["nbins"])
        sn_order.append(prof_info["sn"])
        chi_order.append(prof_info["chi"])
    bin_order, sn_order, chi_order = zip(*sorted(zip(bin_order, sn_order, chi_order)))
    bin_order = bin_order[::-1]
    sn_order = sn_order[::-1]
    chi_order = chi_order[::-1]

    #now find the one with the most bins that meet the sn and chi conditions
    best_i = None
    bin_lim = bin_sampling_limit(pulsar)
    for i in range(len(bin_order)):
        if bin_order[i]<=bin_lim: #only consider profiles where the number of bins used is lower than the bin upper limit
            if sn_order[i]>10. and chi_order[i]>4. == True:
                best_i = i
                break

    if best_i is None:
        logger.info("No profiles fit the threshold parameter")
        return None
    else:
        logger.info("Adequate profile found with {0} bins".format(bin_order[best_i]))
        prof_name = glob.glob("*{0}_bins*{1}*.bestprof".format(bin_order[best_i], pulsar[1:]))[0]
        return prof_name

#----------------------------------------------------------------------
def sn_chi_test(bestprof, sn_thresh=10., chi_thresh=4.):
    """
    Checks whether the input signal to noise ratio and chi values are good or bad

    Parameters:
    -----------
    bestprof: string
        The name of the bestprof file to check
    sn_thresh: float
        OPTIONAL - The lower limit on what is an acceptable sn. Default=10.
    ch_thresh: float
        OPTIONAL - The lower limit on what is an acceptable reduced chi squared value. Defaultl=4.
    
    Returns:
    --------
    test: boolean
        Whether or not the input sn and chi combination is acceptable
    """
    test = False
    info_dict = bestprof_info(bestprof)
    sn = info_dict["sn"]
    chi = info_dict["chi"]
    if sn >= sn_thresh and chi >= chi_thresh:
        test = True
    elif sn == 0. and chi >= chi_thresh:
        test = True
    return test

#----------------------------------------------------------------------
def submit_multiple_pointings(run_params):
    """
    Submits many pointings on either 50 or 100 bins depending on the pulsar's period. Launches the next step of binfinder

    Parameters:
    -----------
    run_params: object
        The run_params object defined in data_process_pipeline
    """
    job_ids = []
    comp_config=config.load_config_file()

    #Check number of bins to use
    bin_limit=bin_sampling_limit(run_params.pulsar)
    if bin_limit<100:
        nbins=50
    else:
        nbins=100
   
    logger.info("Submitting multiple prepfold jobs:")
    logger.info("Pulsar name: {}".format(run_params.pulsar))
    logger.info("Number of bins to fold on: {}".format(nbins))
 
    #submit a job for each pointing
    for i, pointing in enumerate(run_params.pointing_dir):
        
        #os.chdir(pointing)
        #create slurm job:
        commands = []
        commands = add_prepfold_to_commands(commands, pointing, run_params.pulsar, run_params.obsid,\
                    run_params.beg, run_params.end, nbins)
        name = "bf_multi_{0}_{1}_{2}_bins_{3}".format(run_params.pulsar, run_params.obsid, nbins, i)
        batch_dir = os.path.join(comp_config['base_product_dir'], run_params.obsid, "batch")
        logger.info("Submitting pointing: {0}".format(pointing))
        logger.info("Job name: {}".format(name))
        myid = submit_slurm(name, commands,\
                    batch_dir=batch_dir,\
                    slurm_kwargs={"time": "2:00:00"},\
                    module_list=['mwa_search/{0}'.format(run_params.mwa_search),\
                                'presto/no-python'],\
                    submit=True, vcstools_version="{0}".format(run_params.vcs_tools))


        job_ids.append(myid)

    p = ""
    for pointing in run_params.pointing_dir:
        p += p + " "

    commands=[]
    commands.append("binfinder.py -d {0} -O {1} -p {2} -o {3} -L {4} {5} --vcs_tools {6}\
                    --mwa_search {7} -p {8} -b {9} -e {10}"\
                    .format(p, run_params.cal_id, run_params.pulsar, run_params.obsid,\
                    run_params.loglvl, run_params.vcs_tools, run_params.mwa_search, run_params.pulsar,\
                    run_params.beg, run_params.end))

    name="bf_post_multi_{0}".format(run_params.pulsar)
    batch_dir = os.path.join(comp_config['base_product_dir'], run_params.obsid, "batch")
    logger.info("Submitting post-prepfold binfinder:")
    logger.info("Job name: {}".format(name))
    myid = submit_slurm(name, commands,\
            batch_dir=batch_dir,\
            slurm_kwargs={"time": "00:30:00"},\
            module_list=['mwa_search/{0}'.format(run_params.mwa_search),\
                        "presto/no-python"],\
            submit=True, depend=job_ids, depend_type="afterany",\
            vcstools_version="master")


#----------------------------------------------------------------------
def submit_prepfold(run_params, nbins):
    """
    Submits a prepfold job for the given parameters
    
    Parameters:
    -----------
    run_params: object
        The run_params object defined by data_process_pipeline
    nbins: int
        The number of bins to fold on 
    """
    commands = []
    commands.append("echo '############### Prepfolding on {} bins ###############'".format(nbins))
    commands = add_prepfold_to_commands(commands, run_params.pointing_dir, run_params.pulsar, run_params.obsid, run_params.beg, run_params.end, nbins)
    commands.append("echo '############### Relaunching binfinder script ###############'" )

    #binfinder relaunch:
    commands.append("binfinder.py -d {0} -O {1} -p {2} -o {3} -L {4} --vcs_tools {5}\
                    --mwa_search {6} -b {7} -e {8}"\
                    .format(run_params.pointing_dir, run_params.cal_id, run_params.pulsar,\
                    run_params.obsid, run_params.loglvl, run_params.vcs_tools, run_params.mwa_search,\
                    run_params.beg, run_params.end))

    comp_config = config.load_config_file()
    batch_dir = os.path.join(comp_config['base_product_dir'], run_params.obsid, "batch")
    name = "bf_{0}_{1}_{2}_bins".format(run_params.pulsar, run_params.obsid, nbins)   
 
    logger.info("Submitting prepfold and resubmission job:")
    logger.info("Pointing directory: {}".format(run_params.pointing_dir))
    logger.info("Pulsar name: {}".format(run_params.pulsar))
    logger.info("Number of bins to fold on: {}".format(nbins))
    logger.info("Job name: {}".format(name))

    submit_slurm(name, commands,\
                batch_dir=batch_dir,\
                slurm_kwargs={"time": "2:00:00"},\
                module_list=['mwa_search/{0}'.format(run_params.mwa_search),\
                            'presto/no-python'],\
                submit=True, vcstools_version="{0}".format(run_params.vcs_tools))



#----------------------------------------------------------------------
def find_best_pointing(run_params, nbins):

    """
    Finds the pointing directory with the highest S/N out of those given in the pointing directory

    Parameters:
    -----------
    run_params: object
        The run_params object defined in data_process_pipeline
    nbins: int
        The number of bins that were folded on that the function will look through
    """

    bestprof_info_list = []
    for pointing in run_params.pointing_dir:
        os.chdir(pointing)
        logger.info("searching directory: {0}".format(pointing))
        prof_name = glob.glob("*{0}_bins*{1}*.bestprof".format(nbins, run_params.pulsar[1:]))[0]
        bestprof_info_list.append(bestprof_info(filename=prof_name))

    #now we loop through all the info and find the best one
    best_sn = 0.0
    best_i = -1
    for i, info_dict in enumerate(bestprof_info_list):
        if info_dict["chi"]>=4.0 and info_dict["sn"]>best_sn:
            best_sn = info_dict["sn"]
            best_i = i

    if best_i<0 and best_sn<run_params.threshold:
        logger.info("No pulsar found in pointings. Exiting...")
        sys.exit(0)
    else:
        logger.info("Pulsar found in pointings. Running binfinder script on pointing: {0}"\
                    .format(run_params.pointing_dir[best_i]))

    #submit the next iteration
    run_params.set_pointing_dir(run_params.pointing_dir[best_i])
    next_bins = how_many_bins_next(run_params.pulsar, run_params.pointing_dir)
    submit_prepfold(run_params, nbins=next_bins)

#----------------------------------------------------------------------
def how_many_bins_next(pulsar, directory):
    """
    Works out how many bins to use for the next fold operation. 
    The logic fucntion is based only on the pulsar's period and what folds have been done already.
    NOTE: will return None if no further fold operations should be done
    
    Parameters:
    -----------
    pulsar: string
        The J name of the pulsar
    directory: string
        The name of the pointing directory

    Returns:
    --------
    next_bins: int
        The number of bins to fold on next
    """
    bins_in_dir = find_bins_in_dir(directory)
    bin_limit = bin_sampling_limit(pulsar)
    
    #if msp
    if bin_limit<100:
        if 50 not in bins_in_dir:
            return 50
        elif bin_limit not in bins_in_dir:
            return bin_limit
        else:
            return None
    #high period pulsar
    elif bin_limit<1024:
        if 100 not in bins_in_dir:
            return 100
        elif bin_limit not in bins_in_dir:
            return bin_limit
        else:
            next_bins =int(min(bins_in_dir[1:])/2)
    #regular period pulsar
    else:
        if 100 not in bins_in_dir:
            return 100
        elif 1024 not in bins_in_dir:
            return 1024
        else:
            next_bins = int(min(bins_in_dir[1:])/2)
     
    #Check to see if this is an acceptable number
    if next_bins<100:
        next_bins = None

    return next_bins
#----------------------------------------------------------------------
def find_bins_in_dir(directory):
    """
    Works out what folds have been executed in a given directory already. 
    Assumes all files are in the format 'obsid_bins_...'

    Parameters:
    -----------
    directory: string
        The directory to look in

    Returns:
    --------
    bins_in_dir: list
        A list of ints of the bin numbers found in the directory
    """
    fold_files=glob.glob("*.pfd.bestprof")
    bins_in_dir=[]
    for bestprof in fold_files:
        bins_in_dir.append(int(bestprof.split("_")[1]))
    bins_in_dir.sort()
    return bins_in_dir
 
#----------------------------------------------------------------------
def work_out_what_to_do(run_params):
    """
    A logic structure that decides what to do next in the binfinder pipeline

    Parameters:
    -----------
    run_params: object
        The run_params object defined by data_proces_pipeline
    """
    #Multiple pointings?
    if isinstance(run_params.pointing_dir, list):
        #Have these directories been folded on?
        bins_in_dir = find_bins_in_dir(run_params.pointing_dir[0])
        if len(bins_in_dir)==0:
            #No folds done yet, fold on all pointings
            submit_multiple_pointings(run_params)      
        else:
            #Check the pointings and fold on the best one, if any
            find_best_pointing(run_params, nbins=bins_in_dir[0])

    elif isinstance(run_params.pointing_dir, str):
        #Get some info on where we're at
        bin_limit = bin_sampling_limit(run_params.pulsar)
        bins_in_dir = find_bins_in_dir(run_params.pointing_dir)
        if len(bins_in_dir)==0:
            #No folds done
            next_bins = how_many_bins_next(run_params.pulsar, run_params.pointing_dir)
            submit_prepfold(run_params, next_bins)

        elif len(bins_in_dir)==1:
            #only 100/50 bin fold in directory
            test_dir = os.path.join(run_params.pointing_dir, "*bestprof")
            test_dir = glob.glob(test_dir)[0]
            test = sn_chi_test(test_dir)    
            if test==True:
                next_bins = how_many_bins_next(run_params.pulsar, run_params.pointing_dir) 
                submit_prepfold(run_params, next_bins)
                return
            else:
                logger.info("No pulsar found in initial pointing. Exiting...")
                return

        else: #more than one fold done
           #check the last fold done 
            if bin_limit<50:
                last_fold_bins = bins_in_dir[0]
            else:
                last_fold_bins = min(bins_in_dir[1:])
            test_dir = os.path.join(run_params.pointing_dir,"*{0}_bins*bestprof".format(last_fold_bins))
            test_dir = glob.glob(test_dir)[0]
            test = sn_chi_test(test_dir)
            if test==False:
                next_bins = how_many_bins_next(run_params.pulsar, run_params.pointing_dir)
                if next_bins is not None:
                    submit_prepfold(run_params, next_bins)
                else:
                    logger.info("Minimum bin count hit")
                    if bin_limit<50: 
                        sub_bins = bins_in_dir[1]
                    else:
                        sub_bins = bins_in_dir[0]
                    submit_to_db_and_continue(run_params, sub_bins)
            else:
                logger.info("Pulsar detected: {0}".format(run_params.pulsar))
                logger.info("Submitting to database with {0} bins".format(last_fold_bins))
                submit_to_db_and_continue(run_params, last_fold_bins)
                return


#----------------------------------------------------------------------
if __name__ == '__main__':
    #dictionary for choosing log-levels
    loglevels = dict(DEBUG=logging.DEBUG,
                     INFO=logging.INFO,
                     WARNING=logging.WARNING,
                     ERROR = logging.ERROR)

    #Arguments
    parser = argparse.ArgumentParser(description="A script that handles pulsar folding operations")

    required = parser.add_argument_group("Required Inputs:")
    required.add_argument("-d", "--pointing_dir", action="store", nargs="+", help="Pointing directory(s) that contains the spliced fits files.")
    required.add_argument("-o", "--obsid", type=str, help="The observation ID")
    required.add_argument("-O", "--cal_id", type=str, help="The Obs ID of the calibrator")
    required.add_argument("-p", "--pulsar", type=str, help="The name of the pulsar. eg. J2241-5236")
    required.add_argument("-b", "--beg", type=int, help="The beginning of the observation. Will try to find if unsupplied")
    required.add_argument("-e", "--end", type=int, help="The end of the observation. Will try to find if unsupplied")
    

    other = parser.add_argument_group("Other Options:")
    other.add_argument("-t", "--threshold", type=float, default=10.0, help="The signal to noise threshold to stop at. Default = 10.0")
    other.add_argument("-L", "--loglvl", type=str, default="INFO", help="Logger verbosity level. Default: INFO", choices=loglevels.keys())
    other.add_argument("-S", "--stop", action="store_true", help="Use this tag to tell binfinder to launch the next step in the data processing pipleline when finished")
    other.add_argument("--mwa_search", type=str, default="master", help="The version of mwa_search to use. Default: master")
    other.add_argument("--vcs_tools", type=str, default="master", help="The version of vcs_tools to use. Default: master")

    args = parser.parse_args()
    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    ch.setLevel(loglevels[args.loglvl])
    formatter = logging.Formatter('%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.propagate = False

    #Checking required inputs
    if args.pointing_dir == None:
        logger.error("No pointing directory supplied. Please specify the pointing directory path and rerun")
        sys.exit(1)
    elif args.cal_id == None:
        logger.error("No calibrator ID supplied. Please input a cal ID and rerun")
        sys.exit(1)
    elif args.pulsar == None:
        logger.error("No pulsar name supplied. Please input a pulsar and rerun")
        sys.exit(1)


    run_params = data_process_pipeline.run_params_class\
                    (args.pointing_dir, args.cal_id, pulsar=args.pulsar,obsid=args.obsid,\
                    threshold=args.threshold, stop=args.stop,loglvl=args.loglvl,\
                    mwa_search=args.mwa_search, vcs_tools=args.vcs_tools,\
                    beg=args.beg, end=args.end)


    #NOTE: for some reason, you need to run prepfold from the directory it outputs to if you want it to properly make an image. The script will make this work regardless by using os.chdir
    logger.info("Pointing Dir: {} sadfasd".format(run_params.pointing_dir))
    if isinstance(run_params.pointing_dir, str):
        os.chdir(run_params.pointing_dir)
        

    #Try to find the beginning and end using ondisk files
    if run_params.end is None or run_params.beg is None:
        logger.info("Attempting to find beg and end using on-disk files")
        comp_config = config.load_config_file()
        base_dir = comp_config['base_product_dir']
        beg, end = fpio.find_combined_beg_end(run_params.obsid, base_path=base_dir)
        if beg is None or end is None:
            logger.error("Combined files not on disk. Please manually input beginning and end")
            sys.exit(1)
        else:
            run_params.set_beg(beg)
            run_params.set_end(end)

    work_out_what_to_do(run_params)
