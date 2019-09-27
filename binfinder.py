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
import file_maxmin
import psrqpy
import pandas as pd
logger = logging.getLogger(__name__)

#get ATNF db location
try:
    ATNF_LOC = os.environ['PSRCAT_FILE']
except:
    logger.warn("ATNF database could not be loaded on disk. This may lead to a connection failure")
    ATNF_LOC = None


#----------------------------------------------------------------------
def add_prepfold_to_commands(commands, pointing, pulsar, obsid, use_mask=True, start=None, end=None, nbins=100, ntimechunk=120, dmstep=1, period_search_n=1):

    #find the beginning and end of the pulsar's beam coverage for this obs
    if start==None or end==None:
        start, end = pulsar_beam_coverage(obsid, pulsar)
        logger.info("start and end of pulsar beam coverage for on-disk files:{0}, {1}".format(start, end))
        if start>=1. or end<0.:
            logger.error("pulsar is not in beam for any of the on-disk files. Ending...")
            sys.exit(1)

    comp_config = config.load_config_file()
    #Figure out whether or not to input a mask
    if use_mask == True:
        check_mask = glob.glob("{0}{1}/incoh/*.mask".format(comp_config['base_product_dir'], obsid))
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
    variables += "-start {0} -end {1} ".format(start, end) 
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
    pulsar_bash_string = '${pulsar}'

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
def pulsar_beam_coverage(obsid, pulsar, beg=None, end=None):
    #returns the beginning and end time as a fraction that a pulsar is in the primary beam for the obsid files
    #beg and end should only be supplied if the files are not present on the system

    #find the enter and exit times of pulsar normalized with the observing time
    names_ra_dec = fpio.grab_source_alog(pulsar_list=[pulsar])
    beam_source_data, _ = fpio.find_sources_in_obs([obsid], names_ra_dec)

    enter_obs_norm = beam_source_data[obsid][0][1]
    exit_obs_norm = beam_source_data[obsid][0][2]

    if beg is None and end is None:
        #find the beginning and end time of the observation FILES you have on disk
        files_beg, files_end = checks.find_beg_end(obsid)
        files_duration = files_end - files_beg + 1
    else:
        #uses manually input beginning and end times to find beam coverage
        files_beg = beg
        files_end = end
        files_duration = files_end - files_beg + 1

    #find how long the total observation is (because that's what enter and exit uses)
    obs_beg, obs_end, obs_dur = file_maxmin.print_minmax(obsid)
    obs_dur = obs_end-obs_beg

    #times the source enters and exits beam
    time_enter = obs_beg + obs_dur*enter_obs_norm
    time_exit = obs_beg + obs_dur*exit_obs_norm

    #normalised time the source enters/exits the beam in the files
    enter_files = (time_enter-files_beg)/files_duration
    exit_files = (time_exit-files_beg)/files_duration

    if enter_files<0.:
        enter_files=0.
    if exit_files>1.:
        exit_files=1.
    if enter_files>1.:
        logger.warn("source {0} is not in the beam for the files on disk".format(pulsar))
    if exit_files<0.:
        logger.warn("source {0} is not in the beam for the files on the disk".format(pulsar))

    return enter_files, exit_files

#----------------------------------------------------------------------
def bestprof_info(prevbins=None, filename=None):
    #returns a dictionary that includes the relevant information from the .bestprof file
    if filename is not None:
        bestprof_path = filename
    else:
        bestprof_path = glob.glob("*{0}_bins*bestprof".format(prevbins))[0]
    #open the file and read the info into a dictionary
    info_dict = {}
    f = open(bestprof_path, "r")
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
    #returns the minimum number of bins you can use for this pulsar based on MWA sampling rate
    
    query = psrqpy.QueryATNF(params=["P0"], psrs=[pulsar], loadfromdb=ATNF_LOC).pandas
    period = query["P0"][0]
    min_bins = int(period/sampling_rate + 1) #the +1 is to round the limit up every time
    logger.debug("Bin limit: {0}".format(min_bins))
    return min_bins
    
#----------------------------------------------------------------------
def submit_to_db(run_params, ideal_bins):

    logger.info("submitting profile to database: {0}".format(run_params.bestprof))
    #Add path to filenames for submit script
    cwd = os.getcwd()
    ppps = cwd + "/" + glob.glob("*{0}_bins*{1}*.pfd.ps".format(ideal_bins, run_params.pulsar[1:]))[0]
    bestprof_name = cwd + "/" + glob.glob("*{0}_bins*{1}*.pfd.bestprof".format(ideal_bins, run_params.pulsar[1:]))[0]
    png_output = cwd +  "/" + glob.glob("*{0}_bins*{1}*.png".format(ideal_bins, run_params.pulsar[1:]))[0]
    pfd = cwd + "/" + glob.glob("*{0}_bins*{1}*.pfd".format(ideal_bins, run_params.pulsar[1:]))[0]
    
    #do the same for 100 bin profiles
    ppps_100 = cwd + "/" + glob.glob("*_100_bins*{0}*.pfd.ps".format(run_params.pulsar[1:]))[0]
    bestprof_name_100 = cwd + "/" + glob.glob("*_100_bins*{0}*.pfd.bestprof".format(run_params.pulsar[1:]))[0]
    png_output_100 = cwd +  "/" + glob.glob("*_100_bins*{0}*.png".format(run_params.pulsar[1:]))[0]
    pfd_100 = cwd + "/" + glob.glob("*_100_bins*{0}*.pfd".format(run_params.pulsar[1:]))[0]
    
    products = [ppps, bestprof_name, png_output, pfd,\
            ppps_100, bestprof_name_100, png_output_100, pfd_100]
    
    #move all of these data products to a suitable directory
    data_dir = "/group/mwaops/vcs/{0}/data_products/{1}".format(run_params.obsid, run_params.pulsar)
    for product in products:
        data_process_pipeline.copy_data(product, data_dir)

    commands = []
    commands.append('submit_to_database.py -o {0} --cal_id {1} -p {2} --bestprof {3} --ppps {4}'\
    .format(run_params.obsid, run_params.cal_id, run_params.pulsar, bestprof_name, ppps))
    commands.append('echo "submitted profile to database: {0}"'.format(bestprof_name))


    if run_params.stop==False:
        #Run stokes fold
        commands.append("data_process_pipeline.py -d {0} -O {1} -p {2} -o {3} -b {4} -L {5}\
                        --mwa_search {6} --vcs_tools {7} -m s"\
                        .format(run_params.pointing_dir, run_params.cal_id, run_params.pulsar,\
                        run_params.obsid, run_params.best_bins, run_params.loglvl, run_params.mwa_search,\
                        run_params.vcs_tools))

    #commands.append('echo "Searching for pulsar using the pipeline to test the pipelines effectivness"')
    #commands.append('mwa_search_pipeline.py -o {0} -a --search --pulsar {1} -O {2}\
    #                --code_comment "Known pulsar auto test"'.format(run_params.obsid, run_params.pulsar,\
    #                run_params.cal_id))


    name = "Submit_{0}_{1}".format(run_params.pulsar, run_params.obsid)
    comp_config = config.load_config_file()
    batch_dir = "{0}{1}/batch/".format(comp_config['base_product_dir'], run_params.obsid)

    submit_slurm(name, commands,\
                 batch_dir=batch_dir,\
                 slurm_kwargs={"time": "00:05:00"},\
                 module_list=['mwa_search/{0}'.format(run_params.mwa_search)],\
                 submit=True, vcstools_version="{0}".format(run_params.vcs_tools))

#----------------------------------------------------------------------
def get_best_profile(pointing_dir, pulsar, threshold=10):

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
            if sn_chi_test(sn_order[i], chi_order[i]) == True:
                best_i = i
                break

    if best_i is None:
        logger.info("No profiles fit the threshold parameter")
        return None, None
    else:
        logger.info("Adequate profile found with {0} bins".format(bin_order[best_i]))
        prof_name = glob.glob("*{0}_bins*{1}*.bestprof".format(bin_order[best_i], pulsar[1:]))[0]
        return prof_name, bin_order[best_i]

#----------------------------------------------------------------------
def sn_chi_test(sn, chi, sn_thresh=10., chi_thresh=4.):
    test = False
    if sn >= sn_thresh and chi >= chi_thresh:
        test = True
    elif sn == 0. and chi >= chi_thresh:
        test = True
    return test

#----------------------------------------------------------------------
def submit_multifold(run_params, nbins=100):

    job_ids = []
    comp_config=config.load_config_file()

    #Check beam coverage for the pulsar
    start, end = pulsar_beam_coverage(run_params.obsid, run_params.pulsar)
    logger.info("start and end of pulsar beam coverage for on-disk files:{0}, {1}".format(start, end))
    if start>=1. or end<0.:
        logger.error("pulsar is not in beam for any of the on-disk files. Ending...")
        sys.exit(1)

    for i, pointing in enumerate(run_params.pointing_dir):
        logger.info("submitting pointing:{0}".format(pointing))
        #os.chdir(pointing)
        #create slurm job:
        commands = []
        commands = add_prepfold_to_commands(commands, pointing, run_params.pulsar, run_params.obsid,\
                    start=start, end=end, nbins=nbins)

        name = "multifold_binfind_{0}_{1}".format(run_params.pulsar, i)
        batch_dir = "{0}{1}/batch/".format(comp_config['base_product_dir'], run_params.obsid)
        myid = submit_slurm(name, commands,\
                    batch_dir=batch_dir,\
                    slurm_kwargs={"time": "2:00:00"},\
                    module_list=['mwa_search/{0}'.format(run_params.mwa_search),\
                                'presto/no-python'],\
                    submit=True, vcstools_version="{0}".format(run_params.vcs_tools))


        job_ids.append(myid)

    #Now submit the check script
    if run_params.stop==True:
        stop="-S"
    else:
        stop=""

    p = ""
    for pointing in run_params.pointing_dir:
        p += " " + pointing

    commands=[]
    commands.append("binfinder.py -m b -d {0} -O {1} -p {2} -o {3} -L {4} {5} --vcs_tools {6}\
                    --mwa_search {7} --force_initial -p {8}"\
                    .format(p, run_params.cal_id, run_params.pulsar, run_params.obsid, run_params.loglvl,\
                    stop, run_params.vcs_tools, run_params.mwa_search, run_params.pulsar))

    name="best_pointing_{0}".format(run_params.pulsar)
    batch_dir = "{0}{1}/batch/".format(comp_config['base_product_dir'], run_params.obsid)
    myid = submit_slurm(name, commands,\
            batch_dir=batch_dir,\
            slurm_kwargs={"time": "00:30:00"},\
            module_list=['mwa_search/{0}'.format(run_params.mwa_search),\
                        "presto/no-python"],\
            submit=True, depend=job_ids, depend_type="afterany",\
            vcstools_version="master")


#----------------------------------------------------------------------
def submit_prepfold(run_params, nbins=100, finish=False):

    if nbins is not int:
        nbins = int(float(nbins))
    
    comp_config = config.load_config_file()

    #Check to see if there is a 100 bin fold already
    bin_lim = bin_sampling_limit(run_params.pulsar)
    if len(glob.glob("*_100_bins**{0}*bestprof".format(run_params.pulsar)))>0 and bin_lim>100 and nbins is not 100:
        #add a prepfold command for 100 bins
        logger.info("Folding on 100 bins for pointing {0}".format(run_params.pointing))
        commands = []
        commands = add_to_prepfold(commands, run_params.pointing_dir, run_params.pulsar, run_params.obsid, nbins=100)
 
        name = "binfinder_prepfold_only_{0}_100".format(run_params.pulsar)
        batch_dir = "{0}{1}/batch/".format(comp_config['base_product_dir'], run_params.obsid)
        submit_slurm(name, commands,\
                    batch_dir=batch_dir,\
                    slurm_kwargs={"time": "2:00:00"},\
                    module_list=['mwa_search/{0}'.format(run_params.mwa_search),\
                                'presto/no-python'],\
                    submit=True, vcstools_version="{0}".format(run_params.vcs_tools))
        logger.info("Prepfold job successfully submitted: {0}".format(name))

    
    launch_line = "binfinder.py -d {0} -t {1} -O {2} -o {3} -L {4} --prevbins {5} --vcs_tools {6}\
                    --mwa_search {7} -p {8}"\
                    .format(run_params.pointing_dir, run_params.threshold, run_params.cal_id,\
                    run_params.obsid, run_params.loglvl, nbins, run_params.vcs_tools,\
                    run_params.mwa_search, run_params.pulsar)

    if run_params.stop==True:
        launch_line += " -S"

    logger.info("Submitting job for {0} bins".format(nbins))
    #create slurm job:
    commands = []
    if firstrun==True:
        commands = add_to_prepfold(commands, run_params.pointing_dir, run_params.pulsar, run_params.obsid, nbins=100)
    commands = add_prepfold_to_commands(commands, run_params.pointing_dir, run_params.pulsar, run_params.obsid, nbins=nbins)

    if finish==False:
        #Rerun this script
        commands.append('echo "Running script again. Passing prevbins = {0}"'.format(nbins))
        launch_line += " -m f"
    else:
        #Run again only once and without prepfold
        commands.append('echo "Running script again without folding. Passing prevbins = {0}"'.format(nbins))
        launch_line += " -m e"

    commands.append(launch_line)

    comp_config = config.load_config_file()
    name = "binfinder_{0}_{1}".format(run_params.pulsar, nbins)
    batch_dir = "{0}{1}/batch/".format(comp_config['base_product_dir'], run_params.obsid)
    submit_slurm(name, commands,\
                batch_dir=batch_dir,\
                slurm_kwargs={"time": "2:00:00"},\
                module_list=['mwa_search/{0}'.format(run_params.mwa_search),\
                            'presto/no-python'],\
                submit=True, vcstools_version="{0}".format(run_params.vcs_tools))
    logger.info("Job successfully submitted: {0}".format(name))



#----------------------------------------------------------------------
def find_best_pointing(run_params, nbins=100, submit_next=True):
    
    """
    Finds the pointing directory with the highest S/N out of those given in run_params.pointing_dir
    The nbins value is the bin value to search over. Default is 100 
    submit_next tells the function whether or not to continue with the binfinder pipeline. Default=True
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

    if submit_next is not True:
        return run_params.pointing_dir[best_i]
    else:
        run_params.set_pointing_dir(run_params.pointing_dir[best_i])
        #check sampling rate and submit prepfold
        bin_limit = bin_sampling_limit(run_params.pulsar)
        if bin_limit >= 1024: 
            submit_prepfold(run_params, nbins=1024)
        elif bin_limit >=100:
            submit_prepfold(run_params, nbins=bin_limit)
        else:
            logger.info("This pulsar has a low period. Folding with a smaller number of bins")
            submit_prepfold(run_params, nbins=bin_limit, finish=True)
        

#----------------------------------------------------------------------
def iterate_bins(run_params):

    #If this is not the first run:
    if run_params.prevbins is not None:
        #Ensuring prevbins is in the correct int format
        run_params.set_prevbins(int(float(run_params.prevbins)))
        #get information of the previous run
        info_dict = bestprof_info(prevbins=run_params.prevbins)

        #Check to see if SN and chi are above threshold
        #If finish == True, prepfold will run again
        info_dict = bestprof_info(prevbins=run_params.prevbins)         
        sn = info_dict["sn"]
        chi = info_dict["chi"]
        #decide whether to continue based on previous run's SN and chi
        finish = sn_chi_test(sn, chi)
            
        if finish is False:
            #Choosing the number of bins to use
            nbins = info_dict["nbins"]/2
            
            #Comparing nbins to sampling rate
            bin_upper_lim = bin_sampling_limit(run_params.pulsar)

            if nbins < bin_upper_lim and bin_upper_lim<1024:
                logger.info("Time sampling limit reached. Bins will be reduced once more")
                nbins = bin_upper_lim
                finish=True
                
            if nbins<=100:
                logger.warn("Minimum number of bins hit. Script will run once more")
                finish=True

            #create slurm job:
            submit_prepfold(run_params, nbins=nbins, finish=finish)

        else:
            #Threshold reached, find the best profile and submit it to DB
            logger.info("S/N iand Chi above threshold at {0} bins in pointing directory: {1}".format(info_dict["nbins"], run_params.pointing_dir))
            logger.info("Finding best profile in directory and submitting to database")
            bestprof, _ = get_best_profile(run_params.pointing_dir, run_params.threshold)

            if bestprof==None:
                logger.info("Non-detection on pulsar {0}".format(run_params.pulsar))
                logger.info("Exiting....")
                sys.exit(0)

            run_params.set_best_bins(info_dict["nbins"])
            #Plot the bestprof nicely
            prof_path = run_params.pointing_dir + "/" + bestprof
            nice_prof_path = plotting_toolkit.plot_bestprof(prof_path, run_params.pointing_dir, nocrop=run_params.nocrop)
            #copy to data products directory
            data_process_pipeline.copy_data(nice_prof_path, "/group/mwaops/vcs/{0}/data_products/{1}"\
                                            .format(run_params.obsid, run_params.pulsar))
            #submit
            submit_to_db(run_params, bestprof)


    else:
        #This is the first run
        bin_lim = bin_sampling_limit(run_params.pulsar)
        if bin_lim < 1024:
            submit_prepfold(run_params, nbins=bin_lim)
        else:
            submit_prepfold(run_params, nbins=1024)


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
    required.add_argument("-d", "--pointing_dir", nargs="+", help="Pointing directory(s) that contains\
                            the spliced fits files. If mode='m', more than one argument may be supplied.")
    required.add_argument("-O", "--cal_id", type=str, help="The Obs ID of the calibrator")
    required.add_argument("-p", "--pulsar", type=str, default=None, help="The name of the pulsar. eg. J2241-5236")


    other = parser.add_argument_group("Other Options:")
    other.add_argument("-t", "--threshold", type=float, default=10.0, help="The signal to noise threshold to stop at. Default = 10.0")
    other.add_argument("-o", "--obsid", type=str, default=None, help="The observation ID")
    other.add_argument("-L", "--loglvl", type=str, default="INFO", help="Logger verbosity level. Default: INFO", choices=loglevels.keys())
    other.add_argument("--force_initial", action="store_true", help="Use this tag to force the script to treat this as the first run.")
    other.add_argument("-S", "--stop", action="store_true", help="Use this tag to tell binfinder to launch the next step in the data processing pipleline when finished")
    other.add_argument("--mwa_search", type=str, default="master", help="The version of mwa_search to use. Default: master")
    other.add_argument("--vcs_tools", type=str, default="master", help="The version of vcs_tools to use. Default: master")

    modeop = parser.add_argument_group("Mode Options:")
    modeop = required.add_argument("-m", "--mode", type=str, help="""The mode in which to run binfinder\n\
                        User Options:\n\
                        'f' - Finds an adequate number of bins to fold on and submits to database if\n\
                        a detection is found.\n\
                        'm' - Use this mode if this is part of a multi-beam observation. This will\
                        fold on many input pointings and then find the best detection, if any, and\
                        find an adequate number of bins for it.\n\
                        
                        Not for typical user inputs:\n\
                        'e' - Folds once on the default number of bins and submits the result to\
                        the database.\n\
                        'b' - Finds the best detection out of a set of pointing directories.""")


    non_user = parser.add_argument_group("Non-User input Options:")
    non_user.add_argument("--prevbins", type=int, default=None, help="The number of bins used in prepfold on the previous run. Not necessary for initial runs")

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
    elif args.mode == None:
        logger.error("Mode not supplied. Please input a mode from the list of modes and rerun")


    run_params = data_process_pipeline.run_params_class\
                    (args.pointing_dir, args.cal_id,\
                    prevbins=args.prevbins, pulsar=args.pulsar,\
                    obsid=args.obsid, threshold=args.threshold,\
                    stop=args.stop, force_initial=args.force_initial,\
                    mode=args.mode, loglvl=args.loglvl,\
                    mwa_search=args.mwa_search, vcs_tools=args.vcs_tools)


    #NOTE: for some reason, you need to run prepfold from the directory it outputs to if you want it to properly make an image. The script will make this work regardless by using os.chdir

    if run_params.mode is not "m" and run_params.mode is not "b":
        os.chdir(run_params.pointing_dir)


    if run_params.mode=="e":
        logger.info("Submitting to database")
        prof_name, ideal_bins = get_best_profile(run_params.pointing_dir, run_params.pulsar, run_params.threshold)

        if prof_name==None:
            logger.info("Non detection - no adequate profiles. Exiting....")
            sys.exit(0)
        #plot the profile properly
        plotting_toolkit.plot_bestprof("{0}/{1}".format(run_params.pointing_dir, prof_name),\
                                        out_dir=run_params.pointing_dir)
        mydict = bestprof_info(filename=prof_name)
        run_params.set_best_bins(int(float(mydict["nbins"])))
        #Plot the bestprof nicely
        plotting_toolkit.plot_bestprof(prof_name, out_dir=run_params.pointing_dir)
        #submit
        submit_to_db(run_params, ideal_bins)

    elif run_params.mode=="m":
        submit_multifold(run_params)
    elif run_params.mode=="b":
        find_best_pointing(run_params, nbins=100)
    elif run_params.mode=="f":
        iterate_bins(run_params)
    else:
        logger.error("Unreognized mode. Please run again with a proper mode selected.")
