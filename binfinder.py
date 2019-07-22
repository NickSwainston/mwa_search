#!/usr/bin/env python3

import os
import glob
import logging
import argparse
import sys
import data_process_pipeline
from job_submit import submit_slurm
import plotting_toolkit

logger = logging.getLogger(__name__)


#----------------------------------------------------------------------
def bestprof_info(prevbins=None, filename=None):
    #returns a dictionary that includes the relevant information from the .bestprof file
    if filename is not None:
        bestprof_path = filename
    else:
        bestprof_path = glob.glob("*{0}*bestprof".format(prevbins))[0]

    #open the file and read the info into a dictionary
    info_dict = {}
    f = open(bestprof_path, "r")
    lines = f.read()
    lines = lines.split("\n")
    #info:
    info_dict["obsid"] = lines[0].split()[4].split("_")[0]
    info_dict["pulsar"] = lines[1].split()[3].split("_")[1]
    info_dict["nbins"] = lines[9].split()[4]
    info_dict["chi"] = lines[12].split()[4]
    info_dict["sn"] = lines[13].split()[4][2:]
    info_dict["dm"] = lines[14].split()[4]
    info_dict["period"] = lines[15].split()[4] #in ms
    info_dict["period_error"] = lines[15].split()[6]
    f.close()
    return info_dict

#----------------------------------------------------------------------
def submit_to_db(run_params, prof_name):

    logger.info("submitting profile to database: {0}".format(prof_name))

    #Add path to filenames for submit script
    mydict = bestprof_info(filename = prof_name)
    ppps = os.getcwd() + "/" + glob.glob("*{0}*{1}*.pfd.ps".format(mydict["nbins"], run_params.pulsar[1:]))[0]
    prof_name = os.getcwd() + "/" + prof_name

    commands = []
    commands.append('submit_to_database.py -o {0} --cal_id {1} -p {2} --bestprof {3} --ppps {4}'\
    .format(run_params.obsid, run_params.cal_id, run_params.pulsar, prof_name, ppps))
    commands.append('echo "submitted profile to database: {0}"'.format(prof_name))


    if run_params.stop==False:
        #Run stokes fold
        commands.append("data_process_pipeline.py -d {0} -O {1} -p {2} -o {3} -b {4} -L {5}\
                        --mwa_search {6} --vcs_tools {7} -m s"\
                        .format(run_params.pointing_dir, run_params.cal_id, run_params.pulsar,\
                        run_params.obsid, run_params.best_bins, run_params.loglvl, run_params.mwa_search,\
                        run_params.vcs_tools))

    commands.append('echo "Searching for pulsar using the pipeline to test the pipelines effectivness"')
    commands.append('mwa_search_pipeline.py -o {0} -a --search --pulsar {1} -O {2}\
                    --code_comment "Known pulsar auto test"'.format(run_params.obsid, run_params.pulsar,\
                    run_params.cal_id))


    name = "Submit_{0}_{1}".format(run_params.pulsar, run_params.obsid)
    batch_dir = "/group/mwaops/vcs/{0}/batch/".format(run_params.obsid)

    submit_slurm(name, commands,\
                 batch_dir=batch_dir,\
                 slurm_kwargs={"time": "00:05:00"},\
                 module_list=['mwa_search/{0}'.format(run_params.mwa_search)],\
                 submit=True, vcstools_version="{0}".format(run_params.vcs_tools))

#----------------------------------------------------------------------
def check_conditions(threshold, prevbins):

    #returns a dictionary of a bunch of stuff that decides if and how to run prepfold
    info_dict = bestprof_info(prevbins=prevbins)
    condition_dict = {}
    if float(info_dict["sn"]) < threshold:
        condition_dict["sn_good"] = False
        logger.info("Signal to noise ratio of the previous run was below the threshold")
    else:
        condition_dict["sn_good"] = True

    if float(info_dict["chi"]) < 4.0:
        condition_dict["chi_good"] = False
        logger.info("Chi value of the previous run was below 4")
    else:
        condition_dict["chi_good"] = True

    if int(float(info_dict["sn"])) == 0:
        condition_dict["sn_nonzero"] = False
        logger.info("The singal to noise ratio for this file is zero. Using chi for evalutation")
    else:
        condition_dict["sn_nonzero"] = True

    if int(info_dict["nbins"]) > int(float(info_dict["period"]))/1000 * 10000: #the 10k is the 10khz time res of MWA
        condition_dict["sampling_good"] = False
        logger.info("The maximum sampling frequency for this pulsar has been reached")
    else:
        condition_dict["sampling_good"] = True


    return condition_dict

#----------------------------------------------------------------------
def get_best_profile(pointing_dir, threshold):

    #find all of the relevant bestprof profiles in the pointing directory
    bestprof_names = glob.glob("*bins*{0}*.bestprof".format(run_params.pulsar[1:]))
    if len(bestprof_names)==0:
        logger.error("No bestprofs found in directory! Exiting")
        sys.exit(1)

    #throw all of the information from each bestprof into an array
    bin_order = []
    sn_order = []
    chi_order = []
    for prof in bestprof_names:
        prof_info = bestprof_info(filename=prof)
        bin_order.append(int(prof_info["nbins"]))
        sn_order.append(float(prof_info["sn"]))
        chi_order.append(float(prof_info["chi"]))
    bin_order, sn_order, chi_order = zip(*sorted(zip(bin_order, sn_order, chi_order)))
    bin_order = bin_order[::-1]
    sn_order = sn_order[::-1]
    chi_order = chi_order[::-1]

    #now find the one with the most bins that meet the sn and chi conditions
    best_i = None
    prof_name = None
    for i in range(len(bin_order)):
        if float(sn_order[i])>=threshold and float(chi_order[i])>=4.0:
            best_i = i
            break
    if best_i is None:
        logger.info("No profiles fit the threshold parameter")
    else:
        logger.info("Adequate profile found with {0} bins".format(bin_order[best_i]))
        prof_name = glob.glob("*{0}*{1}*.bestprof".format(bin_order[best_i], run_params.pulsar[1:]))[0]

    return prof_name

#----------------------------------------------------------------------
def submit_multifold(run_params, nbins=64):

    job_ids = []
    for i, pointing in enumerate(run_params.pointing_dir):
        logger.info("submitting pointing:{0}".format(pointing))
        os.chdir(pointing)
        #create slurm job:
        commands = []
        #load presto module here because it uses python 2
        commands.append('echo "Folding on known pulsar {0}"'.format(run_params.pulsar))
        commands.append('psrcat -e {0} > {0}.eph'.format(run_params.pulsar))
        commands.append("sed -i '/UNITS           TCB/d' {0}.eph".format(run_params.pulsar))
        commands.append("prepfold -o {0}_{2}_bins -noxwin -nosearch -runavg -noclip -timing {1}.eph\
                        -nsub 256 1*.fits -n {2}".format(run_params.obsid, run_params.pulsar, nbins))
        commands.append('errorcode=$?')
        commands.append('pulsar={}'.format(run_params.pulsar[1:]))
        pulsar_bash_string = '${pulsar}'
        #Some old ephems don't have the correct ra and dec formating and
        #causes an error with -timing but not -psr
        commands.append('if [ "$errorcode" != "0" ]; then')
        commands.append('   echo "Folding using the -psr option"')
        commands.append('   prepfold -o {0}_{2}_bins -noxwin -nosearch -runavg -noclip -psr {1}\
                        -nsub 256 1*.fits -n {2}'.format(run_params.obsid, run_params.pulsar, nbins))
        commands.append('   pulsar={}'.format(run_params.pulsar))
        commands.append('fi')
        commands.append('rm {0}.eph'.format(run_params.pulsar))


        name = "multifold_{0}_{1}".format(run_params.pulsar, i)
        batch_dir = "/group/mwaops/vcs/{0}/batch/".format(run_params.obsid)
        myid = submit_slurm(name, commands,\
                    batch_dir=batch_dir,\
                    slurm_kwargs={"time": "1:00:00"},\
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
        p = p + pointing + " "

    commands=[]
    commands.append("binfinder.py -m b -d {0} -O {1} -p {2} -o {3} -L {4} {5} --vcs_tools {6}\
                    --mwa_search {7} --force_initial -p {8}"\
                    .format(p, run_params.cal_id, run_params.pulsar, run_params.obsid, run_params.loglvl,\
                    stop, run_params.vcs_tools, run_params.mwa_search, run_params.pulsar))

    name="best_fold_{0}".format(run_params.pulsar)
    batch_dir = "/group/mwaops/vcs/{0}/batch/".format(run_params.obsid)
    myid = submit_slurm(name, commands,\
            batch_dir=batch_dir,\
            slurm_kwargs={"time": "00:10:00"},\
            module_list=["mwa_search/k_smith",\
                        "presto/no-python"],\
            submit=True, depend=job_ids,\
            vcstools_version="multi-pixel_beamform")


#----------------------------------------------------------------------
def submit_prepfold(run_params, nbins=32, finish=False):

    if nbins is not int:
        nbins = int(float(nbins))

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
    #load presto module here because it uses python 2
    commands.append('echo "Folding on known pulsar {0}"'.format(run_params.pulsar))
    commands.append('psrcat -e {0} > {0}.eph'.format(run_params.pulsar))
    commands.append("sed -i '/UNITS           TCB/d' {0}.eph".format(run_params.pulsar))
    commands.append("prepfold -o {0}_{2}_bins -noxwin -nosearch -runavg -noclip -timing {1}.eph\
                    -nsub 256 1*.fits -n {2}".format(run_params.obsid, run_params.pulsar, nbins))
    commands.append('errorcode=$?')
    commands.append('pulsar={}'.format(run_params.pulsar[1:]))
    pulsar_bash_string = '${pulsar}'
    #Some old ephems don't have the correct ra and dec formating and
    #causes an error with -timing but not -psr
    commands.append('if [ "$errorcode" != "0" ]; then')
    commands.append('   echo "Folding using the -psr option"')
    commands.append('   prepfold -o {0}_{2}_bins -noxwin -nosearch -runavg -noclip -psr {1}\
                    -nsub 256 1*.fits -n {2}'.format(run_params.obsid, run_params.pulsar, nbins))
    commands.append('   pulsar={}'.format(run_params.pulsar))
    commands.append('fi')
    commands.append('rm {0}.eph'.format(run_params.pulsar))

    if finish==False:
        #Rerun this script
        commands.append('echo "Running script again. Passing prevbins = {0}"'.format(nbins))
        launch_line += " -m f"
    else:
        #Run again only once and without prepfold
        commands.append('echo "Running script again without folding. Passing prevbins = {0}"'.format(nbins))
        launch_line += " -m e"

    commands.append(launch_line)

    name = "binfinder_{0}_{1}".format(run_params.pulsar, nbins)
    batch_dir = "/group/mwaops/vcs/{0}/batch/".format(run_params.obsid)
    submit_slurm(name, commands,\
                batch_dir=batch_dir,\
                slurm_kwargs={"time": "2:00:00"},\
                module_list=['mwa_search/{0}'.format(run_params.mwa_search),\
                            'presto/no-python'],\
                submit=True, vcstools_version="{0}".format(run_params.vcs_tools))
    logger.info("Job successfully submitted")



#----------------------------------------------------------------------
def find_best_pointing(run_params, nbins=64):

    bestprof_info_list = []
    for pointing in run_params.pointing_dir:
        os.chdir(pointing)
        logger.info("searching directory: {0}".format(pointing))
        prof_name = glob.glob("*{0}*{1}*.bestprof".format(nbins, run_params.pulsar[1:]))[0]
        bestprof_info_list.append(bestprof_info(filename=prof_name))

    #now we loop through all the info and find the best one
    best_sn = 0.0
    best_i = -1
    for i, info_dict in enumerate(bestprof_info_list):
        if float(info_dict["chi"])>=4.0 and float(info_dict["sn"])>best_sn:
            best_sn = float(info_dict["sn"])
            best_i = i

    if best_i<0 and best_sn<5.0:
        logger.info("No pulsar found in pointings. Exiting...")
    else:
        logger.info("Pulsar found in pointings. Running binfinder script on pointing: {0}"\
                    .format(run_params.pointing_dir[best_i]))
    
        if run_params.stop==True:
            stop = "-S"
        else:
            stop = ""
        commands = []
        commands.append("binfinder.py -d {0} -t {1} -O {2} -o {3} -L {4} {5} --vcs_tools {6}\
                        --mwa_search {7} -p {8} -m f"\
                        .format(run_params.pointing_dir[best_i], run_params.threshold, run_params.cal_id,\
                        run_params.obsid, run_params.loglvl, stop, run_params.vcs_tools,\
                        run_params.mwa_search, run_params.pulsar))

    name = "binfinder_{0}_{1}".format(run_params.pulsar, nbins)
    batch_dir = "/group/mwaops/vcs/{0}/batch/".format(run_params.obsid)
    submit_slurm(name, commands,\
                batch_dir=batch_dir,\
                slurm_kwargs={"time": "2:00:00"},\
                module_list=['mwa_search/{0}'.format(run_params.mwa_search),\
                            'presto/no-python'],\
                submit=True, vcstools_version="{0}".format(run_params.vcs_tools))
    logger.info("Job successfully submitted")




#----------------------------------------------------------------------
def iterate_bins(run_params):

    #If this is not the first run:
    if run_params.prevbins is not None:
        #Ensuring prevbins is in the correct int format
        run_params.set_prevbins(int(float(run_params.prevbins)))
        #get information of the previous run
        info_dict = bestprof_info(prevbins=run_params.prevbins)

        #Check to see if SN and chi are above threshold
        #If continue == True, prepfold will run again
        cont = False
        condition_dict = check_conditions(run_params.threshold, run_params.prevbins)
        if condition_dict["sn_nonzero"] is False:
            if condition_dict["sn_good"] is False:
                cont = True
        elif condition_dict["chi_good"] is False:
            cont = True


        finish = False
        if cont is True:
            #Choosing the number of bins to use
            nbins = int(float(info_dict["nbins"]))/2
            while nbins>int(float(info_dict["period"])/1000 * 10000):
                logger.info("Time sampling limit reached. Bins will be reduced")
                nbins = nbins/2
                if nbins<=32:
                    break
            if nbins<=32:
                logger.warn("Minimum number of bins hit. Script will run once more")
                finish=True

            #create slurm job:
            submit_prepfold(run_params, nbins=nbins, finish=finish)

        else:
            #Threshold reached, find the best profile and submit it to DB
            logger.info("Signal to noise or Chi above threshold at {0} bins".format(info_dict["nbins"]))
            logger.info("Finding best profile in directory and submitting to database")
            bestprof = get_best_profile(run_params.pointing_dir, run_params.threshold)
            if bestprof==None:
                logger.info("No profiles found with threshold parameters. Attempting threshold=5.0")
                bestprof = get_best_profile(run_params.pointing_dir, 5.0)
                if bestprof==None:
                    logger.info("Non-detection on pulsar {0}".format(run_params.pulsar))
                    logger.info("Exiting....")
                    sys.exit(0)

            run_params.set_best_bins(int(float(info_dict["nbins"])))
            submit_to_db(run_params, bestprof)

    else:
        #This is the first run
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
    other.add_argument("--vcs_tools", type=str, default="multi-pixel_beamform", help="The version of vcs_tools to use. Default: multi-pixel_beamform")

    modeop = parser.add_argument_group("Mode Options:")
    modeop = required.add_argument("-m", "--mode", type=str, help="""The mode in which to run binfinder\n\
                        'f' - Finds an adequate number of bins to fold on\n\
                        'e' - Folds once on the default number of bins and submits the result to\
                        the database. NOT RECOMMENDED FOR MANUAL INPUT\n\
                        'm' - Use this mode if this is part of a multi-beam observation. This will\
                        find the best detection, if any, out of many pointings\n\
                        'b' - Finds the best detection out of a set of pointing directories""")


    non_user = parser.add_argument_group("Non-User input Options:")
    non_user.add_argument("--prevbins", type=int, default=None, help="The number of bins used in prepfold on the previous run. Not necessary for initial runs")

    args = parser.parse_args()

    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    ch.setLevel(loglevels[args.loglvl])
    formatter = logging.Formatter('%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

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

    """
    NOTE: for some reason, you need to run prepfold from the directory it outputs to if you want it to properly make an image. The script will make this work regardless by using os.chdir
    """
    if run_params.mode is not "m" and run_params.mode is not "b":
        os.chdir(run_params.pointing_dir)


    if run_params.mode=="e":
        logger.info("Submitting to database")
        prof_name = get_best_profile(run_params.pointing_dir, run_params.threshold)

        if prof_name==None:
            logger.info("No profile found for input threshold. Trying again with Threshold=5.0")
            prof_name = get_best_profile(run_params.pointing_dir, 5.0)
            run_params.stop_now()

        if prof_name==None:
            logger.info("Non detection - no adequate profiles. Exiting....")
            sys.exit(0)

        #plot the profile properly
        plotting_toolkit.plot_bestprof("{0}/{1}".format(run_params.pointing_dir,run_params.prof_name),\
                                        out_dir=run_params.pointing_dir)
        mydict = bestprof_info(filename=prof_name)
        run_params.set_best_bins(int(float(mydict["nbins"])))
        submit_to_db(run_params, prof_name)

    elif run_params.mode=="m":
        submit_multifold(run_params)
    elif run_params.mode=="b":
        find_best_pointing(run_params, nbins=64)
    elif run_params.mode=="f":
        iterate_bins(run_params)
    else:
        logger.error("Unreognized mode. Please run again with a proper mode selected.")
