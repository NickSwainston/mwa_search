#!/usr/bin/env python3

import os
import glob
import logging
import argparse
import sys
import logging
from job_submit import submit_slurm

logger = logging.getLogger(__name__)

def bestprof_info(prevbins=None, filename=None):
    #returns a dictionary that includes the relevant information from the .bestprof file
    if filename is not None
        bestprof_path = filename
    #elif prevbins == None:
    #    bestprof_path = glob.glob("*PSR**bestprof")[0]
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

def find_obsid_pulsar(pointing_dir):
    
    #given a pointing directory, returns a dictionary containing the obsid and pulsar name
    mydict = {}
    mydict["obsid"] = pointing_dir.split("/")[4]
    pulsar_coords = pointing_dir.split("/")[6].split("_")
    pulsar_1 = pulsar_coords[0].split(":")[0] + pulsar_coords[0].split(":")[1]
    pulsar_2 = pulsar_coords[1].split(":")[0] + pulsar_coords[1].split(":")[1]
    mydict["pulsar"] = pulsar_1 + pulsar_2
    
    return mydict

def submit_to_db(pointing_dir, prof_name, cal_id)


    mydict = bestprof_info(prof_name)
    obsid = mydict["obsid"]
    pulsar = mydict["pulsar"]
    ppps = glob.glob("*{0}*.pfd.ps".format(mydict[nbins]))

    commands = []
    commands.append('submit_to_database.py -o {0} --cal_id {1} -p {2} --bestprof {3} --ppps {4}'.format(obsid, cal_id, pulsar, prof_name, ppps))

    commands.append('   echo "Searching for pulsar using the pipeline to test the pipelines effectivness"')
    commands.append('   mwa_search_pipeline.py -o {0} -b {1} -e {2} --search --pulsar {3} -O {4} --code_comment "Known pulsar auto test {3}"'.format(obsid, begin, end, pulsar, cal_id))

    logger.info("Profile submitted to database: {0}".format(prof_name))



def find_best_profile(pointing_dir, threshold)

    #find all of the bestprof profiles in the pointing directory
    bestprof_names = glob.glob("*.bestprof")
    if len(bestprofs)==0:
        logger.error("No bestprofs found in directory! Exiting")
        sys.exit(1)

    #throw all of the information from each bestprof into an array   
    bin_order = []
    sn_order = []
    chi_order = []
    for prof in bestprof_names: 
        prof_info = bestprof_info(prof)
        bin_order.append(prof_info["nbins"])
        sn_order.append(prof_info["sn"])
        chi_order.append(prof_info["chi"])
    bin_order, sn_order, chi_order = zip(*sorted(zip(bin_order, sn_order, chi_order)))
    bin_order = [::-1]
    sn_order = [::-1]
    chi_order = [::-1]

    #now find the one with the most bins that meet the sn and chi conditions
    best_i = None    
    for i,bins in enumerate(bin_order):
        if sn_order[i]>=threshold and chi[i]>=4.0:
            best_i = i
            break
    if best_i is None:
        logger.info("No profiles fit the threshold parameter. Exiitng...")
        sys.exit(0)
    
    profile_name = glob.glob("*{0}*.bestprof".format(bin_order[i]))

    
def submit_prepfold(pulsar, obsid, threshold, nbins=32, rerun_pf=True, loglvl="INFO")

    logger.info("Submitting job for {0} bins".format(nbins))
    #create slurm job:
    commands = []
    #load presto module here because it uses python 2
    commands.append('echo "Folding on known pulsar"'.format(pulsar))
    commands.append('psrcat -e {0} > {0}.eph'.format(pulsar))
    commands.append("sed -i '/UNITS           TCB/d' {0}.eph".format(pulsar))
    commands.append("prepfold -o {0}_{2}_bins -noxwin -nosearch -runavg -noclip -timing {1}.eph -nsub 256 1*fits -n {2}".format(obsid, pulsar, nbins))
    commands.append('errorcode=$?')
    commands.append('pulsar={}'.format(pulsar[1:]))
    pulsar_bash_string = '${pulsar}'
    #Some old ephems don't have the correct ra and dec formating and
    #causes an error with -timing but not -psr
    commands.append('if [ "$errorcode" != "0" ]; then')
    commands.append('   echo "Folding using the -psr option"')
    commands.append('   prepfold -o {0}_{2}_bins -noxwin -nosearch -runavg -noclip -psr {1} -nsub 256 1*fits -n {2}'.format(obsid, pulsar, nbins))
    commands.append('   pulsar={}'.format(pulsar))
    commands.append('fi')
    commands.append('rm {0}.eph'.format(pulsar))

    if rerun_pf == True:
        #Rerun this script
        commands.append('echo "Running script again. Passing prevbins = {0}"'.format(nbins))
        commands.append('binfinder.py -p {0} -t {1} -L {2} --prevbins {3}'.format(pointing_dir, threshold, loglvl, nbins))
    else:
        #submit either this run or the previous run to slurm by rerunning this script with "--submit"
        commands.append()   

 
    name = "binfinder_{0}_{1}".format(pulsar, nbins)
    batch_dir = "/group/mwaops/vcs/{0}/batch".format(obsid)
    submit_slurm(name, commands,
                batch_dir=batch_dir,
                slurm_kwargs={"time": "5:00:00"},
                module_list=['mwa_search/k_smith',
                            'presto/no-python'],
                submit=True, vcstools_version="multi-pixel_beamform")
    logger.info("Job successfully submitted")


def find_optimal_prof(pointing_dir, prevbins, threshold, loglvl):

    #Ensuring prevbins is in the correct int format
    if prevbins is not None:
        prevbins = int(float(prevbins))
    
    #get information of the previous run
    info_dict = bestprof_info(prevbins=prevbins)

    logger.info("previous bins:    {0}".format(info_dict["nbins"]))
    logger.info("previous signal/noise:     {0}".format(info_dict["sn"]))
    logger.info("previous chi:              {0}".format(info_dict["chi"]))
    
    #Check to see if the threshold has been reached
    #If continue == True, prepfold will run again
    cont = True
    condition_dict = check_conditions(threshold, prevbins)
    if condition_dict["sn_nonzero"]==False:
        if condition_dict["sn_good"]==False:
            cont = False
    if condition_dict["chi_good"]==False:
        cont = False    

    rerun_pf == True 
    if cont==True:
        #Choosing the number of bins to use
        if condition_dict["sampling_good"] and condition_dict["nbins_good"]:
            nbins = int(float(info_dict["nbins"]))*2
        else:
            nbins = int(float(info_dict["period"])/1000 * 10000)
            rerun_pf = False
            logger.info("Prepfold will run once more and then terminate")

        pulsar = info_dict["pulsar"]
        obsid = info_dict["obsid"]
       
        #create slurm job:
        submit_prepfold(pulsar, obsid, nbins, terminate=False, loglvl)


    else:
        logger.info("Signal to noise or Chi limit breached at {0} bins".format(info_dict["nbins"]))
        logger.info("Exiting....")
        submit_bestprof()

    if rerun_pf == False:
        #Job done. Exit
        logger.info("Maximum number of bins found: {0}".format(nbins))
        logger.info("Exiting....") 

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

    if int(info_dict["nbins"])*2 > 16384:
        condition_dict["nbins_good"] = False
        logger.info("The number of bins is above the maximum allowed of 16384")
    else: 
        condition_dict["nbins_good"] = True

    if int(float(info_dict["sn"])) == 0:
        condition_dict["sn_nonzero"] = False
        logger.info("The singal to noise ratio for this file is zero. Using chi for evalutation")
    else:
        condition_dict["sn_nonzero"] = True
    
    if int(info_dict["nbins"])*2 > int(float(info_dict["period"]))/1000 * 10000: #the 10k is the 10khz time res of MWA
        condition_dict["sampling_good"] = False
        logger.info("The maximum sampling frequency for this pulsar has been reached")
    else:
        condition_dict["sampling_good"] = True


    return condition_dict


if __name__ == '__main__':
    #dictionary for choosing log-levels
    loglevels = dict(DEBUG=logging.DEBUG,
                    INFO=logging.INFO,
                    WARNING=logging.WARNING,
                    ERROR = logging.ERROR)

    #Arguments
    parser = argparse.ArgumentParser(description="""Finds the optimal number of bins for a given pointing""")

    parser.add_argument_group("User Input Options:")
    parser.add_argument("-d", "--pointing_dir", type=str, help="Pointing directory that contains the spliced fits files")
    parser.add_argument("-t", "--threshold", type=float, default=10.0, help="The signal to noise threshold to stop at. Default = 10.0")
    parser.add_argument("-L", "--loglvl", type=str, default="INFO", help="Logger verbosity level. Default: INFO", choices=loglevels.keys())

    parser.add_argument_group("Batch Script Options:")
    parser.add_argument("--prevbins", type=str, default=None, help="The number of bins used in prepfold on the previous run. Not necessary for initial runs")
    parser.add_argument("--submit", action="store_false", help="Use this tag if you are ready to submit a .bestprof to the database")

    args = parser.parse_args()

    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    ch.setLevel(loglevels[args.loglvl])
    formatter = logging.Formatter('%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    if args.pointing_dir == None:
        logger.error("No pointing directory supplied. Please specify the pointing directory path and rerun")
        sys.exit(1)

    """
    NOTE: for some reason, you need to run prepfold from the directory it outputs to if you want it to properly make an image. The script will make this work regardless by using os.chdir
    """
    os.chdir(args.pointing_dir)
    
    if args.submit==True:
        submit_db(args.pointing_dir)
        sys.exit(0)

    if len(glob.glob(*.bestprof))==0:
        initial == True

    if initial == False:
        find_optimal_prof(args.pointing_dir, args.prevbins, args.threshold, args.loglvl)
    else:
        o_p = find_obsid_pulsar(pointing_dir)
        submit_prepfold(o_p["pulsar"], o_p["obsid"],  args.threshold, loglvl=args.loglvl)
