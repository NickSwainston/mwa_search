#!/usr/bin/env python3

import numpy as np
import os
import glob
import logging
import argparse
import sys
import logging
from job_submit import submit_slurm

logger = logging.getLogger(__name__)

def bestprof_info(pointing_dir, prevbins):
    #returns a dictionary that includes the relevant information from the .bestprof file
    if prevbins == None:
        bestprof_path = glob.glob(pointing_dir + "/*PSR**bestprof")[0]
    else:
        bestprof_path = glob.glob(pointing_dir + "/*{0}**.pfd.bestprof".format(prevbins))[0]

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
    info_dict["period"] = lines[15].split()[4]
    info_dict["period_error"] = lines[15].split()[6]
    f.close()
    return info_dict

def find_optimal_prof(pointing_dir, prevbins, threshold, loglvl):

    info_dict = bestprof_info(pointing_dir, prevbins)

    logger.debug("signal/noise: {0}".format(info_dict["sn"]))
    logger.debug("chi:          {0}".format(info_dict["chi"]))
    #check to see if the previous run hit the threshold
    if float(info_dict["sn"]) > threshold and float(info_dict["chi"]) > 4.0:
        #double the number of bins
        nbins = int(info_dict["nbins"])*2
        pulsar = info_dict["pulsar"]
        obsid = info_dict["obsid"]
        
        #create slurm job:
        commands = []
        #load presto module here because it uses python 2
        commands.append('echo "Folding on known pulsar"'.format(pulsar))
        commands.append('psrcat -e {0} > {0}.eph'.format(pulsar))
        commands.append("sed -i '/UNITS           TCB/d' {0}.eph".format(pulsar))
        commands.append("prepfold -o {2}/{0}_{3}_bins.pfd.bestprof -noxwin -runavg -noclip -timing {1}.eph -nsub 256 {2}/1*fits -n {3}".format(obsid, pulsar, pointing_dir, nbins))
        commands.append('errorcode=$?')
        commands.append('pulsar={}'.format(pulsar[1:]))
        pulsar_bash_string = '${pulsar}'
        #Some old ephems don't have the correct ra and dec formating and
        #causes an error with -timing but not -psr
        commands.append('if [ "$errorcode" != "0" ]; then')
        commands.append('   echo "Folding using the -psr option"')
        commands.append('   prepfold -o {2}/{0}_{3}_bins.pfd.bestprof -noxwin -runavg -noclip -psr {1} -nsub 256 {2}/1*fits -n {3}'.format(obsid, pulsar, pointing_dir, nbins))
        commands.append('   pulsar={}'.format(pulsar))
        commands.append('fi')
        commands.append('rm {0}.eph'.format(pulsar))

        #Rerun this script
        commands.append('echo "Running script again using {0} bins"'.format(nbins))
        commands.append('binfinder.py -p {0} -t {1} -L {2} --prevbins {3}'.format(pointing_dir, nbins, threshold, loglvl))

        name = "binfinder_{0}_{1}".format(pulsar, nbins)
        batch_dir = "/group/mwaops/vcs/{0}/batch".format(obsid)
        submit_slurm(name, commands,
                     batch_dir=batch_dir,
                     slurm_kwargs={"time": "5:00:00",
                                   "nice":"90"},
                     module_list=['mwa_search/k_smith',
                                  'presto/no-python'],
                     submit=True, vcstools_version="multi-pixel_beamform")



    else:
        logger.info("Maximum number of bins reached for this threshold")
        logger.info("The last successful script folded with {0} bins".format(int(prevbins)/2))
        logger.info("The unsuccessful fold attempt will be deleted and this script will now end.")
        #delete the unsuccessful bestprof file
        os.remove(glog.glob("*{0}**.pfd.bestprof".format(prevbins))[0])
        #done


if __name__ == '__main__':
    #dictionary for choosing log-levels
    loglevels = dict(DEBUG=logging.DEBUG,
                    INFO=logging.INFO,
                    WARNING=logging.WARNING,
                    ERROR = logging.ERROR)

    #Arguments
    parser = argparse.ArgumentParser(description="""Finds the optimal number of bins for a given pointing""")

    parser.add_argument_group("User Input Options:")
    parser.add_argument("-p", "--pointing_dir", type=str, help="Pointing directory that contains the spliced fits files")
    parser.add_argument("-t", "--threshold", type=float, default=5.0, help="The signal to noise threshold to stop at. Default = 5.0")
    parser.add_argument("-L", "--loglvl", type=str, default="INFO", help="Logger verbosity level. Default: INFO", choices=loglevels.keys())

    parser.add_argument_group("Batch Script Options:")
    parser.add_argument("--prevbins", type=int, default=None, help="The number of bins used in prepfold on the previous run. Not necessary for initial runs")

    args = parser.parse_args()

    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    ch.setLevel(loglevels[args.loglvl])
    formatter = logging.Formatter('%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    if args.pointing_dir == None:
        logger.error("No pointing directory supplied. Please specify the pointing directory path and rerun")
        sys.exit()


    find_optimal_prof(args.pointing_dir, args.prevbins, args.threshold, args.loglvl)

