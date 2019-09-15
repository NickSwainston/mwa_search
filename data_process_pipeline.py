#!/usr/bin/env python3

import os
import logging
import argparse
import socket
from job_submit import submit_slurm
import config

logger = logging.getLogger(__name__)


#----------------------------------------------------------------------
class run_params_class:

    def __init__(self, pointing_dir=None, cal_id=None,obsid=None, pulsar=None,\
                threshold=10.0, stop=False, next_mode=True, loglvl="INFO",\
                mode=None, mwa_search="master", vcs_tools="multi-pixel_beamform",\
                nbins=None, subint=10.0, RM=None, RM_err=None, prevbins=None,\
                best_bins=None, force_initial=False, nocrop=False, bestprof=None,\
                archive=None, out_dir=None, epndb_dir=None):

        #Obs inormation
        self.pointing_dir   = pointing_dir
        self.cal_id         = cal_id
        self.obsid          = obsid
        self.pulsar         = pulsar

        #Versions
        self.mwa_search     = mwa_search
        self.vcs_tools      = vcs_tools

        #Run Options
        self.stop           = stop
        self.loglvl         = loglvl
        self.force_initial  = force_initial
        self.mode           = mode

        #Plotting Options
        self.nocrop         = nocrop
        self.bestprof       = bestprof
        self.archive        = archive
        self.out_dir        = out_dir
        self.epndb_dir      = epndb_dir        

        #Other Parameters
        self.threshold      = threshold
        self.nbins          = nbins
        self.subint         = subint
        self.RM             = RM
        self.RM_err         = RM_err
        self.prevbins       = prevbins
        self.best_bins      = best_bins


        if self.obsid==None:
            self.obsid=info_from_dir(self.pointing_dir)["obsid"]
        if self.pointing_dir is not None:
            if len(self.pointing_dir)==1:
                self.pointing_dir=self.pointing_dir[0]


    def set_prevbins(self, prevbins):
        self.prevbins = prevbins

    def set_best_bins(self, bins):
        self.best_bins = bins

    def set_RM_and_err(self, RM, RM_err):
        self.RM = RM
        self.RM_err = RM_err

    def set_pointing_dir(self, new_dir):
        self.pointing_dir = new_dir

    def stop_now(self):
        self.stop=True

#----------------------------------------------------------------------
def copy_data(data_path, target_directory):
    #copies the data_path file to target_directory
    #Make the target directory if needed
    os.makedirs(target_directory, exist_ok=True)
    try:
        os.popen("cp {0} {1}".format(data_path, target_directory))
    except RuntimeError as error:
        logger.warning("File:{0} could not be copied to {1}".format(data_path, target_directory))
        logger.warning("Error message: {0}".format(error)) 

#----------------------------------------------------------------------
def info_from_dir(pointing_dir):

    #given a pointing directory, returns a dictionary containing the obsid and pulsar name
    mydict = {}
    pointing_dir = pointing_dir.split("/")
    for i in pointing_dir:
        if i == "":
            pointing_dir.remove(i)

    mydict["obsid"] = pointing_dir[3]

    idx = pointing_dir.index("pointings")
    pulsar_coords = pointing_dir[idx+1]
    pulsar_1 = pulsar_coords.split(":")[0] + pulsar_coords.split(":")[1]
    pulsar_2 = pulsar_coords.split(":")[2].split("_")[1] + pulsar_coords.split(":")[3]
    mydict["pulsar"] = "J" + pulsar_1 + pulsar_2

    return mydict

#----------------------------------------------------------------------
def stokes_fold(run_params):

    launch_line = "stokes_fold.py -m i -d {0} -p {1} -b {2} -s {3} -L {4} --vcs_tools {5} --mwa_search {6}"\
                .format(run_params.pointing_dir, run_params.pulsar, run_params.nbins, run_params.subint,\
                run_params.loglvl, run_params.vcs_tools, run_params.mwa_search)
    if run_params.stop==True:
        launch_line += " -S"

    commands=[]
    commands.append(launch_line)

    name="Stokes_Fold_{0}_{1}".format(run_params.pulsar, run_params.obsid)
    comp_config = config.load_config_file()
    batch_dir = "{0}{1}/batch/".format(comp_config['base_product_dir'], run_params.obsid)

    submit_slurm(name, commands,\
                batch_dir=batch_dir,\
                slurm_kwargs={"time": "00:02:00"},\
                module_list=["mwa_search/{0}".format(run_params.mwa_search),\
                            "dspsr/master", "psrchive/master"],\
                submit=True, vcstools_version="{0}".format(run_params.vcs_tools))


#----------------------------------------------------------------------
def binfind(run_params):

    launch_line = "binfinder.py -O {0} -t {1} -p {2} -o {3} -L {4} --mwa_search {5}\
                --vcs_tools {6}"\
                .format(run_params.cal_id, run_params.threshold, run_params.pulsar,\
                run_params.obsid, run_params.loglvl, run_params.mwa_search,\
                run_params.vcs_tools)

    if run_params.stop==True:
        launch_line += " -S"

    #Run binfinder.py
    if run_params.mode=='f':
        launch_line += " -d {0}".format(run_params.pointing_dir)
        launch_line += " -m f"
    elif run_params.mode=='p':
        launch_line += " -d {0}".format(run_params.pointing_dir)
        launch_line += " -m c"
    elif run_params.mode=="m":
        pointing_string=""
        for p in run_params.pointing_dir:
            logger.info("folding on: {0}".format(p))
            pointing_string = pointing_string + p + " "
        launch_line += " -d {0}".format(pointing_string)
        launch_line += " -m m"

    commands = []
    commands.append("echo 'Submitting binfinder in mode {0}'".format(run_params.mode))
    commands.append(launch_line)

    #decide how much time to allocate based on number of poitnigns
    n_pointings = len(run_params.pointing_dir)
    if n_pointings<100:
        time = "00:30:00"
    elif n_pointings<400:
        time = "02:00:00"
    elif n_pointings<1000:
        time = "05:00:00"
    else:
        time = "10:00:00"

    name = "binfind_initiate_{0}_{1}".format(run_params.pulsar, run_params.obsid)
    comp_config = config.load_config_file()
    batch_dir = "{0}{1}/batch/".format(comp_config['base_product_dir'], run_params.obsid)
    submit_slurm(name, commands,\
                batch_dir=batch_dir,\
                slurm_kwargs={"time": time},\
                module_list=['mwa_search/{0}'.format(run_params.mwa_search),\
                            'presto/no-python'],\
                submit=True, vcstools_version="{0}".format(run_params.vcs_tools))

    logger.info("Job successfully submitted")

#----------------------------------------------------------------------
if __name__ == '__main__':
    #dictionary for choosing log-levels
    loglevels = dict(DEBUG=logging.DEBUG,
                    INFO=logging.INFO,
                    WARNING=logging.WARNING,
                    ERROR = logging.ERROR)


    #Arguments

    parser = argparse.ArgumentParser(description="""A pipeline for processing calibrated VCS data""")

    obsop = parser.add_argument_group("Observation Options")
    obsop.add_argument("-d", "--pointing_dir", nargs='+', help="The location of the pointing directory/s")
    obsop.add_argument("-o", "--obsid", type=str, help="The obs ID of the data")
    obsop.add_argument("-O", "--cal_id", type=str, help="The ID of the calibrator used to calibrate the data")
    obsop.add_argument("-p", "--pulsar", type=str, help="The J name of the pulsar. e.g. J2241-5236")

    binfindop = parser.add_argument_group("Binfinder Options")
    binfindop.add_argument("-t", "--threshold", type=float, default=10.0, help="The presto sigma value\
                             above which is deemed a detection. If this value is not exceeded in any\
                             of the folds, the pipeline will terminate. However, a value of 5.0 will\
                             always be checked and uploaded to the pulsar database if it exists.\
                             Default: 10.0")

    stokesop = parser.add_argument_group("Stokes Fold Options")
    stokesop.add_argument("-b", "--nbins", type=int, default=128, help="The number of bins for to fold over for the stokes folding script. Default: 128")
    stokesop.add_argument("-s", "--subint", type=float, default=10.0, help="The length of the integrations (in seconds) used for dspsr. Default: 10.0")

    otherop = parser.add_argument_group("Other Options")
    otherop.add_argument("-L", "--loglvl", type=str, default="INFO", help="Logger verbosity level. Default: INFO", choices=loglevels.keys())
    otherop.add_argument("-S", "--stop", action="store_true", help="Use this mode to tell the pipeline not to continue processing data after finishing the desired task")
    otherop.add_argument("--mwa_search", type=str, default="master",  help="The version of mwa_search to use. Default: master")
    otherop.add_argument("--vcs_tools", type=str, default="multi-pixel_beamform", help="The version of vcs_tools to use. Default: master")
    modeop = parser.add_argument_group("Mode Options")
    otherop.add_argument("-m", "--mode", type=str, help="The mode in which to run this script:\n\
                        Binfinder Options:\n\
                        'p' - Folds on a small number of bins in order to check if a pulsar is\n\
                         detected in the given pointing directory (runs the b mode after by default)\n\
                        'm' - Folds on  any numer of pointings and finds the one with the best detection\
                        (runs the b mode after by default)\n\
                        'b' - attempts to find an adequate number of bins to fold the pulsar with\n\
                         and outputs a bestprof file (runs the s mode after by default)\n\
                        Stokes Fold Options:\n\
                        's' - will fold across stokes IQUV and attempt to find the rotation measure")


    args = parser.parse_args()

    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    ch.setLevel(loglevels[args.loglvl])
    formatter = logging.Formatter('%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.propagate = False

    run_params = run_params_class(pointing_dir=args.pointing_dir, cal_id=args.cal_id,\
                                pulsar=args.pulsar, obsid=args.obsid, stop=args.stop,\
                                mode=args.mode, mwa_search=args.mwa_search,\
                                vcs_tools=args.vcs_tools, loglvl=args.loglvl,\
                                threshold=args.threshold, nbins=args.nbins,\
                                subint=args.subint)

    if run_params.mode=="b" or run_params.mode=="p" or run_params.mode=="m":
        binfind(run_params)
    elif run_params.mode=="s":
        stokes_fold(run_params)
    else:
        logger.error("Mode not recognized. Please rerun with a valid mode identifer")
