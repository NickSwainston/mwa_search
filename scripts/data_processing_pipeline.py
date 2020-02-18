#!/usr/bin/env python3

import os
import logging
import argparse
import config
import glob
import sys
import datetime

from job_submit import submit_slurm
from mwa_metadb_utils import get_common_obs_metadata
import stokes_fold
import binfinder

logger = logging.getLogger(__name__)


#----------------------------------------------------------------------
class run_params_class:

    def __init__(self, pointing_dir=None, cal_id=None, obsid=None, pulsar=None,\
                threshold=10.0, stop=False, loglvl="INFO",\
                mwa_search="master", vcs_tools="master",\
                subint=10.0, RM=None, RM_err=None, stokes_bins=None,\
                nocrop=False, bestprof=None, archive=None, out_dir=None,\
                epndb_dir=None, beg=None, end=None, freq=None, stokes_dep=None,
                no_ephem=False, dspsr_ops="", prep_ops=""):

        #Obs inormation
        self.pointing_dir   = pointing_dir
        self.cal_id         = cal_id
        self.obsid          = obsid
        self.pulsar         = pulsar
        self.beg            = beg
        self.end            = end
        self.freq           = freq

        #Versions
        self.mwa_search     = mwa_search
        self.vcs_tools      = vcs_tools

        #Run Options
        self.stop           = stop
        self.loglvl         = loglvl
        self.stokes_dep     = stokes_dep
        self.no_ephem       = no_ephem
        self.dspsr_ops      = dspsr_ops
        self.prep_ops       = prep_ops

        #Plotting Options
        self.nocrop         = nocrop
        self.bestprof       = bestprof
        self.archive        = archive
        self.out_dir        = out_dir
        self.epndb_dir      = epndb_dir

        #Other Parameters
        self.threshold      = threshold
        self.subint         = subint
        self.RM             = RM
        self.RM_err         = RM_err
        self.stokes_bins    = stokes_bins

        if self.pointing_dir:
            if isinstance(self.pointing_dir, list) and len(self.pointing_dir)==1:
                self.pointing_dir=self.pointing_dir[0]
        if not self.freq and self.obsid:
            self.set_freq_from_metadata(obsid)


    def set_beg(self, beg):
        self.beg = beg

    def set_end(self, end):
        self.end = end

    def set_best_bins(self, bins):
        self.best_bins = bins

    def set_stokes_bins(self, bins):
        self.stokes_bins = bins

    def set_RM_and_err(self, RM, RM_err):
        self.RM = RM
        self.RM_err = RM_err

    def set_pointing_dir(self, new_dir):
        self.pointing_dir = new_dir

    def stop_now(self):
        self.stop=True

    def set_freq(self, new_freq):
        self.freq = new_freq

    def set_freq_from_metadata(self, obsid):
        self.freq = get_common_obs_metadata(obsid)[5]

    def clear_dependencies(self):
        self.stokes_dep = None

#----------------------------------------------------------------------
def binfinder_launch_line(run_params, dpp=False, single_pointing=None):

    if dpp:
        launch_line = "data_processing_pipeline.py"
    else:
        launch_line = "binfinder.py"

    if single_pointing:
        p = single_pointing
    elif isinstance(run_params.pointing_dir, list):
        p = ""
        for pointing in run_params.pointing_dir:
                p += "{} ".format(pointing)
    else:
        p=run_params.pointing_dir

    launch_line += " -d {0} -O {1} -p {2} -o {3} -L {4} --vcs_tools {5}\
                    --mwa_search {6}"\
                    .format(p, run_params.cal_id, run_params.pulsar,\
                    run_params.obsid, run_params.loglvl, run_params.vcs_tools, run_params.mwa_search)
    if run_params.freq:
        launch_line += " -f {}".format(run_params.freq)
    if run_params.beg:
       launch_line += " --beg {}".format(run_params.beg)
    if run_params.end:
        launch_line += " --end {}".format(run_params.end)
    if run_params.stokes_dep:
        launch_line += " --stokes_dep {}".format(run_params.stokes_dep)
    if run_params.no_ephem:
        launch_line += " --no_ephem"
    if run_params.dspsr_ops != "":
        launch_line += " --dspsr_ops '{}'".format(run_params.dspsr_ops)
    if run_params.prep_ops != "":
        launch_line += " --prep_ops '{}'".format(run_params.prep_ops)

    return launch_line

#----------------------------------------------------------------------
def prepfold_time_alloc(duration, nbins, npfact=2, ndmfact=3, nosearch=False):
    """
    Estimates the jobtime for prepfold jobs based on inputs

    Parameters:
    -----------
    duration: float
        The duration of the folding time in seconds
    nbins: int
        Then number of bins to fold over
    npfact: int
        OPTIONAL - prepfold's npfact option. Default: 2
    ndmfact: int
        OPTIONAL - prepfold's ndmfact option. Default: 3
    nosearch: boolean
        OPTIONAL - true if prepfold's nosearch tag is used. Default: False

    Returns:
    --------
    time: string
        The allocated time for the fold job as a string that can be passed to the slurm batch handler
    """
    time = 600
    time += nbins
    time += duration

    if not nosearch:
        time += (10*npfact)**2 * (10*ndmfact)**2 + 100*(nbins)**(1/2)
    if time > 86400:
        logger.warn("Estimation for prepfold time greater than one day")
        time = 86400
    time = str(datetime.timedelta(seconds = int(time)))

    return time

#----------------------------------------------------------------------
def stokes_launch_line(run_params, dpp=False, custom_pointing=None):

    if dpp:
        launch_line = "data_processing_pipeline.py"
    else:
        launch_line = "stokes_fold.py"

    if custom_pointing:
        p = custom_pointing
    else:
        p = run_params.pointing_dir

    launch_line += " -d {0} -p {1} -b {2} -s {3} -o {4} -L {5}\
                    --mwa_search {6} --vcs_tools {7}"\
                    .format(p, run_params.pulsar, run_params.stokes_bins,\
                    run_params.subint, run_params.obsid, run_params.loglvl,\
                    run_params.mwa_search, run_params.vcs_tools)
    if run_params.freq:
        launch_line += " -f {}".format(run_params.freq)
    if run_params.beg:
       launch_line += " --beg {}".format(run_params.beg)
    if run_params.end:
        launch_line += " --end {}".format(run_params.end)
    if run_params.stop:
        launch_line += " -S"
    if run_params.no_ephem:
        launch_line += " --no_ephem"
    if run_params.dspsr_ops != "":
        launch_line += " --dspsr_ops {}".format(run_params.dspsr_ops)

    return launch_line

#----------------------------------------------------------------------
def copy_data(data_path, target_directory):
    """
    Copies the data path file to the targeted directory. Will make targeted directory if it doesn't exist

    Parameters:
    -----------
    data_path: str
        The path to copy
    target_directory: str
        The path of the directory to copy the data to
    """
    os.makedirs(target_directory, exist_ok=True)
    try:
        os.popen("cp {0} {1}".format(data_path, target_directory))
    except RuntimeError as error:
        logger.warning("File:{0} could not be copied to {1}".format(data_path, target_directory))
        logger.warning("Error message: {0}".format(error))

#----------------------------------------------------------------------
def stokes_fold(run_params):
    """
    Launches the stokes_fold part of the data processing pipeling

    Parameters:
    -----------
    run_params: object
        The run_params object defined by data_processing_pipeline
    """
    launch_line = stokes_launch_line(run_params)
    commands=[launch_line]
    name="Stokes_Fold_init_{0}_{1}".format(run_params.pulsar, run_params.obsid)
    comp_config = config.load_config_file()
    batch_dir = "{0}{1}/batch/".format(comp_config['base_product_dir'], run_params.obsid)

    job_id = submit_slurm(name, commands,\
                        batch_dir=batch_dir,\
                        slurm_kwargs={"time": "00:10:00"},\
                        module_list=["mwa_search/{0}".format(run_params.mwa_search),\
                                    "dspsr/master", "psrchive/master"],\
                        submit=True, vcstools_version="{0}".format(run_params.vcs_tools))
    logger.info("Job successfully submitted: {0}".format(name))
    return job_id

#----------------------------------------------------------------------
def binfind(run_params):
    """
    Launches the binfinding part of the data processing pipeline

    Parameters:
    -----------
    run_params: object
        The run_params object defined by data_proces_pipeline
    """
    launch_line = binfinder_launch_line(run_params)
    commands = [launch_line]
    #decide how much time to allocate based on number of pointings
    n_pointings = len(run_params.pointing_dir)
    if n_pointings<100:
        time = "00:30:00"
    elif n_pointings<400:
        time = "02:00:00"
    elif n_pointings<1000:
        time = "05:00:00"
    else:
        time = "10:00:00"

    name = "bf_initiate_{0}_{1}".format(run_params.pulsar, run_params.obsid)
    logger.info("Submitting binfinder script:")
    logger.info("")
    logger.info("Job Name: {}".format(name))
    comp_config = config.load_config_file()
    batch_dir = "{0}{1}/batch/".format(comp_config['base_product_dir'], run_params.obsid)
    job_id = submit_slurm(name, commands,\
                        batch_dir=batch_dir,\
                        slurm_kwargs={"time": time},\
                        module_list=['mwa_search/{0}'.format(run_params.mwa_search),\
                                    'presto/no-python'],\
                        submit=True, vcstools_version="{0}".format(run_params.vcs_tools))

    logger.info("Job successfully submitted: {0}".format(name))
    return job_id

def work_out_what_to_do(run_params):
    """

    Parameters:
    -----------
    run_params: object
        The run_params object defined by data_proces_pipeline
    """
    if isinstance(run_params.pointing_dir, list):
        logger.info("Multiple pointing directories")
        binfind(run_params)
    else:
        fits_files_in_dir = glob.glob("{}/*.fits".format(run_params.pointing_dir))
        bestprof_files_in_dir = glob.glob("{}/*.bestprof".format(run_params.pointing_dir))
        archive_files_in_dir = glob.glob("{}/*.ar2".format(run_params.pointing_dir))
        pol_plot_in_dir = glob.glob("{}/*Polarimetry*.png".format(run_params.pointing_dir))
        if fits_files_in_dir and not bestprof_files_in_dir:
            logger.info("Launching binfinder")
            binfind(run_params)
        elif fits_files_in_dir and (not archive_files_in_dir or not pol_plot_in_dir):
            logger.info("Launching stokes_fold")
            stokes_fold(run_params)
        else:
            logger.info("Nothing left to do...")

#----------------------------------------------------------------------
if __name__ == '__main__':

    loglevels = dict(DEBUG=logging.DEBUG,
                    INFO=logging.INFO,
                    WARNING=logging.WARNING,
                    ERROR = logging.ERROR)

    parser = argparse.ArgumentParser(description="""A pipeline for processing calibrated VCS data""",\
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    obsop = parser.add_argument_group("Observation Options")
    obsop.add_argument("-d", "--pointing_dir", nargs='+', type=str, help="The location of the pointing directory/s")
    obsop.add_argument("-o", "--obsid", type=str, help="The obs ID of the data")
    obsop.add_argument("-O", "--cal_id", type=str, help="The ID of the calibrator used to calibrate the data")
    obsop.add_argument("-p", "--pulsar", type=str, help="The J name of the pulsar. e.g. J2241-5236")
    obsop.add_argument("--beg", type=int, help="The beginning of the observation")
    obsop.add_argument("--end", type=int, help="The end of the observation")
    obsop.add_argument("-f", "--freq", type=float, help="The central frequency of the observation in MHz")

    binfindop = parser.add_argument_group("Binfinder Options")
    binfindop.add_argument("-t", "--threshold", type=float, default=10.0, help="The presto sigma value\
                             above which is deemed a detection. If this value is not exceeded in any\
                             of the folds, the pipeline will terminate")
    binfindop.add_argument("--prep_ops", type=str, default="", help="Provide as a string in quotes any prepfold command you would like to use for folding.\
                        eg: ' -dm 50.0 -p 0.50625' (make sure there is a space before the first argument))")

    stokesop = parser.add_argument_group("Stokes Fold Options")
    stokesop.add_argument("-b", "--nbins", type=int, help="The number of bins for to fold over for the stokes folding script")
    stokesop.add_argument("-s", "--subint", type=float, default=10.0, help="The length of the integrations (in seconds) used for dspsr.")
    stokesop.add_argument("--dspsr_ops", type=str, default="", help="Provide as a string in quotes any dspsr command you would like to use for folding.\
                        eg: ' -D 50.0 -c 506.25' (make sure there is a space before the first argument)")

    otherop = parser.add_argument_group("Other Options")
    otherop.add_argument("--no_ephem", action="store_true", help="Use this to override the use of an ephemeris for foldign the pulsar")
    otherop.add_argument("-L", "--loglvl", type=str, default="INFO", help="Logger verbosity level", choices=loglevels.keys())
    otherop.add_argument("-S", "--stop", action="store_true", help="Use this mode to tell the pipeline not to continue processing data after finishing the desired task")
    otherop.add_argument("--mwa_search", type=str, default="master",  help="The version of mwa_search to use")
    otherop.add_argument("--vcs_tools", type=str, default="master", help="The version of vcs_tools to use")

    args = parser.parse_args()

    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    ch.setLevel(loglevels[args.loglvl])
    formatter = logging.Formatter('%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.propagate = False

    #Make assertions
    if not args.pointing_dir:
        logger.error("Pointing directory must be supplied")
        sys.exit(1)
    if not args.obsid:
        logger.error("Observtion ID must be supplied")
        sys.exit(1)
    if not args.pulsar:
        logger.error("Pulsar name must be supplied")
        sys.exit(1)
    if not args.beg or not args.end:
        logger.error("Beginning and end times must be supplied")
        sys.exit(1)

    run_params = run_params_class(pointing_dir=args.pointing_dir, cal_id=args.cal_id,\
                                pulsar=args.pulsar, obsid=args.obsid, stop=args.stop,\
                                mwa_search=args.mwa_search,\
                                vcs_tools=args.vcs_tools, loglvl=args.loglvl,\
                                threshold=args.threshold, stokes_bins=args.nbins,\
                                subint=args.subint, beg=args.beg, end=args.end, freq=args.freq,
                                dspsr_ops=args.dspsr_ops, prep_ops=args.prep_ops,\
                                no_ephem=args.no_ephem)

    work_out_what_to_do(run_params)

