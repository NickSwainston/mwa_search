#!/usr/bin/env python3

import os
import logging
import argparse
from config_vcs import load_config_file
import glob
import sys
import shutil

from job_submit import submit_slurm
from mwa_metadb_utils import get_common_obs_metadata
import stokes_fold
import binfinder
import submit_to_database as std

logger = logging.getLogger(__name__)

#----------------------------------------------------------------------
def binfinder_launch_line(run_params, dpp=False):
    """
    Creates a launch command using the run_params class

    Parameters:
    -----------
    run_params: object
        The run_params object
    dpp: boolean
        OPTIONAL - If True, will launch the data_processing_pipeline with the run_params variables instead of binfinder.py. Default: False

    Returns:
    --------
    launch_line: str
        The launch command
    """
    if dpp:
        launch_line = "data_processing_pipeline.py"
    else:
        launch_line = "binfinder.py"

    if isinstance(run_params.pointing_dir, list):
        p = ""
        for pointing in run_params.pointing_dir:
                p += "{} ".format(pointing)
    else:
        p=run_params.pointing_dir

    launch_line += " -d {0} -O {1} -p {2} -o {3} -L {4} --vcs_tools {5} --mwa_search {6}"\
                    .format(p, run_params.cal_id, run_params.pulsar,\
                    run_params.obsid, run_params.loglvl, run_params.vcs_tools, run_params.mwa_search)
    if run_params.freq:
        launch_line += " -f {}".format(run_params.freq)
    if run_params.beg:
       launch_line += " --beg {}".format(run_params.beg)
    if run_params.end:
        launch_line += " --end {}".format(run_params.end)
    if run_params.dm:
        launch_line += " --dm {}".format(run_params.dm)
    if run_params.period:
        launch_line += " --period {}".format(run_params.period)
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
def stokes_launch_line(run_params, dpp=False, custom_pointing=None):
    """
    Creates a launch command using the run_params class

    Parameters:
    -----------
    run_params: object
        The run_params object
    dpp: boolean
        OPTIONAL - If True, will launch the data_processing_pipeline with the run_params variables instead of stokes_fold.py. Default: False
    custom_pointing: str
        OPTIONAL - A custom pointing directory to run from

    Returns:
    --------
    launch_line: str
        The launch command
    """
    if dpp:
        launch_line = "data_processing_pipeline.py"
    else:
        launch_line = "stokes_fold.py"

    if custom_pointing:
        p = custom_pointing
    else:
        p = run_params.pointing_dir

    launch_line += " -d {0} -p {1} -L {2} -o {3} -O {4} --mwa_search {5} --vcs_tools {5}"\
                    .format(p, run_params.pulsar, run_params.loglvl, run_params.obsid,
                    run_params.cal_id, run_params.mwa_search, run_params.vcs_tools)

    if run_params.stokes_bins:
        launch_line += " -b {}".format(run_params.stokes_bins)
    if run_params.subint:
        launch_line += " -s {}".format(run_params.subint)
    if run_params.freq:
        launch_line += " -f {}".format(run_params.freq)
    if run_params.beg:
       launch_line += " --beg {}".format(run_params.beg)
    if run_params.end:
        launch_line += " --end {}".format(run_params.end)
    if run_params.dm:
        launch_line += " --dm {}".format(run_params.dm)
    if run_params.period:
        launch_line += " --period {}".format(run_params.period)
    if run_params.stop:
        launch_line += " -S"
    if run_params.no_ephem:
        launch_line += " --no_ephem"
    if run_params.dspsr_ops != "":
        launch_line += " --dspsr_ops {}".format(run_params.dspsr_ops)
    if run_params.cand:
        launch_line += " --cand"
    if run_params.rvmres:
        launch_line += " --rvmres {}".format(run_params.rvmres)

    return launch_line

#----------------------------------------------------------------------
def prepfold_time_alloc(prepfold_dict, beg, end):
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
    nopsearch = False
    nopdsearch = False
    nodmsearch = False
    nosearch = False
    if "nopsearch" in prepfold_dict:
        nopsearch = True
    if "nodmsearch" in prepfold_dict:
        nodmsearch = True
    if "nopdsearch" in prepfold_dict:
        nopdsearch = True
    if "nosearch" in prepfold_dict:
        nosearch = True
    npfact = prepfold_dict["npfact"]
    ndmfact = prepfold_dict["ndmfact"]
    nbins = prepfold_dict["n"]
    duration = (prepfold_dict["end"] - prepfold_dict["start"]) * (end - beg)

    time = 600
    time += nbins
    time += duration

    if not nosearch:
        ptime = 1
        pdtime = 1
        dmtime = 1
        if not nopsearch:
            ptime = npfact*nbins
        if not nopdsearch:
            pdtime = npfact*nbins
        if not nodmsearch:
            dmtime = ndmfact*nbins
        time += ((ptime * pdtime * dmtime)/1e4)
    time = time*2 #compute time is very sporadic so just give double the allocation time

    return time

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
        shutil.copy(data_path, target_directory)
    except RuntimeError as error:
        logger.warning("File:{0} could not be copied to {1}".format(data_path, target_directory))
        logger.warning("Error message: {0}".format(error))

def upload_formatted_file(filename, obsid, pulsar, bins, cal_id, filetype, name_info="", extension=""):
    """
    Creates a new filename and uploads an archive file to the pulsar database if the same named
    file does not already exist on the database

    Parameters:
    -----------
    archive: string
        The name of the archive file to upload
    obsid: int
        The observation ID
    pulsar: string
        Then name of the pulsar
    bins: int
        The Number of bins
    cal_id: int
        The calibrator ID
    filetype: int
        The type of file to upload: 1: Archive, 2: Timeseries, 3: Diagnostics, 4: Calibration Solution, 5: Bestprof
    name_info: str
        OPTIONAL - additional info to add to the name of the uploaded file. Default: ''
    extention: str
        OPTIONAL - The file extension of the uploaded file. Default: ''
    """
    all_ftypes = std.get_filetypes_from_db(obsid, pulsar, filetype)
    fname_pref = std.filename_prefix(obsid, pulsar, bins=bins, cal=cal_id)
    upname = "{}".format(fname_pref)
    upname += name_info
    upname += extension

    metadata = get_common_obs_metadata(obsid)
    subbands = std.get_subbands(metadata)

    if os.path.basename(upname) not in all_ftypes:
        logger.info("Archive file not on databse. Uploading...")
        shutil.copy(filename, upname)
        std.upload_file_to_db(obsid, pulsar, upname, filetype, metadata=metadata, coh=True)
        os.remove(upname)
    else:
        logger.info("file on database. Not uploading")

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
    comp_config = load_config_file()
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
    comp_config = load_config_file()
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

    foldop = parser.add_argument_group("Folding/processing Options")
    foldop.add_argument("-t", "--threshold", type=float, default=8.0, help="The presto sigma value\
                             above which is deemed a detection. If this value is not exceeded in any\
                             of the folds, the pipeline will terminate")
    foldop.add_argument("--prep_ops", type=str, default="", help="Provide as a string in quotes any prepfold command you would like to use for folding.\
                        eg: ' -dm 50.0 -p 0.50625' (make sure there is a space before the first argument))")
    foldop.add_argument("--dm", type=float, default=None, help="The dispersion measure to fold around")
    foldop.add_argument("--period", type=float, default=None, help="The period to fold around in milliseconds")
    foldop.add_argument("-b", "--nbins", type=int, help="The number of bins for to fold over for the stokes folding script")
    foldop.add_argument("-s", "--subint", type=float, default=None, help="The length of the integrations (in seconds) used for dspsr.")
    foldop.add_argument("--dspsr_ops", type=str, default="", help="Provide as a string in quotes any dspsr command you would like to use for folding.\
                        eg: ' -D 50.0 -c 506.25' (make sure there is a space before the first argument)")
    foldop.add_argument("--rvmres", type=int, default=90, help="The number of degree samples to try for alpha and beta.")

    otherop = parser.add_argument_group("Other Options")
    otherop.add_argument("--no_ephem", action="store_true", help="Use this to override the use of an ephemeris for foldign the pulsar")
    otherop.add_argument("-L", "--loglvl", type=str, default="INFO", help="Logger verbosity level", choices=loglevels.keys())
    otherop.add_argument("-S", "--stop", action="store_true", help="Use this mode to tell the pipeline not to continue processing data after finishing the desired task")
    otherop.add_argument("--mwa_search", type=str, default="master",  help="The version of mwa_search to use")
    otherop.add_argument("--vcs_tools", type=str, default="master", help="The version of vcs_tools to use")
    otherop.add_argument("--cand", action="store_true", help="use this tag if this is not a kown pulsar")

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

    rp={}
    rp["pointing_dir"]      = args.pointing_dir
    rp["cal_id"]            = args.cal_id
    rp["pulsar"]            = args.pulsar
    rp["obsid"]             = args.obsid
    rp["stop"]              = args.stop
    rp["mwa_search"]        = args.mwa_search
    rp["vcs_tools"]         = args.vcs_tools
    rp["loglvl"]            = args.loglvl
    rp["threshold"]         = args.threshold
    rp["stokes_bins"]       = args.nbins
    rp["beg"]               = args.beg
    rp["end"]               = args.end
    rp["freq"]              = args.freq
    rp["dspsr_ops"]         = args.dspsr_ops
    rp["prep_ops"]          = args.prep_ops
    rp["no_ephem"]          = args.no_ephem
    rp["dm"]                = args.dm
    rp["period"]            = args.period
    rp["cand"]              = args.cand
    rp["rvmres"]            = args.rvmres
    if args.subint:
        rp["subint"] = args.subint
    run_params = run_params_class(**rp)

    work_out_what_to_do(run_params)

