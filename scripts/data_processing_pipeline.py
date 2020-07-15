#!/usr/bin/env python3

import os
from os import path.join as ospj
import logging
import argparse
from config_vcs import load_config_file
import glob
import sys
import shutit

from job_submit import submit_slurm
from mwa_metadb_utils import get_common_obs_metadata
import stokes_fold
from binfinder import bf_main
import submit_to_database as std
import pipe_helper
import yaml_helper
import dpp_check_args

logger = logging.getLogger(__name__)


def get_next_name(pipe):
    if not pipe["completed"]["bf"]:
        if not pipe["completed"]["post"]:
            name = "post"
        elif not pipe["completed"]["submit_move"]:
            name = "submit_move"
    elif not polarimetry["completed"]:
        name = "Polarimetry"
    return name

def resubmit_self(pipe, dependencies=None):
    batch_dir = os.path.join(comp_config['base_product_dir'], pipe["obs"]["id"], "batch")
    func = name = get_next_name(pipe)
    batch_name = f"dpp_{name}_{pipe['source']['name']}_{pipe['obs']['id']}"
    yaml_name = ospj(pipe["run_ops"]["dir"], f"{pipe['obs']['id']}_{pipe['source']['name']}.yaml")
    yaml_helper.dump_to_yaml(pipe, yaml_name)
    commands = [f"data_processing_pipeline.py --yaml {yaml_name}"]
    job_id = submit_slurm(batch_name, commands,\
                batch_dir=batch_dir,
                slurm_kwargs={"time": "00:30:00"},
                module_list=[f"mwa_search/{pipe['run_ops']['mwa_search']}", f"vcstools/{pipe['run_ops']['vcstools']}"],
                submit=True, depend=dependencies, depend_type="afterok")

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


def main(kwargs):
    kwargs = dpp_check_args.yaml_check_args(kwargs)
    os.chdir(kwargs["run_dir"])
    pipe = pipe_helper.initiate_pipe(kwargs)
    if not pipe["completed"]["bf"]:
        bf_main(pipe)
    elif not pipe["completed"]["submit_move"]:
        submit_move_main(pipe)
    elif not pipe["completed"]["polarimetry"]:
        pass #TODO: make the pol pipe


if __name__ == '__main__':

    loglevels = dict(DEBUG=logging.DEBUG,
                    INFO=logging.INFO,
                    WARNING=logging.WARNING,
                    ERROR = logging.ERROR)

    parser = argparse.ArgumentParser(description="""A pipeline for processing calibrated VCS data""",\
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    yaml = parser.add_argument_group("YAML File")
    yaml.add_argument("--yaml", type="str", help="The pathname of the .yaml file to work from")

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
    foldop.add_argument("--dm", type=float, default=None, help="The dispersion measure to fold around")
    foldop.add_argument("--period", type=float, default=None, help="The period to fold around in milliseconds")
    foldop.add_argument("-b", "--nbins", type=int, help="The number of bins for to fold over for the stokes folding script")
    foldop.add_argument("-s", "--subint", type=float, default=None, help="The length of the integrations (in seconds) used for dspsr.")
    foldop.add_argument("--dspsr_ops", type=str, default="", help="Provide as a string in quotes any dspsr command you would like to use for folding.\
                        eg: ' -D 50.0 -c 506.25' (make sure there is a space before the first argument)")
    foldop.add_argument("--rvmres", type=int, default=90, help="The number of degree samples to try for alpha and beta.")
    foldop.add_argument("--mask", type=str, help="The pathname of the mask to use for folding")

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
    kwargs = vars(args)
    main(kwargs)