import logging
from os.path import join
from glob import glob

from vcstools.config import load_config_file
from dpp.helper_config import from_yaml
from dpp.helper_checks import check_file_dir_exists

comp_config = load_config_file()
logger = logging.getLogger(__name__)


def cfg_status(psr_dir):
    """
    Checks a cfg to see how it ended and returns an int:
    100: Run completed with detection
    101: Run completed but RVM could not be fit
    200: No detections in the run
    201: cfg does not exist
    300: Something went wrong
    400: Something else happened
    """
    try:
        cfg = glob(join(psr_dir, "*.yaml"))[-1]
    except IndexError as e:
        status=201
    else:
        cfg = from_yaml(cfg)
        # Expected terminations
        milestones = [cfg["completed"][i] for i in cfg["completed"].keys()]
        if cfg["completed"]["classify"] and not cfg["completed"]["post_folds"]: # No detections
            status=200
        elif all(milestones): # Completed run
            status=100
        elif cfg["completed"]["RM"] and not cfg["completed"]["RVM_initial"]: # Bad paswing file
            status=101
        elif not cfg["run_ops"]["expected_finish"]: # Something went wrong
            status=300
        else: # I don't know what happened
            status=400
    return status


def opp_status(obsid):
    """Looks through all cfg files in obsid directory and returns dictionary on their status"""
    dpp_dir = join(comp_config["base_data_dir"], obsid, "dpp")
    check_file_dir_exists(dpp_dir)
    glob_cmd = join(dpp_dir, f"{obsid}*")
    psr_dirs = glob(glob_cmd)
    status_dict = {"100":[], "101":[], "200":[], "201":[], "300":[], "400":[]}
    for _dir in psr_dirs:
        status = cfg_status(_dir)
        psr = _dir.split("/")[-1].split("_")[-1]
        status_dict[str(status)].append(psr)
    return status_dict