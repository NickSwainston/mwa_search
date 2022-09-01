import logging
from os.path import join, isdir
from glob import glob

from vcstools.config import load_config_file
from vcstools.prof_utils import ProfileLengthError, NoFitError
from dpp.helper_config import from_yaml
from dpp.helper_checks import check_file_dir_exists
from dpp.helper_bestprof import NoUsableFoldsError
from dpp.helper_checks import InvalidPAFileError, FitsNotFoundError, PFDNotFoundError, PointingNotFoundError
from dpp.helper_classify import ClassifierFilesNotFoundError, InvalidClassifyBinsError

logger = logging.getLogger(__name__)

STATUS_DICT = {
        "100": "Run completed with detection",
        "101": "Run completed but RVM could not be fit due to unsifficient paswing file",
        "102": "Profile has insufficient bins to continue with Gaussian fit",
        "103": "Profile is too noisy to continue with Gaussian fit",
        "200": "No detections in the run",
        "201": "A PFD file is missing",
        "202": "Classifier files are missing",
        "203": "Bin count for classifer inputs is incorrect",
        "204": "Fits file is missing",
        "205": "Pointing directory is missing",
        "300": "Pipeline has not started",
        "400": "Something unknown went wrong"}

def message_from_status(status):
    """Takes a status code and returns the defualt error message for it"""
    if status not in STATUS_DICT.keys():
        raise ValueError(f"Invalid status: {status}. Valid statuses include: {STATUS_DICT.keys()}")
    return(STATUS_DICT[status])

def status_from_error(e):
    """Takes an error message and returns the designated status from it"""
    if type(e)   == InvalidPAFileError:
        status = "101"
    elif type(e) == ProfileLengthError:
        status = "102"
    elif type(e) == NoFitError:
        status = "103"
    elif type(e) == NoUsableFoldsError:
        status = "200"
    elif type(e) == PFDNotFoundError:
        status = "201"
    elif type(e) == ClassifierFilesNotFoundError:
        status = "202"
    elif type(e) == InvalidClassifyBinsError:
        status = "203"
    elif type(e) == FitsNotFoundError:
        status = "204"
    elif type(e) == PointingNotFoundError:
        status = "205"
    else: # Not sure what happened
        status = "400"
    return status


def cfg_status(psr_dir):
    """Checks a cfg to see how it ended and returns the status code"""
    cfg = glob(join(psr_dir, "*.yaml"))[-1]
    cfg = from_yaml(cfg)
    return cfg["run_ops"]["exit_status"], cfg["run_ops"]["detection"]


def opp_status(obsid):
    """Looks through all cfg files in obsid directory and returns dictionary on their status"""
    comp_config = load_config_file()
    dpp_dir = join(comp_config["base_data_dir"], obsid, "dpp")
    check_file_dir_exists(dpp_dir)
    glob_cmd = join(dpp_dir, f"{obsid}*")
    psr_dirs = glob(glob_cmd)
    run_status = {}
    detections = []
    for code in STATUS_DICT.keys():
        run_status[code] = []
    for _dir in psr_dirs:
        if not isdir(_dir):
            continue
        status, det = cfg_status(_dir)
        if det:
            detections.append(_dir)
        psr_dir = _dir.split("/")[-1]
        try:
            print(psr_dir)
            run_status[str(status)].append(psr_dir)
        except IndexError as e:
            logger.warn(f"Config file not found in directory: {_dir}")
            continue
    return run_status, detections
