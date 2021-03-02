import logging
from os.path import exists, join
from glob import glob
import numpy as np

from dpp.helper_prepfold import generate_prep_name
from dpp.helper_terminate import finish_unsuccessful

logger = logging.getLogger(__name__)


class InvalidFileError(Exception):
    """Raise when a file is not valid for any reason"""
    pass


def check_pipe_integrity(cfg):
    """Runs the required checks for the pipeline's current task"""
    # Unconditional:
    check_pfds(cfg)

    # Confitional
    if cfg["completed"]["init_folds"] == False:
        check_all_beamformed_fits(cfg)
    elif cfg["completed"]["classify"] == False:
        check_all_beamformed_fits(cfg)
        check_file_dir_exists(cfg["files"]["classify_dir"])
    elif cfg["completed"]["post_folds"] == False:
        check_all_beamformed_fits(cfg)
    elif cfg["completed"]["upload"] == False:
        check_all_beamformed_fits(cfg)
    elif cfg["completed"]["debase"] == False:
        check_file_dir_exists(cfg["files"]["archive"])
        check_file_dir_exists(cfg["files"]["converted_fits"])
    elif cfg["completed"]["RM"] == False:
        check_file_dir_exists(cfg["files"]["debased_fits"])
    elif cfg["completed"]["RVM_initial"] == False:
        check_file_dir_exists(cfg['files']['paswing'])
        try:
            check_paswing(cfg['files']['paswing'])
        except InvalidFileError as e:
            finish_unsuccessful(cfg, e)
    elif cfg["completed"]["RVM_final"] == False:
        check_file_dir_exists(cfg["files"]["RVM_fit_initial"])


def check_all_beamformed_fits(cfg):
    """
    Checks that fits files exist in all pointings
    Raises FileNotFoundError
    """
    psr_dir = cfg["files"]["psr_dir"]
    for pointing in cfg["folds"].keys():
        pointing_dir = join(psr_dir, pointing)
        fits_in_pointing_dir = join(pointing_dir, "*.fits")
        if not glob(fits_in_pointing_dir):
            raise FileNotFoundError(f".fits files not found in pointing directory: {pointing_dir}")


def check_pfds(cfg):
    """
    Checks that all expected .pfd files exist
    Raises FileNotFoundError
    """
    psr_dir = cfg["files"]["psr_dir"]
    names = []
    if cfg["completed"]["init_folds"] == True:
        for pointing in cfg["folds"].keys():
            for bins in cfg["folds"][pointing]["init"].keys():
                names.append(f"{generate_prep_name(cfg, bins, pointing)}*.pfd")
    if cfg["completed"]["post_folds"] == True:
        my_pointing = cfg["source"]["my_pointing"]
        for bins in cfg["folds"][my_pointing]["post"].keys():
            names.append(f"{generate_prep_name(cfg, bins, my_pointing)}*.pfd")
    for name in names:
        if not glob(name):
            raise FileNotFoundError(f"Expected pfd file not found {name}")


def check_file_dir_exists(file_dir):
    """
    Checks if a file or directory exists
    Raises FileNotFoundError
    """
    if not exists(file_dir):
        raise FileNotFoundError(f"Expected file or directory does not exist {file_dir}")


def check_paswing(paswing_file):
    """
    Checks existence of paswing file and if all PA points in a paswing file are identical.
    Raises InvalidFileError
    """
    paswing = np.loadtxt(paswing_file)
    pa = [i[-2] for i in paswing] # All PA points
    if len(list(set(pa))) < 5: # Check if there are enough usable (nonzero) PA points
        raise InvalidFileError(f"Not enough usable PA points in file {paswing_file}")