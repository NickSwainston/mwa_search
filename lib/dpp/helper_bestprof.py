#!/usr/bin/env python3
import logging
import numpy as np
from glob import glob
from os.path import join

from dpp.helper_config import from_yaml, dump_to_yaml
from dpp.helper_files import glob_pfds


logger = logging.getLogger(__name__)


class NoUsableFolds(Exception):
    """Raise when no usable folds are found in a pipe"""
    def __init__(self, *args):
        if args:
            self.message = args[0]
        else:
            self.message = ""

def bestprof_info(filename):
    """
    Finds various information on a .bestprof file
    Parameters:
    filename: string
        The path of the bestprof file
    Returns:
    info_dict: dictionary
        A dictionary consisting of the following:
        obsid: int
            The ID of the observation
        puslar: string
            The J name of the pulsar
        nbins: int
            The number of bins used to fold this profile
        chi: float
            The reduced Chi squared value of the fold
        sn: float
            The signal to noise ratio of the fold
        dm: float
            The pulsar's dispersion measure
        period: float
            The pulsar's period
        period_error: float
            The error in the pulsar's period measurement
    """
    #open the file and read the info into a dictionary
    info_dict = {}
    f = open(filename, "r")
    lines = f.read()
    f.close()
    lines = lines.split("\n")
    #info:
    #info_dict["obsid"] = int(lines[0].split()[4].split("_")[0])
    #info_dict["pulsar"] = lines[1].split()[3].split("_")[1]
    info_dict["nbins"] = int(lines[9].split()[4])
    info_dict["chi"] = float(lines[12].split()[4])
    info_dict["sn"] = float(lines[13].split()[4][2:])
    info_dict["dm"] = float(lines[14].split()[4])
    info_dict["period"] = float(lines[15].split()[4])/1e3 #in seconds
    info_dict["period_error"] = float(lines[15].split()[6])/1e3
    info_dict["profile"] = list(np.genfromtxt(filename)[:,1])
    f.close()
    return info_dict


def find_best_pointing(cfg):
    """Decides the best folding solution from bestprofs and updates cfg"""
    # Populate cfg with initial fold info
    for pointing in cfg["folds"].keys():
        nbins = list(cfg["folds"][pointing]["init"].keys())[0]
        try:
            bestprof_name = glob_pfds(cfg, pointing, nbins, pfd_type="pfd.bestprof")[0]
        except IndexError as e:
            raise IndexError(f"No .bestprofs found: {cfg['run_ops']['psr_dir']}")
        cfg["folds"][pointing]["init"][nbins] = bestprof_info(bestprof_name)

    # Search for pointings with positive classifications (>=3 out of 5 is positive classification)
    positive_pointings = [pointing for pointing in cfg["folds"].keys() if cfg["folds"][pointing]["classifier"]>=3]

    # Throw exception if there aren't any positive detections
    if len(positive_pointings) == 0:
        raise NoUsableFolds(f"No positive classifications found in pulsar directory {cfg['run_ops']['psr_dir']}")
    else:
        best_chi = 0
        for pointing in positive_pointings:
            nbins = list(cfg["folds"][pointing]["init"].keys())[0]
            if cfg["folds"][pointing]["init"][nbins]["chi"] > best_chi:
                cfg["source"]["my_pointing"] = pointing
        logger.info(f"Best pointing found with chi value of {best_chi}: {cfg['source']['my_pointing']}")


def populate_post_folds(cfg):
    """Fills the cfg with info on all of the post folds"""
    my_pointing = cfg["source"]["my_pointing"]
    for bins in cfg["folds"][my_pointing]["post"].keys():
        try:
            bestprof_name = glob_pfds(cfg, my_pointing, bins, pfd_type="pfd.bestprof")[0]
        except IndexError as _:
            raise IndexError(f"No .bestprofs found: {cfg['run_ops']['psr_dir']}")
        cfg["folds"][my_pointing]["post"][bins] = bestprof_info(bestprof_name)


def best_post_fold(cfg):
    """Finds the best fold to use in the cfg"""
    pointing = cfg["source"]["my_pointing"]
    min_chi = cfg["run_ops"]["thresh_chi"]
    min_sn = cfg["run_ops"]["thresh_sn"]
    good_chi = cfg["run_ops"]["good_chi"]
    good_sn = cfg["run_ops"]["good_sn"]
    post_folds = [int(i) for i in cfg["folds"][pointing]["post"].keys()]
    post_folds = sorted(post_folds, reverse=True) # Sort from highest to lowest bins

    # "good" loop
    for bin_count in post_folds:
        info = cfg["folds"][pointing]["post"][str(bin_count)]
        if info["sn"] >= good_sn and info["chi"] >= good_chi:
            cfg["source"]["my_bins"] = bin_count
            break

    # "minimum requirements" loop
    if cfg["source"]["my_bins"] == None:
        logger.info("No folds meet 'good' criteria")
        for bin_count in post_folds:
            info = cfg["folds"][pointing]["post"][str(bin_count)]
            if info["sn"] >= min_sn and info["chi"] >= min_chi:
                cfg["source"]["my_bins"] = bin_count
                break

    if cfg["source"]["my_bins"] == None:
        # The classifier still gave a positive detection somewhere. Continue with lowest bin post fold
        logger.warn("No folds meet minimum requirements. Will continue with lowest bin fold")
        cfg["source"]["my_bins"] = post_folds[-1]
    else:
        logger.info(f"Continuing with bin count: {cfg['source']['my_bins']}")