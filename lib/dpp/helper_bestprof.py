#!/usr/bin/env python3
import logging

from dpp.helper_yaml import from_yaml, dump_to_yaml

logger = logging.getLogger(__name__)


class NoUsableFolds(Exception):
    """Raise when no usable folds are found in a pipe"""
    pass


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
    info_dict["obsid"] = int(lines[0].split()[4].split("_")[0])
    info_dict["pulsar"] = lines[1].split()[3].split("_")[1]
    info_dict["nbins"] = int(lines[9].split()[4])
    info_dict["chi"] = float(lines[12].split()[4])
    info_dict["sn"] = float(lines[13].split()[4][2:])
    info_dict["dm"] = float(lines[14].split()[4])
    info_dict["period"] = float(lines[15].split()[4])/1e3 #in seconds
    info_dict["period_error"] = float(lines[15].split()[6])/1e3
    f.close()
    return info_dict


def _populate_master_pointings(master, kwargs):
    # Populate with yamls
    for f in kwargs["yamls"]:
        pipe = from_yaml(f)
        pointing = pipe["run_ops"]["pointing"]
        master[pointing] = {}
        master[pointing]["pipe"] = pipe
    # Populate with bestprof info
    bestprofs = [i for i in kwargs["pfds"] if ".pfd.bestprof" in i]
    pointings = master_dict.keys()
    for p in pointings:
        for b in bestprofs:
            if b.find(p)>=0:
                master_dict[p]["bestprof"] = {}
                master_dict[p]["bestprof"]["name"] = b
                info = bestprof_info(b)
                master_dict[p]["bestprof"]["info"] = info
                # Add bprof info to pipe
                master_dict[p]["pipe"]["folds"]["init"][str(info["nbins"])] = info
                break


def _eval_master_init(master):
    best_eval = 0
    for p in master_dict:
        this_eval = master_dict[p]["bestprof"]["info"]["sn"] * master_dict[p]["bestprof"]["info"]["chi"]
        if this_eval > best_eval:
            best_eval = this_eval
            best_pointing = p
    return best_pointing


def _remove_bad_pointings(best_pointing, kwargs):
    import os
    for f in kwargs["yamls"] + kwargs["pfds"]:
        if f.find(best_pointing) == -1:
            os.remove(f)


def _populate_master_post_folds(master, kwargs):
    # Populate with yamls
    for f in kwargs["yamls"]:
        pipe = from_yaml(f)
        pointing = pipe["run_ops"]["pointing"]
        master[pointing] = {}
        master[pointing]["pipe"] = pipe
    # Populate with bestprof info
    bestprofs = [i for i in kwargs["pfds"] if ".pfd.bestprof" in i]
    pointings = master_dict.keys()
    for p in pointings:
        for b in bestprofs:
            if b.find(p)>=0:
                master_dict[p]["bestprof"] = {}
                info = bestprof_info(b)
                master_dict[p]["bestprof"][str(info["nbins"])] = info
                # Add bprof info to pipe
                master_dict[p]["pipe"]["folds"]["post"][str(info["nbins"])] = info


def _eval_post_folds(master):
    """Finds the bin count to use for the rest of the dpp pipeline for each pointing"""
    for p in master.keys():
        try:
            best = best_post_fold(master[p]["pipe"][:])
        except NoUsableFolds as e:
            logger.warn(f"""Exception encountered: {e.message}
                        Will use initial fold for this pointing""")
            best = [int(i) for i in master[p]["pipe"]["folds"]["init"].keys()]
            best = max(best)
        master[p]["pipe"]["folds"]["best"] = best


def best_post_fold(pipe):
    """Finds the best fold to use in the pipe and returns the bin count"""
    min_chi = pipe["run_ops"]["thresh_chi"]
    min_sn = pipe["run_ops"]["thresh_sn"]
    good_chi = pipe["run_ops"]["good_chi"]
    good_sn = pipe["run_ops"]["good_sn"]
    post_folds = [int(i) for i in pipe["folds"]["post"].keys()]
    post_folds = sorted(post_folds, reverse=True)
    # "good" loop
    best = None
    for bin_count in post_folds:
        info = pipe["folds"]["post"][str(bin_count)]
        if info["sn"] >= good_sn and info["chi"] >= good_chi:
            best = bin_count
            break
    # "minimum requirements" loop
    if best == None:
        for bin_count in post_folds:
            info = pipe["folds"]["post"][str(bin_count)]
            if info["sn"] >= min_sn and info["chi"] >= min_chi:
                best = bin_count
                break
    if best == None:
        raise NoUsableFolds(f"""No folds meeting the minumum requirements found for pointing {pipe['run_ops']['pointing']}
                            Minimum requirements:
                            S/N: {pipe['run_ops']['sn_thresh']}
                            Chi: {pipe['run_ops']['chi_thresh']}""")
    return best


def _dump_master_pointings(master, label=""):
    """Dumps all of the pipes in master to yaml files"""
    for p in master.keys():
        dump_to_yaml(master[p]["pipe"], label=label)


def find_best_pointing_main(kwargs):
    """Decides the best folding solution from bestprofs"""
    master_dict = {}
    _populate_master(kwargs["yamls"], pfds, master_dict)
    best_pointing = _eval_master_init(master_dict)
    _remove_bad_pointings(best_pointing, kwargs)
    # Update yaml file
    dump_to_yaml(master_dict[best_pointing], label=kwargs["label"])


def post_fold_filter_main(kwargs):
    """Decides the best post-fold detection from bestprofs for each pointing supplied"""
    # Master dict heirarchy:
    #   Pointing
    #       bestprof
    #           info
    #       pipe
    #           *pipe_info
    master_dict = {}
    _populate_master_post_folds(master_dict, kwargs)
    _eval_post_folds(master_dict)
    _dump_master_pointings(master_dict, label=kwargs["label"])