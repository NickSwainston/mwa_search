#!/usr/bin/env python3
import logging

from dpp.helper_yaml import from_yaml, dump_to_yaml

logger = logging.getLogger(__name__)

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


def _populate_master(master, kwargs):
    # Populate with yamls
    for f in kwargs["yamls"]:
        pipe = from_yaml(f)
        pointing = pipe["run_ops"]["pointing"]
        master[pointing] = {}
        master[pointing]["pipe"] = pipe
    # Populate with bestprof info
    bestprofs = [i for i in kwargs["pfds"] if "bestprof" in i]
    pointings = master_dict.keys()
    for p in pointings:
        for b in bestprof:
            if b.find(p)>=0:
                master_dict[p]["bestprof"] = {}
                master_dict[p]["bestprof"]["name"] = b
                info = bestprof_info(b)
                master_dict[p]["bestprof"]["info"] = info
                # Add bprof info to pipe
                master_dict[p]["pipe"]["folds"]["init"][info["nbins"]] = info
                break


def _eval_master(master):
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


def find_best_pointing_main(kwargs):
    """Decides the best folding solution from bestprofs"""
    master_dict = {}
    _populate_master(kwargs["yamls"], pfds, master_dict)
    best_pointing = _eval_master(master_dict)
    _remove_bad_pointings(best_pointing, kwargs)
    # Update yaml file
    dump_to_yaml(master_dict[best_pointing], label="find_best_pointing")
    
        