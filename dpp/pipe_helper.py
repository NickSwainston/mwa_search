#!/usr/bin/env python3
from vcstools import data_load
import psrqpy
from misc_helper import bin_sampling_limit, is_binary, required_folds
from pulsar_obs_helper import find_fold_times

logger = logging.getLogger(__name__)

def initiate_pipe(kwargs):
    """Adds all available keys to the pipe dictionary"""
    pipe = dict(kwargs)
    md = get_common_obs_metadata(pipe["obsid"])
    pipe["obs_ra"] = md[1]
    pipe["obs_dec"] = md[2]
    pipe["obs_dur"] = md[3]
    pipe["obs_freq"] = md[5]

    if pipe["cand"] == False:
        query = psrqpy.QueryATNF(psrs=pipe["pulsar"], loadfromdb=data_load.ATNF_LOC).pandas
        pipe["source_ra"] = query["RAJ"][0]
        pipe["source_dec"] = query["DECJ"][0]
        pipe["ATNF_RM"] = query["RM"][0]
        pipe["ATNF_RM_e"] = query["RM_ERR"][0]
        pipe["binary"] = is_binary(pulsar, query=query)
        pipe["sampling_limit"] = bin_sampling_limit(pulsar, query=query)
        pipe["source_enter_frac"], pipe["source_exit_frac"], pipe["power"] = find_fold_times(pipe["pulsar"], pipe["obsid"], pipe["obs_beg"], pipe["end"])
        pipe["init_folds"], pipe["required_folds"] = required_folds(pipe["pulsar"], query=query)
        
    pipe["init_folds_complete"] = False
    pipe["required_folds_complete"] = False
    pipe["bf_complete"] = False
    pipe["polarimetry_complete"] = False

    return pipe
