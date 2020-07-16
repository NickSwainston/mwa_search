#!/usr/bin/env python3
from vcstools import data_load
import psrqpy
from misc_helper import bin_sampling_limit, is_binary, required_bin_folds
from pulsar_obs_helper import find_fold_times

logger = logging.getLogger(__name__)


def initiate_pipe(kwargs):
    """Adds all available keys to the pipe dictionary and figures out some useful constants"""
    pipe = {"obs": {}, "source": {},
            "completed": {}, "folds": {}, "run_ops": {}}

    pipe["run_ops"]["dirs"] = kwargs["run_dirs"]
    pipe["run_ops"]["my_dir"] = None
    pipe["run_ops"]["loglvl"] = kwargs["loglvl"]
    pipe["run_ops"]["mwa_search"] = kwargs["mwa_search"]
    pipe["run_ops"]["vcstools"] = kwargs["vcstools"]
    pipe["run_ops"]["mask"] = kwargs["mask"]
    pipe["run_ops"]["thresh_chi"] = 4.0
    pipe["run_ops"]["thresh_sn"] = 8.0

    md = get_common_obs_metadata(pipe["obsid"])
    pipe["obs"]["ra"] = md[1]
    pipe["obs"]["dec"] = md[2]
    pipe["obs"]["dur"] = md[3]
    pipe["obs"]["freq"] = md[5]
    pipe["obs"]["id"] = kwargs["obsid"]
    pipe["obs"]["cal"] = kwargs["cal_id"]
    pipe["obs"]["beg"] = kwargs["beg"]
    pipe["obs"]["end"] = kwargs["end"]

    pipe["source"]["cand"] = kwargs["cand"]

    if pipe["cand"] == False:
        query = psrqpy.QueryATNF(
            psrs=pipe["pulsar"], loadfromdb=data_load.ATNF_LOC).pandas
        pipe["source"]["name"] = kwargs["pulsar"]
        pipe["source"]["ra"] = query["RAJ"][0]
        pipe["source"]["dec"] = query["DECJ"][0]
        pipe["source"]["ATNF_RM"] = query["RM"][0]
        pipe["source"]["ATNF_RM_e"] = query["RM_ERR"][0]
        pipe["source"]["synth_RM"] = None
        pipe["source"]["synth_RM_e"] = None
        pipe["source"]["fit_RM"] = None
        pipe["source"]["fit_RM_e"] = None
        pipe["source"]["ATNF_P"] = query["P0"][0]
        pipe["source"]["ATNF_DM"] = query["DM"][0]
        pipe["source"]["my_DM"] = None
        pipe["source"]["my_P"] = None
        pipe["source"]["binary"] = is_binary(pulsar, query=query)
        pipe["source"]["sampling_limit"] = bin_sampling_limit(
            pulsar, query=query)
        pipe["source"]["enter_frac"], pipe["source"]["exit_frac"], pipe["source"]["power"] = find_fold_times(
            pipe["pulsar"], pipe["obsid"], pipe["obs_beg"], pipe["obs_end"])
        init, post = = required_bin_folds(pipe["pulsar"], query=query)
        pipe["folds"] = {}
        for _, i in enumerate(init):
            pipe["folds"]["init"][str(i)] = {}
        for _, i in enumerate(post):
            pipe["folds"]["post"][str(i)] = {}

    pipe["completed"] = {}
    pipe["completed"]["init_folds"] = False
    pipe["completed"]["post_folds"] = False
    pipe["compelted"]["submit_move"] = False
    pipe["completed"]["bf"] = False
    pipe["completed"]["polarimetry"] = False

    return pipe
