#!/usr/bin/env python3
#from vcstools import data_load
import psrqpy
import logging
import os
import glob

from misc_helper import bin_sampling_limit, is_binary, required_bin_folds
from pulsar_obs_helper import find_fold_times
from mwa_metadb_utils import get_common_obs_metadata


logger = logging.getLogger(__name__)


try:  # get ATNF db location
    ATNF_LOC = os.environ['PSRCAT_FILE']
except:
    logger.warn(
        "ATNF database could not be loaded on disk. This may lead to a connection failure")
    ATNF_LOC = None


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
    pipe["run_ops"]["vdif"] = None

    md = get_common_obs_metadata(kwargs["obsid"])
    pipe["obs"]["ra"] = md[1]
    pipe["obs"]["dec"] = md[2]
    pipe["obs"]["dur"] = md[3]
    pipe["obs"]["freq"] = md[5]
    pipe["obs"]["id"] = kwargs["obsid"]
    pipe["obs"]["cal"] = kwargs["cal_id"]
    pipe["obs"]["beg"] = kwargs["obs_beg"]
    pipe["obs"]["end"] = kwargs["obs_end"]

    pipe["source"]["cand"] = kwargs["cand"]

    if pipe["source"]["cand"] == False:
        pipe["source"]["name"] = kwargs["pulsar"]
        query = psrqpy.QueryATNF(
            psrs=pipe["source"]["name"], loadfromdb=ATNF_LOC).pandas
        pipe["source"]["ATNF"] = dict(query)
        pipe["source"]["RM_type"] = None
        pipe["source"]["synth_RM"] = None
        pipe["source"]["synth_RM_e"] = None
        pipe["source"]["fit_RM"] = None
        pipe["source"]["fit_RM_e"] = None
        pipe["source"]["my_RM"] = None
        pipe["source"]["my_RM_e"] = None
        pipe["source"]["my_DM"] = None
        pipe["source"]["my_P"] = None
        pipe["source"]["my_bins"] = None
        pipe["source"]["binary"] = is_binary(pipe["source"]["name"], query=query)
        pipe["source"]["sampling_limit"] = int(bin_sampling_limit(
            pipe["source"]["name"], query=query))
        pipe["source"]["enter_frac"], pipe["source"]["exit_frac"], pipe["source"]["power"] = find_fold_times(
            pipe["source"]["name"], pipe["obs"]["id"], pipe["obs"]["beg"], pipe["obs"]["end"])
        pipe["source"]["enter_frac"] = float(pipe["source"]["enter_frac"])
        pipe["source"]["exit_frac"] = float(pipe["source"]["exit_frac"])
        pipe["source"]["seek"] = pipe["source"]["enter_frac"] * (pipe["obs"]["end"] - pipe["obs"]["beg"])
        pipe["source"]["total"] = (pipe["source"]["exit_frac"] - pipe["source"]["enter_frac"]) * (pipe["obs"]["end"] - pipe["obs"]["beg"])
        pipe["source"]["power"] = float(pipe["source"]["power"])
        init, post = required_bin_folds(pipe["source"]["name"], query=query)
        pipe["folds"] = {"init":{}, "post":{}}
        for _, i in enumerate(init):
            pipe["folds"]["init"][str(i)] = {}
        for _, i in enumerate(post):
            pipe["folds"]["post"][str(i)] = {}

    pipe["pol"]["archive1"] = None
    pipe["pol"]["archive2"] = None
    pipe["pol"]["rmfit"] = None
    pipe["pol"]["rmsynth"] = None
    pipe["pol"]["rvmfit"] = None

    pipe["completed"] = {}
    pipe["completed"]["init_folds"] = False
    pipe["completed"]["post_folds"] = False
    pipe["completed"]["upload_and_move"] = False
    pipe["completed"]["bf"] = False
    pipe["completed"]["polarimetry"] = False
    pipe["completed"]["init_dspsr"] = False

    return pipe
