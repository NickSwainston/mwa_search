import psrqpy
import logging
import os
import yaml

from dpp.misc_helper import bin_sampling_limit, is_binary, required_bin_folds
from dpp.pulsar_obs_helper import find_fold_times
from mwa_metadb_utils import get_common_obs_metadata
from vcstools import data_load

logger = logging.getLogger(__name__)


def initiate_pipe(kwargs, psr, metadata=None, full_meta=None, query=None):
    """Adds all available keys to the pipe dictionary and figures out some useful constants"""
    pipe = {"obs": {}, "source": {},
            "completed": {}, "folds": {}, "run_ops": {}, "pol": {}}

    pipe["run_ops"]["dirs"] = kwargs["run_dir"]
    pipe["run_ops"]["loglvl"] = kwargs["loglvl"]
    pipe["run_ops"]["mwa_search"] = kwargs["mwa_search"]
    pipe["run_ops"]["vcstools"] = kwargs["vcstools"]
    pipe["run_ops"]["thresh_chi"] = 4.0
    pipe["run_ops"]["thresh_sn"] = 8.0
    pipe["run_ops"]["vdif"] = None
    pipe["run_ops"]["mask"] = None

    if metadata is None:
        metadata = get_common_obs_metadata(kwargs["obsid"])
    pipe["obs"]["ra"] = metadata[1]
    pipe["obs"]["dec"] = metadata[2]
    pipe["obs"]["dur"] = metadata[3]
    pipe["obs"]["freq"] = metadata[5]
    pipe["obs"]["id"] = kwargs["obsid"]
    pipe["obs"]["cal"] = kwargs["cal_id"]
    pipe["obs"]["beg"] = kwargs["obs_beg"]
    pipe["obs"]["end"] = kwargs["obs_end"]

    pipe["source"]["cand"] = kwargs["cand"]

    if pipe["source"]["cand"] == False:
        pipe["source"]["name"] = psr
        if query is None:
            query = psrqpy.QueryATNF(
                psrs=pipe["source"]["name"], loadfromdb=data_load.ATNF_LOC).pandas
        pipe["source"]["ATNF"] = dict(query)
        pipe["source"]["ATNF_P"] = query["P0"][0]
        pipe["source"]["ATNF_DM"] = query["DM"][0]
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
        pipe["source"]["sampling_limit"] = int(bin_sampling_limit(
            pipe["source"]["name"], query=query))
        pipe["source"]["enter_frac"], pipe["source"]["exit_frac"], pipe["source"]["power"] = find_fold_times(
            pipe["source"]["name"], pipe["obs"]["id"], pipe["obs"]["beg"], pipe["obs"]["end"], metadata=metadata, full_meta=full_meta)
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
        pipe["source"]["binary"] = is_binary(pipe["source"]["name"], query=query)
        pipe["source"]["edited_eph"] = None
        pipe["source"]["edited_eph_name"] = None
        #create an edited epehemris for binary folding if necessary
        if pipe["source"]["binary"]:
            pipe["source"]["edited_eph_name"] = f"{pipe['run_ops']['file_precursor']}.eph"
            pipe["source"]["edited_eph"] = create_edited_eph(pipe["source"]["name"], pipe["source"]["edited_eph_name"])


    pipe["pol"]["archive1"] = None
    pipe["pol"]["archive2"] = None
    pipe["pol"]["rmfit"] = None
    pipe["pol"]["rmsynth"] = None
    pipe["pol"]["rvmfit"] = None
    pipe["pol"]["rvmres"] = kwargs["rvmres"]

    pipe["completed"] = {}
    pipe["completed"]["init_folds"] = False
    pipe["completed"]["post_folds"] = False
    pipe["completed"]["upload_and_move"] = False
    pipe["completed"]["bf"] = False
    pipe["completed"]["polarimetry"] = False
    pipe["completed"]["init_dspsr"] = False

    return pipe


def create_edited_eph(pulsar_name, eph_name):
    """Created a string version of 'psrcat -e' and removes the last line"""
    eph = subprocess.check_output(["psrcat", "-e", pulsar_name])
    eph = eph.decode("utf-8")
    eph = "\n".join(tuple(eph.split("\n")[:-2]))
    return eph


def from_yaml(filepath):
    with open(filepath) as f:
        my_dict = yaml.load(f, Loader=yaml.Loader)
    return my_dict


def dump_to_yaml(pipe, label=""):
    name = f"{pipe['obs']['id']}_{pipe['run_ops']['pointing']}_{pipe['source']['name']}"
    if label:
        name += f"_{label}"
    name += ".yaml"
    with open(name, 'w') as f:
        yaml.dump(pipe, f, default_flow_style=False)