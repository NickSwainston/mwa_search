#!/usr/bin/env python3
import psrqpy
import logging
from os.path import join
import yaml
import subprocess

from dpp.helper_source_info import bin_sampling_limit, is_binary, required_bin_folds
from dpp.helper_obs_info import find_fold_times
from dpp.helper_files import file_precursor
from vcstools.metadb_utils import get_common_obs_metadata
from vcstools.progress_bar import progress_bar
from vcstools import data_load
from vcstools.config import load_config_file


comp_config = load_config_file()
logger = logging.getLogger(__name__)


def initiate_cfg(kwargs, psr, pointings, metadata=None, full_meta=None, query=None):
    """
    Adds all available keys to the cfg dictionary and figures out some useful constants
    Takes kwargs from observation_processing_pipeline
    """
    cfg = {"obs": {}, "source": {}, "completed": {}, "folds": {}, "run_ops": {}, "pol": {}}
    if metadata is None or full_meta is None:
        metadata, full_meta = get_common_obs_metadata(kwargs["obsid"], return_all=True)
    if query is None:
        query = psrqpy.QueryATNF(loadfromdb=data_load.ATNF_LOC).pandas

    query_index = list(query["JNAME"]).index(psr)

    cfg["obs"]["ra"] = metadata[1]
    cfg["obs"]["dec"] = metadata[2]
    cfg["obs"]["dur"] = metadata[3]
    cfg["obs"]["freq"] = metadata[5]
    cfg["obs"]["id"] = kwargs["obsid"]
    cfg["obs"]["cal"] = kwargs["calid"]
    cfg["obs"]["beg"] = kwargs["beg"]
    cfg["obs"]["end"] = kwargs["end"]

    cfg["run_ops"]["loglvl"] = kwargs["loglvl"]
    cfg["run_ops"]["mwa_search"] = kwargs["mwa_search"]
    cfg["run_ops"]["vcstools"] = kwargs["vcstools"]
    cfg["run_ops"]["label"] = kwargs["label"]
    cfg["run_ops"]["thresh_chi"] = 3.5
    cfg["run_ops"]["thresh_sn"] = 8.0
    cfg["run_ops"]["good_chi"] = 4.0
    cfg["run_ops"]["good_sn"] = 20.0
    cfg["run_ops"]["vdif"] = None
    cfg["run_ops"]["mask"] = None
    cfg["run_ops"]["file_precursor"] = file_precursor(kwargs, psr)
    cfg["run_ops"]["psr_dir"] = join(comp_config["base_data_dir"], str(cfg["obs"]["id"]), "dpp", cfg["run_ops"]["file_precursor"])
    cfg["run_ops"]["myname"] = join(cfg["run_ops"]["psr_dir"], f"{cfg['run_ops']['file_precursor']}_cfg.yaml")
    cfg["run_ops"]["logfile"] = join(cfg["run_ops"]["psr_dir"], f"{cfg['run_ops']['file_precursor']}_logs.txt")
    cfg["run_ops"]["batch_dir"] = join(comp_config['base_data_dir'], cfg["obs"]["id"], "batch")
    cfg["run_ops"]["classify_dir"] = join(cfg["run_ops"]["psr_dir"], "classifier_ppp")

    cfg["source"]["name"] = psr
    cfg["source"]["sampling_limit"] = int(bin_sampling_limit(cfg["source"]["name"], query=query))
    cfg["source"]["ATNF_P"] = float(query["P0"][query_index])
    cfg["source"]["ATNF_DM"] = float(query["DM"][query_index])
    cfg["source"]["synth_RM"] = None
    cfg["source"]["synth_RM_e"] = None
    cfg["source"]["my_DM"] = None
    cfg["source"]["my_P"] = None
    cfg["source"]["my_bins"] = None
    cfg["source"]["enter_frac"] = None
    cfg["source"]["exit_frac"] = None
    cfg["source"]["power"] = None
    cfg["source"]["my_pointing"] = None
    cfg["source"]["my_bins"] = None

    cfg["source"]["enter_frac"], cfg["source"]["exit_frac"], cfg["source"]["power"] = find_fold_times(
        cfg["source"]["name"], cfg["obs"]["id"], cfg["obs"]["beg"], cfg["obs"]["end"], metadata=metadata, full_meta=full_meta)
    try:
        cfg["source"]["enter_frac"] = float(cfg["source"]["enter_frac"])
        cfg["source"]["exit_frac"] = float(cfg["source"]["exit_frac"])
        cfg["source"]["power"] = float(cfg["source"]["power"])
    except ValueError as _: #If any of these are nones, the pulsar isn't in the beam for given beg, end
        raise ValueError(f"{cfg['source']['name']} not in beam for given times")
    cfg["source"]["binary"] = is_binary(cfg["source"]["name"], query=query)
    cfg["source"]["seek"] = cfg["source"]["enter_frac"] * (cfg["obs"]["end"] - cfg["obs"]["beg"])
    cfg["source"]["total"] = (cfg["source"]["exit_frac"] - cfg["source"]["enter_frac"]) * (cfg["obs"]["end"] - cfg["obs"]["beg"])
    if cfg["source"]["binary"]:
        cfg["source"]["edited_eph_name"] = join(cfg["run_ops"]["psr_dir"], f"{cfg['run_ops']['file_precursor']}.eph")
        cfg["source"]["edited_eph"] = create_edited_eph(cfg["source"]["name"])
        with open(cfg["source"]["edited_eph_name"], "w") as f:
            f.write(cfg["source"]["edited_eph"])
    else:
        cfg["source"]["edited_eph"] = None
        cfg["source"]["edited_eph_name"] = None
    init, post = required_bin_folds(cfg["source"]["name"], query=query)
    for pointing in pointings:
        cfg["folds"] = {pointing: {"init":{}, "post":{}}}
        cfg["folds"][pointing]["classifier"] = 0
        cfg["folds"][pointing]["dir"] = join(cfg["run_ops"]["psr_dir"], pointing)
        for _, i in enumerate(init):
            cfg["folds"][pointing]["init"][str(i)] = {}
        for _, i in enumerate(post):
            cfg["folds"][pointing]["post"][str(i)] = {}

    cfg["pol"]["alpha"] = None
    cfg["pol"]["beta"] = None
    cfg["pol"]["l0"] = None
    cfg["pol"]["pa0"] = None
    cfg["pol"]["chi"] = None

    cfg["completed"] = {}
    cfg["completed"]["init_folds"] = False
    cfg["completed"]["classify"] = False
    cfg["completed"]["post_folds"] = False
    cfg["completed"]["upload"] = False
    cfg["completed"]["polarimetry_1"] = False
    cfg["completed"]["polarimetry_2"] = False
    cfg["completed"]["polarimetry_3"] = False
    cfg["completed"]["polarimetry_4"] = False
    cfg["completed"]["polarimetry_5"] = False
    cfg["completed"]["polarimetry_6"] = False
    return cfg


def create_edited_eph(pulsar_name):
    """Created a string version of 'psrcat -e2' and removes the last line"""
    eph = subprocess.check_output(["psrcat", "-e2", pulsar_name])
    eph = eph.decode("utf-8")
    eph = "\n".join(tuple(eph.split("\n")[:-2]))
    return eph


def reset_cfg(cfg):
    """Sets all progress to incomplete to force a rerun"""
    cfg["completed"]["init_folds"] = False
    cfg["completed"]["classify"] = False
    cfg["completed"]["post_folds"] = False
    cfg["completed"]["upload"] = False
    cfg["completed"]["bf"] = False
    cfg["completed"]["polarimetry_1"] = False
    cfg["completed"]["polarimetry_2"] = False
    cfg["completed"]["polarimetry_3"] = False
    cfg["completed"]["polarimetry_4"] = False
    cfg["completed"]["polarimetry_5"] = False
    cfg["completed"]["polarimetry_6"] = False


def from_yaml(filepath):
    with open(filepath) as f:
        my_dict = yaml.load(f, Loader=yaml.Loader)
    return my_dict


def dump_to_yaml(cfg):
    with open(cfg["run_ops"]["myname"], 'w') as f:
        yaml.dump(cfg, f, default_flow_style=False)
    return cfg["run_ops"]["myname"]


def create_cfgs_main(kwargs, psrs_pointing_dict):
    """
    uses kwargs from observation_processing_pipeline.py
    """
    metadata, full_meta = get_common_obs_metadata(kwargs["obsid"], return_all=True)
    query = psrqpy.QueryATNF(loadfromdb=data_load.ATNF_LOC).pandas
    cfgs = []
    for psr in progress_bar(psrs_pointing_dict.keys(), "Initiating pulsar configs: "):
        logger.info(psr)
        pointing_list = psrs_pointing_dict[psr]
        try:
            cfgs.append(initiate_cfg(kwargs, psr, pointing_list, metadata=metadata, full_meta=full_meta, query=query))
        except ValueError as e:
            logger.info(e)
            continue
    return cfgs