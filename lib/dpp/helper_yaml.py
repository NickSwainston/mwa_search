import psrqpy
import logging
import os
import yaml
import subprocess

from dpp.helper_source_info import bin_sampling_limit, is_binary, required_bin_folds
from dpp.helper_obs_info import find_fold_times
from mwa_metadb_utils import get_common_obs_metadata
from vcstools import data_load

logger = logging.getLogger(__name__)


def initiate_pipe(kwargs, psr, pointing, metadata=None, full_meta=None, query=None, enter=None, exit=None, power=None):
    """Adds all available keys to the pipe dictionary and figures out some useful constants"""
    pipe = {"obs": {}, "source": {},
            "completed": {}, "folds": {}, "run_ops": {}, "pol": {}}

    pipe["run_ops"]["loglvl"] = kwargs["loglvl"]
    pipe["run_ops"]["mwa_search"] = kwargs["mwa_search"]
    pipe["run_ops"]["vcstools"] = kwargs["vcstools"]
    pipe["run_ops"]["thresh_chi"] = 3.5
    pipe["run_ops"]["thresh_sn"] = 8.0
    pipe["run_ops"]["good_chi"] = 4.0
    pipe["run_ops"]["good_sn"] = 20.0
    pipe["run_ops"]["vdif"] = None
    pipe["run_ops"]["mask"] = None
    pipe["run_ops"]["pointing"] = pointing
    check_run_ops(pipe)

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
    check_obs_inputs(pipe)

    pipe["source"]["cand"] = kwargs["cand"]
    if pipe["source"]["cand"] == False:
        pipe["source"]["name"] = psr
        if query is None:
            query = psrqpy.QueryATNF(
                psrs=pipe["source"]["name"], loadfromdb=data_load.ATNF_LOC).pandas
        pipe["source"]["sampling_limit"] = int(bin_sampling_limit(
            pipe["source"]["name"], query=query))
        pipe["source"]["ATNF"] = dict(query)
        pipe["source"]["ATNF_P"] = query["P0"][0]
        pipe["source"]["ATNF_DM"] = query["DM"][0]
        pipe["source"]["RM_type"] = None
        pipe["source"]["synth_RM"] = None
        pipe["source"]["synth_RM_e"] = None
        pipe["source"]["my_DM"] = None
        pipe["source"]["my_P"] = None
        pipe["source"]["my_bins"] = None
        pipe["source"]["edited_eph"] = None
        pipe["source"]["edited_eph_name"] = None
        pipe["source"]["enter_frac"] = enter
        pipe["source"]["exit_frac"] = exit
        pipe["source"]["power"] = power
        if None in (pipe["source"]["enter_frac"], pipe["source"]["exit_frac"], pipe["source"]["power"]):
            pipe["source"]["enter_frac"], pipe["source"]["exit_frac"], pipe["source"]["power"] = find_fold_times(
                pipe["source"]["name"], pipe["obs"]["id"], pipe["obs"]["beg"], pipe["obs"]["end"], metadata=metadata, full_meta=full_meta)
        init, post = required_bin_folds(pipe["source"]["name"], query=query)
        pipe["folds"] = {"init":{}, "post":{}, "best":{}}
        for _, i in enumerate(init):
            pipe["folds"]["init"][str(i)] = {}
        for _, i in enumerate(post):
            pipe["folds"]["post"][str(i)] = {}
        pipe["source"]["binary"] = is_binary(pipe["source"]["name"], query=query)
        check_source_inputs(pipe)
        pipe["source"]["seek"] = pipe["source"]["enter_frac"] * (pipe["obs"]["end"] - pipe["obs"]["beg"])
        pipe["source"]["total"] = (pipe["source"]["exit_frac"] - pipe["source"]["enter_frac"]) * (pipe["obs"]["end"] - pipe["obs"]["beg"])
        pipe["run_ops"]["file_precursor"] = f"{pipe['obs']['id']}_{pipe['run_ops']['pointing']}_{pipe['source']['name']}"
        if pipe["source"]["binary"]:
            pipe["source"]["edited_eph_name"] = f"{pipe['run_ops']['file_precursor']}.eph"
            pipe["source"]["edited_eph"] = create_edited_eph(pipe["source"]["name"], pipe["source"]["edited_eph_name"])

    pipe["pol"]["alpha"] = None
    pipe["pol"]["beta"] = None
    pipe["pol"]["l0"] = None
    pipe["pol"]["pa0"] = None
    pipe["pol"]["chi"] = None

    pipe["completed"] = {}
    pipe["completed"]["init_folds"] = False
    pipe["completed"]["post_folds"] = False
    pipe["completed"]["upload_and_move"] = False
    pipe["completed"]["bf"] = False
    pipe["completed"]["polarimetry_1"] = False
    pipe["completed"]["polarimetry_2"] = False
    pipe["completed"]["polarimetry_3"] = False
    pipe["completed"]["polarimetry_4"] = False
    pipe["completed"]["polarimetry_5"] = False
    pipe["completed"]["polarimetry_6"] = False


    return pipe

def check_run_ops(pipe):
    """Checks that the 'run_ops' information given to the pipe is suitable"""
    #pointing
    if not isinstance(pipe["run_ops"]["pointing"], str):
        raise TypeError(f"Pointing not valid: {pipe['run_ops']['pointing']}")


def check_obs_inputs(pipe):
    """Checks that the 'obs' information given to the pipe is suitable"""
    # Obsid
    if not isinstance(pipe["obs"]["id"], int):
        try:
            pipe["obs"]["id"] = int(pipe["obs"]["id"])
            logger.warn("Obsid had to be converted to int. This may be evidence of a bug")
        except (ValueError, TypeError) as e:
            e.message = f"Invalid Observation ID: {pipe['obs']['id']}. Cannot be converted to int"
            raise
    # Beg and end
    if not isinstance(pipe["obs"]["beg"], int):
        try:
            pipe["obs"]["beg"] = int(pipe["obs"]["beg"])
            logger.warn("Begin time had to be converted to int. This may be evidence of a bug")
        except (ValueError, TypeError) as e:
            e.message = f"Invalid begin time: {pipe['obs']['beg']}. Cannot be converted to int"
            raise
    if not isinstance(pipe["obs"]["end"], int):
        try:
            pipe["obs"]["end"] = int(pipe["obs"]["end"])
            logger.warn("End time had to be converted to int. This may be evidence of a bug")
        except (ValueError, TypeError) as e:
            e.message = f"Invalid end time: {pipe['obs']['end']}. Cannot be converted to int"
            raise
    if pipe["obs"]["beg"]>pipe["obs"]["end"]:
        raise ValueError(f"Begining time {pipe['obs']['beg']} greater than end time {pipe['obs']['end']}")


def check_source_inputs(pipe):
    """Checks if the 'source' info given to the pipe is suiutable"""
    #ANTF stuff
    if not isinstance(pipe["source"]["ATNF_P"], float):
        try:
            pipe["source"]["ATNF_P"] = float(pipe["source"]["ATNF_P"])
        except (ValueError, TypeError) as e:
            e.message = f"Invalid ATNF period: {pipe['source']['ATNF_P']}. Cannot be converted to float"
            raise
    if not isinstance(pipe["source"]["ATNF_DM"], float):
        try:
            pipe["source"]["ATNF_DM"] = float(pipe["source"]["ATNF_DM"])
        except (ValueError, TypeError) as e:
            e.message = f"Invalid ATNF dispersion measure: {pipe['source']['ATNF_DM']}. Cannot be converted to float"
            raise
    #Enter/Exit fractions
    if not isinstance(pipe["source"]["enter_frac"], float):
        try:
            pipe["source"]["enter_frac"] = float(pipe["source"]["enter_frac"])
        except (ValueError, TypeError) as e:
            e.message = f"Invalid beam enter fraction: {pipe['source']['enter_frac']}. Cannot be converted to float"
            raise
    if not isinstance(pipe["source"]["exit_frac"], float):
        try:
            pipe["source"]["exit_frac"] = float(pipe["source"]["exit_frac"])
        except (ValueError, TypeError) as e:
            e.message = f"Invalid beam exit fraction: {pipe['source']['exit_frac']}. Cannot be converted to float"
            raise
    if pipe["source"]["enter_frac"] > 1 or pipe["source"]["exit_frac"] < 0:
        msg = f"""Enter/Exit fractions unsuitable
                  Enter: {pipe['source']['enter_frac']}
                  Exit: {pipe['source']['exit_frac']}"""
        raise ValueError(msg)
    #Power
    if not isinstance(pipe["source"]["power"], float):
        try:
            pipe["source"]["power"] = float(pipe["source"]["power"])
        except (ValueError, TypeError) as e:
            e.message = f"Invalid maximum power: {pipe['source']['power']}. Cannot be converted to float"
            raise
    #Sampling limit
    if not isinstance(pipe["source"]["sampling_limit"], int):
        try:
            pipe["source"]["power"] = float(pipe["source"]["power"])
        except (ValueError, TypeError) as e:
            e.message = f"Invalid sampling limit: {pipe['source']['sampling_limit']}. Cannot be converted to float"
            raise
    #Binary
    if not isinstance(pipe["source"]["binary"], bool):
        raise ValueError(f"Invalid binary condition: {pipe['source']['binary']}")


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


def create_yaml_main(kwargs):
    """uses kwargs from make_pulsar_yaml.py"""
    metadata, full_meta = get_common_obs_metadata(kwargs["obsid"], return_all=True)
    query = psrqpy.QueryATNF(loadfromdb=data_load.ATNF_LOC).pandas
    pulsars_pointings_dict = {}
    for psrlist, pointing in zip(kwargs["psr"], kwargs["pointing"]):
        for psr in psrlist.split(":"):
            if psr not in pulsars_pointings_dict.keys():
                pulsars_pointings_dict[psr] = []
            pulsars_pointings_dict[psr].append(pointing)
    for psr in pulsars_pointings_dict.keys():
        logger.info("Processing yaml for PSR: {}".format(psr))
        enter_frac, exit_frac, power = find_fold_times(
                psr, kwargs["obsid"], kwargs["obs_beg"], kwargs["obs_end"], metadata=metadata, full_meta=full_meta)
        for pointing in pulsars_pointings_dict[psr]:
            try:
                pipe = initiate_pipe(kwargs, psr, pointing, metadata=metadata, full_meta=full_meta, query=query[query['PSRJ'] == psr].reset_index(),
                        enter=enter_frac, exit=exit_frac, power=power)
            except (ValueError, TypeError, OSError) as e:
                msg = f"""Exception encountered for pulsar {psr} and pointing {pointing}
                          Error: {e} """
                logger.warn(msg)
                continue
            dump_to_yaml(pipe, label=kwargs["label"])