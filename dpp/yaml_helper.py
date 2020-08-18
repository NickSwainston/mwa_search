#!/usr/bin/env python3
import psrqpy
import logging
import os
import glob
import yaml
import argparse

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


def initiate_pipe(kwargs, psr, pointing, metadata=None, query=None):
    """Adds all available keys to the pipe dictionary and figures out some useful constants"""
    pipe = {"obs": {}, "source": {},
            "completed": {}, "folds": {}, "run_ops": {}, "pol": {}}

    pipe["run_ops"]["dirs"] = kwargs["run_dir"]
    pipe["run_ops"]["pointing"] = pointing
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
                psrs=pipe["source"]["name"], loadfromdb=ATNF_LOC).pandas
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
    pipe["pol"]["rvmres"] = kwargs["rvmres"]

    pipe["completed"] = {}
    pipe["completed"]["init_folds"] = False
    pipe["completed"]["post_folds"] = False
    pipe["completed"]["upload_and_move"] = False
    pipe["completed"]["bf"] = False
    pipe["completed"]["polarimetry"] = False
    pipe["completed"]["init_dspsr"] = False

    return pipe


def from_yaml(filepath):
    with open(filepath) as f:
        my_dict = yaml.load(f, Loader=yaml.Loader)
    return my_dict


def dump_to_yaml(pipe, filepath=None):
    if filepath == None:
        filepath = f"{pipe['run_ops']['pointing']}_{pipe['obs']['id']}_{pipe['source']['name']}.yaml"
    with open(filepath, 'w') as f:
        yaml.dump(pipe, f, default_flow_style=False)


def main(kwargs):
    metadata = get_common_obs_metadata(kwargs["obsid"])
    query = psrqpy.QueryATNF(loadfromdb=ATNF_LOC).pandas
    for psr, pointing in zip(kwargs["psrs"], kwargs["pointings"]):
        pipe = initiate_pipe(kwargs, psr, pointing, metadata=metadata, query=query[query['PSRJ'] == psr].reset_index())
        dump_to_yaml(pipe)


if __name__ == '__main__':
    loglevels = dict(DEBUG=logging.DEBUG,
                     INFO=logging.INFO,
                     WARNING=logging.WARNING,
                     ERROR=logging.ERROR)
    parser = argparse.ArgumentParser(description="""Initialises a .yaml file with all pulsar info required for a DPP run""",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    obsop = parser.add_argument_group("Observation Options")
    obsop.add_argument("-d", "--run_dir", nargs='+',
                       type=str, help="The location of the pointing directory/s")
    obsop.add_argument("-o", "--obsid", type=str,
                       help="The obs ID of the data")
    obsop.add_argument("-O", "--cal_id", type=str,
                       help="The ID of the calibrator used to calibrate the data")
    obsop.add_argument("-p", "--psrs", nargs="+", type=str,
                       help="The J name of the pulsar(s). e.g. J2241-5236")
    obsop.add_argument("--obs_beg", type=int,
                       help="The beginning of the observation")
    obsop.add_argument("--obs_end", type=int, help="The end of the observation")
    obsop.add_argument("--pointings", type=str, nargs="+", help="The pointing(s) location of the source")

    foldop = parser.add_argument_group("Folding/processing Options")
    foldop.add_argument("--sn_min_thresh", type=float, default=8.0, help="The presto sigma value\
                             above which is deemed a detection.")
    foldop.add_argument("--sn_max_thresh", type=float, default=20.0, help="The presto sigma value\
                             above which is deemed a GOOD detection.")
    foldop.add_argument("--chi_thresh", type=float, default=4.0, help="The presto 'chi' value\
                             above which is deemed a detection.")
    foldop.add_argument("--rvmres", type=int, default=90,
                        help="The number of degree samples to try for alpha and beta.")
    foldop.add_argument("--mask", type=str,
                        help="The pathname of the mask to use for folding")

    otherop = parser.add_argument_group("Other Options")
    otherop.add_argument("--no_ephem", action="store_true",
                         help="Use this to override the use of an ephemeris for foldign the pulsar")
    otherop.add_argument("-L", "--loglvl", type=str, default="INFO",
                         help="Logger verbosity level", choices=loglevels.keys())
    otherop.add_argument("--mwa_search", type=str, default="master",
                         help="The version of mwa_search to use")
    otherop.add_argument("--vcstools", type=str, default="master",
                         help="The version of vcs_tools to use")
    otherop.add_argument("--cand", action="store_true",
                         help="use this tag if this is not a kown pulsar")
    args = parser.parse_args()
    logger = logging.getLogger()
    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    formatter = logging.Formatter(
        '%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    kwargs = vars(args)
    main(kwargs)