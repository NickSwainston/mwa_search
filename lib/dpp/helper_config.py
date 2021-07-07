import psrqpy
import logging
from os.path import join
import yaml
import subprocess
import numpy as np

from dpp.helper_source_info import bin_sampling_limit, is_binary, required_bin_folds
from dpp.helper_obs_info import find_fold_times
from dpp.helper_files import file_precursor
from vcstools.metadb_utils import get_common_obs_metadata
from vcstools.progress_bar import progress_bar
from vcstools import data_load
from vcstools.config import load_config_file


comp_config = load_config_file()
logger = logging.getLogger(__name__)


def initiate_cfg(kwargs, psr, pointings, enter, leave, power, query=None, metadata=None):
    """
    Adds all available keys to the cfg dictionary and figures out some useful constants
    Takes kwargs from observation_processing_pipeline
    """
    cfg = {"obs": {}, "source": {}, "completed": {}, "folds": {}, "run_ops": {}, "pol": {}, "files":{}}
    if query is None:
        query = psrqpy.QueryATNF(loadfromdb=data_load.ATNF_LOC).pandas
    if metadata is None:
        metadata = get_common_obs_metadata(kwargs["obsid"])

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

    cfg["files"]["file_precursor"] = file_precursor(kwargs, psr)
    cfg["files"]["psr_dir"] = join(comp_config["base_data_dir"], str(cfg["obs"]["id"]), "dpp", cfg["files"]["file_precursor"])
    cfg["files"]["batch_dir"] = join(comp_config['base_data_dir'], cfg["obs"]["id"], "batch")
    cfg["files"]["classify_dir"] = join(cfg["files"]["psr_dir"], "classifier_ppp")
    cfg["files"]["my_name"] = join(cfg["files"]["psr_dir"], f"{cfg['files']['file_precursor']}_cfg.yaml")
    cfg["files"]["logfile"] = join(cfg["files"]["psr_dir"], f"{cfg['files']['file_precursor']}.log")
    cfg["files"]["archive"] = join(cfg["files"]["psr_dir"], f"{cfg['files']['file_precursor']}_archive.ar")
    cfg["files"]["archive_ascii"] = join(cfg["files"]["psr_dir"], f"{cfg['files']['file_precursor']}_archive.txt")
    cfg["files"]["gfit_plot"] = join(cfg["files"]["psr_dir"], f"{cfg['files']['file_precursor']}_gfit.png")
    cfg["files"]["converted_fits"] = join(cfg["files"]["psr_dir"], f"{cfg['files']['file_precursor']}_archive.fits")
    # debased fits file needs to be the same as archive except for the extension
    arch = cfg["files"]["archive"].split(".ar")[0]
    cfg["files"]["debased_fits"] = f"{arch}.debase.gg"
    # paswing file needs to be the same as debase except for the .gg extension
    debase = cfg["files"]["debased_fits"].split(".gg")[0]
    cfg["files"]["paswing"] = f"{debase}.paswing"
    cfg["files"]["chigrid_initial_ps"] = join(cfg["files"]["psr_dir"], f"{cfg['files']['file_precursor']}_chigrid_initial.ps")
    cfg["files"]["paswing_initial_ps"] = join(cfg["files"]["psr_dir"], f"{cfg['files']['file_precursor']}_paswing_initial.ps")
    cfg["files"]["RVM_fit_initial"] = join(cfg["files"]["psr_dir"], f"{cfg['files']['file_precursor']}_RVM_fit_initial.out")
    cfg["files"]["chigrid_final_ps"] = join(cfg["files"]["psr_dir"], f"{cfg['files']['file_precursor']}_chigrid_final.ps")
    cfg["files"]["paswing_final_ps"] = join(cfg["files"]["psr_dir"], f"{cfg['files']['file_precursor']}_paswing_final.ps")
    cfg["files"]["RVM_fit_final"] = join(cfg["files"]["psr_dir"], f"{cfg['files']['file_precursor']}_RVM_fit_final.out")
    cfg["files"]["ppol_profile_ps"] = join(cfg["files"]["psr_dir"], f"{cfg['files']['file_precursor']}_profile.ps")
    cfg["files"]["ppol_polar_profile_ps"] = join(cfg["files"]["psr_dir"], f"{cfg['files']['file_precursor']}_pol.ps")

    cfg["source"]["enter_frac"] = None
    cfg["source"]["exit_frac"] = None
    cfg["source"]["power"] = None
    cfg["source"]["name"] = psr
    try:
        cfg["source"]["enter_frac"] = float(enter)
        cfg["source"]["exit_frac"] = float(leave)
        cfg["source"]["power"] = float(power)
    except TypeError as _: #If any of these are nones, the pulsar isn't in the beam for given beg, end
        raise TypeError(f"{cfg['source']['name']} not in beam for given times")
    cfg["source"]["sampling_limit"] = int(bin_sampling_limit(cfg["source"]["name"], query=query))
    cfg["source"]["ATNF_P"] = float(query["P0"][query_index])
    cfg["source"]["ATNF_DM"] = float(query["DM"][query_index])
    tau_sc = query["DM"][query_index]
    cfg["source"]["scattering"] = float(tau_sc) if not bool(np.isnan(tau_sc)) else None
    cfg["source"]["my_DM"] = None
    cfg["source"]["my_P"] = None
    cfg["source"]["my_Pdot"] = None
    cfg["source"]["my_bins"] = None
    cfg["source"]["my_pointing"] = None
    cfg["source"]["my_component"] = None
    cfg["source"]["gfit"] = None
    cfg["source"]["binary"] = is_binary(cfg["source"]["name"], query=query)
    cfg["source"]["seek"] = cfg["source"]["enter_frac"] * (cfg["obs"]["end"] - cfg["obs"]["beg"])
    cfg["source"]["total"] = (cfg["source"]["exit_frac"] - cfg["source"]["enter_frac"]) * (cfg["obs"]["end"] - cfg["obs"]["beg"])
    if cfg["source"]["binary"]:
        cfg["source"]["edited_eph_name"] = join(cfg["files"]["psr_dir"], f"{cfg['files']['file_precursor']}.eph")
        cfg["source"]["edited_eph"] = create_edited_eph(cfg["source"]["name"])
    else:
        cfg["source"]["edited_eph"] = None
        cfg["source"]["edited_eph_name"] = None
    init, post = required_bin_folds(cfg["source"]["name"], query=query)
    for pointing in pointings:
        cfg["folds"] = {pointing: {"init":{}, "post":{}}}
        cfg["folds"][pointing]["classifier"] = 0
        cfg["folds"][pointing]["dir"] = join(cfg["files"]["psr_dir"], pointing)
        for _, i in enumerate(init):
            cfg["folds"][pointing]["init"][str(i)] = {}
        for _, i in enumerate(post):
            cfg["folds"][pointing]["post"][str(i)] = {}

    cfg["pol"]["RM"] = None
    cfg["pol"]["RM_e"] = None
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
    cfg["completed"]["debase"] = False
    cfg["completed"]["RM"] = False
    cfg["completed"]["RVM_initial"] = False
    cfg["completed"]["RVM_final"] = False
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
    cfg["completed"]["debase"] = False
    cfg["completed"]["RM"] = False
    cfg["completed"]["RVM_initial"] = False
    cfg["completed"]["RVM_final"] = False


def from_yaml(filepath):
    with open(filepath) as f:
        my_dict = yaml.load(f, Loader=yaml.Loader)
    return my_dict


def dump_to_yaml(cfg):
    with open(cfg["files"]["my_name"], 'w') as f:
        yaml.dump(cfg, f, default_flow_style=False)
    return cfg["files"]["my_name"]


def create_cfgs_main(kwargs, psrs_pointing_dict):
    """
    uses kwargs from observation_processing_pipeline.py
    """
    metadata, full_meta = get_common_obs_metadata(kwargs["obsid"], return_all=True)
    query = psrqpy.QueryATNF(loadfromdb=data_load.ATNF_LOC).pandas
    cfgs = []
    # Make the fold times dictionary (it's done for all pulsars simultaneously for speed)
    psr_list = list(psrs_pointing_dict.keys())
    fold_times_dict = find_fold_times(psr_list, kwargs["obsid"], kwargs["beg"], kwargs["end"], metadata=metadata, full_meta=full_meta, query=query)
    for psr in progress_bar(psr_list, "Initiating pulsar configs: "):
        logger.info(psr)
        pointing_list = psrs_pointing_dict[psr]
        enter = fold_times_dict[psr]["enter"]
        leave = fold_times_dict[psr]["leave"]
        power = fold_times_dict[psr]["power"]
        try:
            cfgs.append(initiate_cfg(kwargs, psr, pointing_list, enter, leave, power, query=query, metadata=metadata))
        except TypeError as e:
            logger.info(e)
            continue
    return cfgs
