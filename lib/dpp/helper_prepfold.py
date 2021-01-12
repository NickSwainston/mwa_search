#!/usr/bin/env python3
import logging
from os.path import join
from os import chdir, getcwd
import datetime

from vcstools.config import load_config_file
from vcstools.job_submit import submit_slurm
from dpp.helper_config import dump_to_yaml
from dpp.helper_relaunch import relaunch_ppp


comp_config = load_config_file()
logger = logging.getLogger(__name__)


def generate_prep_name(cfg, bins, pointing):
    return f"pf_{cfg['files']['file_precursor']}_{pointing}_b{bins}"


def common_kwargs(cfg, bin_count, pointing):
    """Creates a prepfold-friendly dictionary of common arguments to pass to prepfold"""
    name = generate_prep_name(cfg, bin_count, pointing)
    prep_kwargs = {}
    if cfg["run_ops"]["mask"]:
        prep_kwargs["-mask"] = cfg["run_ops"]["mask"]
    prep_kwargs["-o"] = name
    prep_kwargs["-n"] = bin_count
    prep_kwargs["-start"] = cfg["source"]["enter_frac"]
    prep_kwargs["-end"] = cfg["source"]["exit_frac"]
    prep_kwargs["-runavg"] = ""
    prep_kwargs["-noxwin"] = ""
    prep_kwargs["-noclip"] = ""
    prep_kwargs["-nsub"] = 256
    prep_kwargs["-pstep"] = 1
    prep_kwargs["-pdstep"] = 2
    prep_kwargs["-dmstep"] = 1
    prep_kwargs["-npart"] = 120
    prep_kwargs["-npfact"] = 1
    prep_kwargs["-ndmfact"] = 1
    if bin_count >= 300: #greatly reduces search time
        prep_kwargs["-nopdsearch"] = ""
    if bin_count == 100 or bin_count == 50: #init fold - do large search
        prep_kwargs["-npfact"] = 4
        prep_kwargs["-ndmfact"] = 3
    if cfg["source"]["ATNF_P"] < 0.005:  # period less than 50ms
        prep_kwargs["-npfact"] = 4
        prep_kwargs["-ndmfact"] = 3
        prep_kwargs["-dmstep"] = 3
        prep_kwargs["-npart"] = 40
    prep_kwargs["-dm"] = cfg["source"]["ATNF_DM"]
    prep_kwargs["-p"] = cfg["source"]["ATNF_P"]
    if cfg["source"]["my_DM"]:
        prep_kwargs["-dm"] = cfg["source"]["my_DM"]
    if cfg["source"]["my_P"]:
        prep_kwargs["-p"] = cfg["source"]["my_P"]
    return prep_kwargs


def add_prepfold_to_commands(prep_kwargs, psr_name, pointing, eph=None, eph_name=None, presto_container="", binary=False):
    """Adds prepfold commands to a list. If eph is not None, will use -par"""
    commands = []
    kwarg_cmds = "" # Just the important stuff
    for key, val in prep_kwargs.items():
            kwarg_cmds += f" {key} {val}"
    par_cmds = kwarg_cmds # Everything for -par
    psr_cmds = kwarg_cmds # Everything for -psr
    if eph and eph_name: # For -par (binary)
        par_cmds += f" -par {eph_name}"
    psr_cmds += f" -psr {psr_name}"
    par_cmds += f" {join(pointing, '*fits')}"
    psr_cmds += f" {join(pointing, '*fits')}"
    # Make the commands
    commands.append(f"{presto_container} prepfold {par_cmds} ")
    commands.append("errorcode=$?")
    commands.append('if [ "$errorcode" != "0" ]; then') # Sometimes the -par fails so we use -psr
    commands.append('   echo "Attempting to fold using the -psr option"')
    commands.append(f"   {presto_container} prepfold {psr_cmds}")
    commands.append("   fi")
    return commands


def prepfold_time_alloc(cfg, prepfold_kwargs):
    """Estimates the jobtime for prepfold jobs based on cfg and prepfold args"""
    nopsearch = False
    nopdsearch = False
    nodmsearch = False
    nosearch = False
    if "-nopsearch" in prepfold_kwargs:
        nopsearch = True
    if "-nodmsearch" in prepfold_kwargs:
        nodmsearch = True
    if "-nopdsearch" in prepfold_kwargs:
        nopdsearch = True
    if "-nosearch" in prepfold_kwargs:
        nosearch = True
    npfact = prepfold_kwargs["-npfact"]
    ndmfact = prepfold_kwargs["-ndmfact"]
    bin_count = prepfold_kwargs["-n"]
    duration = (cfg["obs"]["end"] - cfg["obs"]["beg"]) * \
        (cfg["source"]["exit_frac"] - cfg["source"]["enter_frac"])

    time = 600
    time += bin_count
    time += duration

    if not nosearch:
        ptime = 1
        pdtime = 1
        dmtime = 1
        if not nopsearch:
            ptime = npfact*bin_count
        if not nopdsearch:
            pdtime = npfact*bin_count
        if not nodmsearch:
            dmtime = ndmfact*bin_count
        time += ((ptime * pdtime * dmtime)/1e4)
    if time > 86399.:
        logger.warn("Estimation for prepfold time greater than one day")
        time = 86399
    time = str(datetime.timedelta(seconds=int(time)))
    return time


def submit_prepfold(cfg, nbins, pointing, psr_dir, depends_on=None, depend_type="afterany"):
    """Creates the commands for a prepfold job and submits it to the queue"""
    # Make the commands for the job
    prep_kwargs = common_kwargs(cfg, int(nbins), pointing)
    cmds = [f"cd {psr_dir}"]
    cmds += add_prepfold_to_commands(prep_kwargs, cfg["source"]["name"], pointing, eph=cfg["source"]["edited_eph"],
        eph_name=cfg["source"]["edited_eph_name"], presto_container="/pawsey/mwa/singularity/presto/presto.sif", binary=cfg["source"]["binary"])
    # TODO: get rid of the container hard-code ^^

    # Work out some things for job submission
    name = generate_prep_name(cfg, nbins, pointing)
    time = prepfold_time_alloc(cfg, prep_kwargs)
    slurm_kwargs = {"time":time}
    modules = ["singularity"]
    mem=8192

    # Submit Job
    jid = submit_slurm(name, cmds,
        slurm_kwargs=slurm_kwargs, module_list=modules, mem=mem, batch_dir=cfg["files"]["batch_dir"], depend=depends_on,
        depend_type=depend_type, vcstools_version=cfg["run_ops"]["vcstools"], submit=True)
    return jid, name


def initial_folds(cfg):
    """Loops through the required initial folds and submits them to slurm"""
    jids = []
    for pointing in cfg["folds"].keys():
        for nbins in cfg["folds"][pointing]["init"].keys():
            jid, name = submit_prepfold(cfg, nbins, pointing, cfg["files"]["psr_dir"])
            jids.append(jid)
            logger.info(f"Submitted prepfold job: {name}")
            logger.info(f"Job ID: {jid}")
    return jids


def post_folds(cfg):
    """Loops through the required post folds and submits them to slurm"""
    jids = []
    pointing = cfg["source"]["my_pointing"]
    for nbins in cfg["folds"][pointing]["post"].keys():
        jid, name = submit_prepfold(cfg, nbins, pointing, cfg["files"]["psr_dir"])
        jids.append(jid)
        logger.info(f"Submitted prepfold job: {name}")
        logger.info(f"Job ID: {jid}")
    return jids


def ppp_prepfold(cfg):
    from dpp.helper_config import dump_to_yaml
    if not cfg["completed"]["init_folds"]:
        logger.info("Creating commands for inital folds")
        dep_jids = initial_folds(cfg)
        cfg["completed"]["init_folds"] = True

    elif not cfg["completed"]["post_folds"]:
        logger.info("Creating command for post folds")
        dep_jids = post_folds(cfg)
        cfg["completed"]["post_folds"] = True
    return dep_jids