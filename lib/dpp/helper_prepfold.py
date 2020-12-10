#!/usr/bin/env python3
import logging

from vcstools.config import load_config_file

comp_config = load_config_file()
logger = logging.getLogger(__name__)


def common_kwargs(pipe, bin_count):
    """Creates a prepfold-friendly dictionary of common arguments to pass to prepfold"""
    name = f"{pipe['obs']['id']}_{pipe['run_ops']['pointing']}_{pipe['source']['name']}_b{bin_count}"
    prep_kwargs = {}
    if pipe["run_ops"]["mask"]:
        prep_kwargs["-mask"] = pipe["run_ops"]["mask"]
    prep_kwargs["-o"] = name
    prep_kwargs["-n"] = bin_count
    prep_kwargs["-start"] = pipe["source"]["enter_frac"]
    prep_kwargs["-end"] = pipe["source"]["exit_frac"]
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
    if pipe["source"]["ATNF_P"] < 0.005:  # period less than 50ms
        prep_kwargs["-npfact"] = 4
        prep_kwargs["-ndmfact"] = 3
        prep_kwargs["-dmstep"] = 3
        prep_kwargs["-npart"] = 40
    prep_kwargs["-dm"] = pipe["source"]["ATNF_DM"]
    prep_kwargs["-p"] = pipe["source"]["ATNF_P"]
    if pipe["source"]["my_DM"]:
        prep_kwargs["-dm"] = pipe["source"]["my_DM"]
    if pipe["source"]["my_P"]:
        prep_kwargs["-p"] = pipe["source"]["my_P"]
    return prep_kwargs


def add_prepfold_to_commands(prep_kwargs, eph=None, eph_name=None):
    """Adds prepfold commands to a list. If eph is not None, will use -par"""
    commands = []
    options = ""
    for key, val in prep_kwargs.items():
            options += f" {key} {val}"
    if eph and eph_name:
        with open(eph_name, "w") as f:
            f.write(eph)
        options += f" -par {eph_name}"
    options += " *fits"
    commands.append(f"prepfold {options}")
    return commands


def write_cmd_to_file(pipe, commands, name):
    """Writes the prepfold command to a text file"""
    with open(name, 'w') as f:
        for cmd in commands:
            f.write(cmd)


def prepfold_cmd_make_main(kwargs):
    """takes kwargs from prepfold_cmd_make"""
    from dpp.helper_yaml import from_yaml, dump_to_yaml
    pipe = from_yaml(kwargs["yaml"])
    folds = []
    if not pipe["completed"]["init_folds"]:
        logger.info("Creating command for inital folds")
        for i in pipe["folds"]["init"].keys():
            folds.append(int(i))
        pipe["completed"]["init_folds"] = True
    elif not pipe["completed"]["post_folds"]:
        logger.info("Creating command for post folds")
        for i in pipe["folds"]["post"].keys():
            folds.append(int(i))
        pipe["completed"]["post_folds"] = True
    for bin_count in folds:
        prep_kwargs = common_kwargs(pipe, bin_count)
        cmd = add_prepfold_to_commands(prep_kwargs, eph=pipe["source"]["edited_eph"], eph_name=pipe["source"]["edited_eph_name"])
        name = f"prepfold_cmd_{pipe['run_ops']['file_precursor']}_{bin_count}b.sh"
        write_cmd_to_file(pipe, cmd, name)
    #update yaml file
    dump_to_yaml(pipe, label=kwargs["label"])