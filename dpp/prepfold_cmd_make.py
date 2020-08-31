#!/usr/bin/env python3
import os
import datetime
import logging
import argparse
import subprocess

import yaml_helper
from config_vcs import load_config_file

comp_config = load_config_file()
logger = logging.getLogger(__name__)


def common_kwargs(pipe, bin_count):
    """Creates a prepfold-friendly dictionary of common arguments to pass to prepfold"""
    name = f"{pipe['obs']['id']}_b{bin_count}_{pipe['source']['name']}"
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


def create_edited_eph(pulsar_name, eph_name):
    """Created a string version of 'psrcat -e' and removes the last line"""
    eph = subprocess.check_output(["psrcat", "-e", pulsar_name])
    eph = eph.decode("utf-8")
    eph = "\n".join(tuple(a.split("\n")[:-2]))
    return eph


def add_prepfold_to_commands(prep_kwargs, eph=None, eph_name=None):
    """Adds prepfold commands to a list. If eph is not None, will use -par"""
    commands = []
    options = ""
    for key, val in prep_kwargs.items():
            options += f" {key} {val}"
    if eph and eph_name:
        with open(eph_name, "w") as f:
            f.write(eph)
        options += f"-par {eph_name}"
    options += " *fits"
    commands.append("prepfold {}".format(options))
    return commands


def write_cmd_to_file(pipe, commands):
    """Writes the prepfold command to a text file"""
    with open(f"prepfold_cmd_{pipe['run_ops']['file_precursor']}.sh", 'w') as f:
        for cmd in commands:
            f.write(cmd)


def main(pipe, label=""):
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
        cmd = add_prepfold_to_commands(prep_kwargs, pipe["source"]["edited_eph"], pipe["source"]["edited_eph"]["name"])
        write_cmd_to_file(pipe, cmd)
    #update yaml file
    yaml_helper.dump_to_yaml(pipe, label=label)


if __name__ == '__main__':
    loglevels = dict(DEBUG=logging.DEBUG,
                     INFO=logging.INFO,
                     WARNING=logging.WARNING,
                     ERROR=logging.ERROR)
    parser = argparse.ArgumentParser(description="""Creates a prepfold command for each pulsar and writes to file""",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--yaml", type=str, required=True, help="The pathname of the yaml file")
    parser.add_argument("--loglvl", type=str, default="INFO", help="Logger verbosity level", choices=loglevels.keys())
    parser.add_argument("--label", type=str, default="", help="A label to apply to the .yaml file")
    args = parser.parse_args()
    logger = logging.getLogger()
    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    formatter = logging.Formatter(
        '%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    pipe = yaml_helper.from_yaml(args.yaml)
    main(pipe, label=args.label)