#!/usr/bin/env python3

import os
from os.path import join as ospj
import glob
import logging
import argparse
import sys
import psrqpy
import datetime
import numpy as np

import data_processing_pipeline as dpp
import plotting_toolkit
import find_pulsar_in_obs as fpio
import sn_flux_est as snfe
import stokes_fold
import check_known_pulsars
import prepfold_launch


logger = logging.getLogger(__name__)


class NoSuitableProfileError(Exception):
    """Raise when no suitable profiles are found"""
    pass


def bestprof_info(filename):
    """
    Finds various information on a .bestprof file

    Parameters:
    filename: string
        The path of the bestprof file
    Returns:
    info_dict: dictionary
        A dictionary consisting of the following:
        obsid: int
            The ID of the observation
        puslar: string
            The J name of the pulsar
        nbins: int
            The number of bins used to fold this profile
        chi: float
            The reduced Chi squared value of the fold
        sn: float
            The signal to noise ratio of the fold
        dm: float
            The pulsar's dispersion measure
        period: float
            The pulsar's period
        period_error: float
            The error in the pulsar's period measurement
    """
    # open the file and read the info into a dictionary
    info_dict = {}
    f = open(filename, "r")
    lines = f.read()
    f.close()
    lines = lines.split("\n")
    # info:
    info_dict["obsid"] = int(lines[0].split()[4].split("_")[0])
    info_dict["pulsar"] = lines[1].split()[3].split("_")[1]
    info_dict["nbins"] = int(lines[9].split()[4])
    info_dict["chi"] = float(lines[12].split()[4])
    info_dict["sn"] = float(lines[13].split()[4][2:])
    info_dict["dm"] = float(lines[14].split()[4])
    info_dict["period"] = float(lines[15].split()[4])/1e3  # in seconds
    info_dict["period_error"] = float(lines[15].split()[6])/1e3
    f.close()
    return info_dict


def find_best_pointing(pipe):
    """Finds the pointing directory with the highest S/N fold"""
    if "64" in pipe["folds"]["init"].keys():
        eval_bins = 64
    elif "50" in pipe["folds"]["init"].keys():
        eval_bins = 50
    bestprof_info_list = []
    for pointing in pipe["run_ops"]["dirs"]:
        os.chdir(pointing)
        prof_name = glob.glob(
            f"*b{eval_bins}**{pipe['source']['name'][1:]}*.bestprof")[0]
        bestprof_info_list.append(bestprof_info(prof_name))
    best_sn = 0.0
    best_i = -1
    for i, info_dict in enumerate(bestprof_info_list):
        if info_dict["chi"] >= pipe["run_ops"]["thresh_chi"] and info_dict["sn"] > best_sn and info_dict["sn"] >= pipe["run_ops"]["thresh_sn"]:
            best_sn = info_dict["sn"]
            best_i = i
    if best_i < 0:
        raise NoSuitableProfileError(
            "No profile above the threshold values were found")
    logger.info(f"Pulsar found in pointings {pipe['run_ops']['dirs'][i]}")

    return pipe['run_ops']['dirs'][i]


def bf_init(pipe):
    """Fold on the initial bin counts"""
    fold_ids = []
    for _dir in pipe["run_ops"]["dirs"]:
        for bin_count in pipe["folds"]["init"]:
            kwargs = prepfold_launch.common_kwargs(pipe, int(bin_count))
            _id = prepfold_launch.submit_prepfold(pipe, _dir, kwargs)
            fold_ids.append(_id)
    pipe["completed"]["init_folds"] = True
    return fold_ids


def bf_post(pipe):
    """Fold on the follow-up bin counts should there be a pulsar in the initial folds"""
    pipe["run_ops"]["my_dir"] = find_best_pointing(pipe)
    os.chdir(pipe["run_ops"]["my_dir"])
    if "100" in pipe["folds"]["init"].keys():
        b = 100
    else:
        b = 50
    prof = glob.glob(f"{pipe['run_ops']['my_dir']}/*b{b}*{pipe['source']['name']}*pfd.bestprof")[0]
    info = bestprof_info(prof)
    pipe["source"]["my_DM"] = float(info["dm"])
    pipe["source"]["my_P"] = float(info["period"])
    fold_ids = []
    for bin_count in pipe["folds"]["post"]:
        kwargs = prepfold_launch.common_kwargs(pipe, int(bin_count))
        _id = prepfold_launch.submit_prepfold(
            pipe, pipe["run_ops"]["my_dir"], kwargs)
        fold_ids.append(_id)
    pipe["completed"]["post_folds"] = True
    return fold_ids


def bf_main(pipe):
    """A logic structure that decides what to do next"""
    if not pipe["completed"]["bf"]:
        if not pipe["completed"]["init_folds"]:
            dep_ids = bf_init(pipe)
            dpp.resubmit_self(pipe, dep_ids=dep_ids)
        elif not pipe["completed"]["post_folds"]:
            dep_ids = bf_post(pipe)
            dpp.resubmit_self(pipe, dep_ids=dep_ids)
    else:
        logger.info("The binfinder pipeline has already been completed")