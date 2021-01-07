#!/usr/bin/env python3
from os.path import join
import logging
from glob import glob

from vcstools.config import load_config_file


comp_config = load_config_file()
logger = logging.getLogger(__name__)


def submit_prepfold_products_db(cfg, dep_id=None, dep_type="afterany"):
    """Submits the best fold profile to the pulsar database. Will also submit .ppps"""
    my_pointing = cfg["source"]["my_pointing"]
    # We will upload the init fold and the best post fold
    bin_list = list(cfg["folds"][my_pointing]["init"].keys()).append(cfg["source"]["my_bins"])
    pointing_dir = join(cfg["run_ops"]["psr_dir"], cfg["run_ops"]["my_pointing"])
    jids = []
    for bin_count in bin_list:
        commands = []
        commands.append(f"cd {pointing_dir}")
        # Get the files to upload
        ppps = glob(f"*{cfg['run_ops']['file_precursor']}*_b{bin_count}*.pfd.ps")[0]
        bestprof = glob(f"*{cfg['run_ops']['file_precursor']}*_b{bin_count}*.pfd.bestprof")[0]
        commands.append(f"echo 'Submitting profile to database with {bin_count} bins'")
        commands.append(f"submit_to_database.py -o {cfg['obs']['id']} --cal_id {cfg['obs']['cal']} -p {cfg['source']['name']} --bestprof {bestprof} --ppps {ppps}")

        # Submit this job
        name = f"Submit_db_{cfg['run_ops']['file_precursor']}"
        batch_dir = join(comp_config['base_data_dir'], cfg['obs']['id'], "batch")
        jids.append(submit_slurm(name, commands,
                        batch_dir=batch_dir, slurm_kwargs={"time": "00:30:00"}, depend=dep_id,
                        module_list=[f"mwa_search/{cfg['run_ops']['mwa_search']}"],
                        vcstools_version=cfg["run_ops"]["vcstools"], submit=True, depend_type=dep_type))

        logger.info(f"Submission script on queue for profile: {bestprof}")
        logger.info(f"Job Name: {name}")
        logger.info(f"Job ID: {this_id}")
    return jids