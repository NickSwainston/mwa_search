#!/usr/bin/env python3
import logging

from job_submit import submit_slurm

logger = logging.getLogger(__name__)

def launch_label(cfg):
    """Returns a label based on how far the piepline has progressed"""
    order = ["initial", "classify", "post", "upload", "polarimetry"]
    counter = sum(cfg["completed"].values()) # Number of True statements
    return order[counter]


def relaunch_ppp(cfg, depends_on=None, depend_type="afterany", fresh_run=False, reset_logs=False):
    """Relaunches the pulsar processing pipeline using the supplied cfg file"""
    label = launch_label(cfg)
    name = f"ppp_{cfg['run_ops']['file_precursor']}_{label}"
    batch_dir = cfg['run_ops']['batch_dir']
    slurm_kwargs = {"time": "00:30:00"}
    mem=8192
    ppp_launch = "pulsar_processing_pipeline.py"
    ppp_launch += f" --cfg {cfg['run_ops']['myname']}"
    if fresh_run:
        ppp_launch += " --force_rerun"
    if reset_logs:
        ppp_launch += " --reset_logs"
    cmds = [f"cd {cfg['run_ops']['psr_dir']}"]
    cmds.append(ppp_launch)
    modules = [f"mwa_search/{cfg['run_ops']['mwa_search']}"]
    jid = submit_slurm(name, cmds,
            slurm_kwargs=slurm_kwargs, module_list=modules, mem=mem, batch_dir=cfg["run_ops"]["batch_dir"], depend=depends_on,
            depend_type=depend_type, vcstools_version=cfg["run_ops"]["vcstools"], submit=True)
    logger.info(f"Submitted relaunch of ppp: {name}")
    logger.info(f"job ID: {jid}")
    if depends_on:
        logger.info(f"Job depends on job id(s): {depends_on}")