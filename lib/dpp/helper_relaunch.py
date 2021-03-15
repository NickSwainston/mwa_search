import logging

from vcstools.job_submit import submit_slurm
from dpp.helper_config import dump_to_yaml

logger = logging.getLogger(__name__)


def launch_label(cfg):
    """Returns a label based on how far the piepline has progressed"""
    order = ["Initial", "Classify", "Post", "Upload", "Debase", "RM", "RVM_initial", "RVM_final", "Finish"]
    counter = sum(cfg["completed"].values()) # Number of True statements
    return order[counter]


def relaunch_ppp(cfg, depends_on=None, depend_type="afterany", fresh_run=False, reset_logs=False, time="00:30:00"):
    """Relaunches the pulsar processing pipeline using the supplied cfg file"""
    # Dump the new cfg
    dump_to_yaml(cfg)
    label = launch_label(cfg)
    name = f"ppp_{label}_{cfg['files']['file_precursor']}"
    slurm_kwargs = {"time": time}
    mem=8192
    ppp_launch = "pulsar_processing_pipeline.py"
    ppp_launch += f" --cfg {cfg['files']['my_name']}"
    if fresh_run:
        ppp_launch += " --fresh_run"
    if reset_logs:
        ppp_launch += " --reset_logs"
    cmds = [f"cd {cfg['files']['psr_dir']}"]
    cmds.append(ppp_launch)
    modules = [f"mwa_search/{cfg['run_ops']['mwa_search']}", "singularity"]
    jid = submit_slurm(name, cmds,
            slurm_kwargs=slurm_kwargs, module_list=modules, mem=mem, batch_dir=cfg["files"]["batch_dir"], depend=depends_on,
            depend_type=depend_type, vcstools_version=cfg["run_ops"]["vcstools"], submit=True)
    logger.info(f"Submitted relaunch of ppp: {name}")
    logger.info(f"job ID: {jid}")
    if depends_on:
        logger.info(f"Job depends on job id(s): {depends_on}")