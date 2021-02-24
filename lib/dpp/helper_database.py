from os.path import join
import logging
from glob import glob

from vcstools.config import load_config_file
from vcstools.job_submit import submit_slurm
from dpp.helper_files import glob_pfds

comp_config = load_config_file()
logger = logging.getLogger(__name__)


def submit_prepfold_products_db(cfg, dep_id=None, dep_type="afterany"):
    """Submits the best fold profile to the pulsar database. Will also submit .ppps"""
    my_pointing = cfg["source"]["my_pointing"]
    # We will upload the init fold and the best post fold
    bin_list = list(cfg["folds"][my_pointing]["init"].keys())
    bin_list.append(cfg["source"]["my_bins"])
    jids = []
    for bin_count in bin_list:
        commands = []
        commands.append(f"cd {cfg['files']['psr_dir']}")
        # Get the files to upload
        try:
            ppps = glob_pfds(cfg, my_pointing, bin_count, pfd_type=".ps")[0]
        except IndexError as e:
            raise IndexError(f"No ppps files found in dir: {cfg['files']['psr_dir']} for pointing {my_pointing} and bin count {bin_count}")
        try:
            bestprof = glob_pfds(cfg, my_pointing, bin_count, pfd_type=".bestprof")[0]
        except IndexError as e:
            raise IndexError(f"No bestprof files found in dir: {cfg['files']['psr_dir']} for pointing {my_pointing} and bin count {bin_count}")
        commands.append(f"echo 'Submitting profile to database with {bin_count} bins'")
        commands.append(f"submit_to_database.py -o {cfg['obs']['id']} --cal_id {cfg['obs']['cal']} -p {cfg['source']['name']} --bestprof {bestprof} --ppps {ppps}")

        # Submit this job
        name = f"Submit_db_{cfg['files']['file_precursor']}_{bin_count}"
        batch_dir = join(comp_config['base_data_dir'], cfg['obs']['id'], "batch")
        this_id = submit_slurm(name, commands,
                        batch_dir=batch_dir, slurm_kwargs={"time": "00:30:00"}, depend=dep_id,
                        module_list=[f"mwa_search/{cfg['run_ops']['mwa_search']}"],
                        vcstools_version=cfg["run_ops"]["vcstools"], submit=True, depend_type=dep_type)

        jids.append(this_id)
        logger.info(f"Submission script on queue for profile: {bestprof}")
        logger.info(f"Job Name: {name}")
        logger.info(f"Job ID: {this_id}")
        cfg["completed"]["upload"] = True
    return jids