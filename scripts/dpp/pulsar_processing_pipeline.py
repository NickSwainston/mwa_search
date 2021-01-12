#!/usr/bin/env python3
import logging
import argparse
import sys
import os

from dpp.helper_config import from_yaml, reset_cfg
from dpp.helper_prepfold import ppp_prepfold
from dpp.helper_classify import classify_main, read_classifications
from dpp.helper_bestprof import find_best_pointing, NoUsableFolds, populate_post_folds, best_post_fold
from dpp.helper_logging import initiate_logs
from dpp.helper_terminate import finish_unsuccessful
from dpp.helper_files import remove_old_results
from dpp.helper_database import submit_prepfold_products_db
from dpp.helper_archive import ppp_file_creation
from dpp.helper_relaunch import relaunch_ppp
from dpp.helper_RM import RM_synth, RM_cor

logger = logging.getLogger(__name__)


def main(kwargs):
    """Initiates the pipeline run for a single pulsar"""
    cfg = from_yaml(kwargs["cfg"])

    # Initiate logging
    writemode = "a"
    if kwargs["reset_logs"]:
        writemode = "w+"
    initiate_logs(cfg["run_ops"]["loglvl"], outfile=cfg["files"]["logfile"], writemode=writemode, stderr=True)

    # Check if a fresh run is forced
    if kwargs["force_rerun"]:
        logger.info("Forcing a fresh run")
        reset_cfg(cfg)
        remove_old_results(cfg) # Remove old files so that they don't interfere with this run

    # Do the next step in the pipeline
    if cfg["completed"]["init_folds"] == False:
        # Do the initial folds
        dep_jids = ppp_prepfold(cfg)
        relaunch_ppp(cfg, depends_on=dep_jids)
    elif cfg["completed"]["classify"] == False:
        # Classify the intial folds
        dep_jid = classify_main(cfg)
        relaunch_ppp(cfg, depends_on=dep_jid)
    elif cfg["completed"]["post_folds"] == False:
        # Read the output of the classifier
        read_classifications(cfg)
        # Decide on next folds
        try:
            find_best_pointing(cfg)
        except NoUsableFolds as e:
            finish_unsuccessful(cfg, e)
        # Submit post folds
        dep_jids = ppp_prepfold(cfg)
        relaunch_ppp(cfg, depends_on=dep_jids)
    elif cfg["completed"]["upload"] == False:
        # Update cfg with fold info
        populate_post_folds(cfg)
        # Find the best post-fold
        best_post_fold(cfg)
        # Upload stuff to database
        submit_prepfold_products_db(cfg)
        # Launch archive/fits creation job
        dep_jid, _ = ppp_file_creation(cfg)
        relaunch_ppp(cfg, depends_on=dep_jid, time="02:00:00") # RM synth might take a while - give it more time
    elif cfg["completed"]["RM"] == False:
        # Perform RM synthesis
        RM_synth(cfg)
        # Correct for RM
        dep_jid, _ = RM_cor(cfg)
        relaunch_ppp(cfg, depends_on=dep_jid)
    elif cfg["completed"]["RVM_initial"] == False:
        logger.info("TODO: RVM fitting")






if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Initialises and launches a pipeline to process beamformed VCS pulsar data for a single pulsar""",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group("Required Options")
    required.add_argument("--cfg", type=str, help="The pathname of the config file for the pulsar")

    runop = parser.add_argument_group("Run Options")
    runop.add_argument("--reset_logs", action="store_true", help="Delete the current log file and make a new one")
    runop.add_argument("--force_rerun", action="store_true", help="Forces a fresh run of the pipeline")
    args = parser.parse_args()
    kwargs = vars(args)
    main(kwargs)