#!/usr/bin/env python3
import logging
import argparse
import sys
import os

from dpp.helper_config import from_yaml, dump_to_yaml, reset_cfg
from dpp.helper_prepfold import ppp_prepfold
from dpp.helper_classify import classify_main, read_classifications
from dpp.helper_bestprof import find_best_pointing, populate_post_folds, best_post_fold
from dpp.helper_logging import initiate_logs
from dpp.helper_terminate import finish_unsuccessful, finish_successful
from dpp.helper_files import remove_old_results
from dpp.helper_database import submit_prepfold_products_db
from dpp.helper_archive import ppp_archive_creation, ppp_baseline_removal
from dpp.helper_relaunch import relaunch_ppp
from dpp.helper_RM import RM_synth, RM_cor
from dpp.helper_RVMfit import RVM_fit, RVM_file_to_cfg
from dpp.helper_checks import check_pipe_integrity
# Custom Errors
from vcstools.prof_utils import ProfileLengthError, NoFitError
from dpp.helper_checks import InvalidPAFileError, FitsNotFoundError, PFDNotFoundError
from dpp.helper_bestprof import NoUsableFoldsError

logger = logging.getLogger(__name__)


def main(kwargs):
    """Initiates the pipeline run for a single pulsar"""
    cfg = from_yaml(kwargs["cfg"])

    cfg["run_ops"]["exit_status"] = "400" # If it breaks unknowingly, this is the correct errorcode
    dump_to_yaml(cfg) # Save the errorcode

    # Initiate logging
    writemode = "a"
    if kwargs["reset_logs"]:
        writemode = "w+"
    initiate_logs(cfg["run_ops"]["loglvl"], outfile=cfg["files"]["logfile"], writemode=writemode, stderr=True)

    # Check if a fresh run is forced
    if kwargs["fresh_run"]:
        logger.info("Forcing a fresh run")
        reset_cfg(cfg)
        remove_old_results(cfg) # Remove old files so that they don't interfere with this run

    # Run cfg through the checks pipeline
    try:
        check_pipe_integrity(cfg)
    except (InvalidPAFileError, FitsNotFoundError, PFDNotFoundError) as e:
        finish_unsuccessful(cfg, e)

    # Do the next step in the pipeline
    if cfg["completed"]["init_folds"] == False:
        # Do the initial folds
        dep_jids = ppp_prepfold(cfg)
        relaunch_ppp(cfg, depends_on=dep_jids)
    elif cfg["completed"]["classify"] == False:
        # Classify the intial folds
        try:
            dep_jid = classify_main(cfg)
        except InvalidClassifyBinsError as e:
            finish_unsuccessful(cfg, e)
        relaunch_ppp(cfg, depends_on=dep_jid)
    elif cfg["completed"]["post_folds"] == False:
        # Read the output of the classifier
        try:
            read_classifications(cfg)
        except ClassifierFilesNotFoundError as e:
            finish_unsuccessful(cfg, e)
        # Decide on next folds
        try:
            find_best_pointing(cfg)
        except NoUsableFoldsError as e:
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
        dep_jid, _ = ppp_archive_creation(cfg)
        relaunch_ppp(cfg, depends_on=dep_jid)
    elif cfg["completed"]["debase"] == False:
        # Baseline RFI removal
        try:
            dep_jid, _ = ppp_baseline_removal(cfg)
        except (ProfileLengthError, NoFitError) as e:
            finish_unsuccessful(cfg, e)
        relaunch_ppp(cfg, depends_on=dep_jid, time="02:00:00") # RM synth might take a while - give it more time
    elif cfg["completed"]["RM"] == False:
        # Perform RM synthesis
        RM_synth(cfg)
        # Correct for RM
        dep_jid, _ = RM_cor(cfg)
        relaunch_ppp(cfg, depends_on=dep_jid)
    elif cfg["completed"]["RVM_initial"] == False:
        # Initial RVM fit
        dep_jid = RVM_fit(cfg)
        relaunch_ppp(cfg, depends_on=dep_jid)
    elif cfg["completed"]["RVM_final"] == False:
        # Read Initial RVM
        RVM_file_to_cfg(cfg)
        # Final RVM fit
        dep_jid = RVM_fit(cfg)
        relaunch_ppp(cfg, depends_on=dep_jid)
    else:
        # Read Initial RVM
        RVM_file_to_cfg(cfg)
        finish_successful(cfg)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Initialises and launches a pipeline to process beamformed VCS pulsar data for a single pulsar""",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group("Required Options")
    required.add_argument("--cfg", type=str, help="The pathname of the config file for the pulsar")

    runop = parser.add_argument_group("Run Options")
    runop.add_argument("--reset_logs", action="store_true", help="Delete the current log file and make a new one")
    runop.add_argument("--fresh_run", action="store_true", help="Forces a fresh run of the pipeline")
    args = parser.parse_args()
    kwargs = vars(args)
    main(kwargs)