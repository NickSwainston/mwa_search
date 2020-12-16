#!/usr/bin/env python3
import logging
import argparse

from vcstools.config import load_config_file
from vcstools.progress_bar import progress_bar
from dpp.helper_obs_info import find_pulsars_in_fov, reformat_psrs_pointings
from dpp.helper_files import setup_cfg_dirs, clean_cfg, find_config_files, create_dpp_dir
from dpp.helper_config import create_cfgs_main, dump_to_yaml, from_yaml
from dpp.helper_relaunch import relaunch_ppp
import pulsar_processing_pipeline as ppp

comp_config = load_config_file()
logger = logging.getLogger(__name__)

def main(kwargs):
    """Initialises the pipeline and begins the run"""
    if kwargs["relaunch"]:
        cfg_names = find_config_files(kwargs["obsid"], kwargs["label"])

    else:
        # Create dpp dir
        create_dpp_dir(kwargs)

        logger.info("Calculating pulsars in beam...")
        psrs_pointing_list = find_pulsars_in_fov(kwargs["obsid"], kwargs["beg"], kwargs["end"])
        psrs_pointing_dict = reformat_psrs_pointings(psrs_pointing_list)

        logger.info("Initialising config files...")
        cfgs = create_cfgs_main(kwargs, psrs_pointing_dict)

        cfg_names = []
        for cfg in progress_bar(cfgs, "Setting up config directories: "):
            setup_dpp_dirs(cfgs)
            clean_cfg(cfg)
            if cfg: # If there are valid pointing directories
                dump_to_yaml(cfg, cfg["run_ops"]["my_name"])
                cfg_names.append(cfg["run_ops"]["my_name"])

    # Launch ppp for each pulsar
    for name in progress_bar(cfg_names, "Launching processing for pulsars: "):
        cfg = from_yaml(name)
        relaunch_ppp(cfg, fresh_run=kwargs["force_rerun"], reset_logs=bssssssssssool(not kwargs["keep_logs"]))


if __name__ == '__main__':

    loglevels = dict(DEBUG=logging.DEBUG,
                     INFO=logging.INFO,
                     WARNING=logging.WARNING,
                     ERROR=logging.ERROR)

    parser = argparse.ArgumentParser(description="""Initialises and launches a pipeline to process beamformed VCS pulsar data for an entire observation""",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group("Required Options")
    required.add_argument("-o", "--obsid", type=str,help="The obs ID of the data")
    required.add_argument("-O", "--calid", type=str, help="The ID of the calibrator used to calibrate the data")
    required.add_argument("--beg", type=int, help="The beginning of the observation")
    required.add_argument("--end", type=int, help="The end of the observation")

    otherop = parser.add_argument_group("Other Options")
    otherop.add_argument("--relaunch", action="store_true", help="Use this tag to relaunch a partially completed opp run.")
    otherop.add_argument("--keep_logs", action="store_true", help="Use this tag to keep the old logs of previous ppp runs")
    otherop.add_argument("--force_rerun", action="store_true", help="Forces a fresh run of the pipeline")
    otherop.add_argument("--label", type=str, default="", help="A label to use to identify the results from this run")
    otherop.add_argument("-L", "--loglvl", type=str, default="INFO", help="Logger verbosity level", choices=loglevels.keys())
    otherop.add_argument("--mwa_search", type=str, default="master", help="The version of mwa_search to use")
    otherop.add_argument("--vcstools", type=str, default="master", help="The version of vcs_tools to use")

    args = parser.parse_args()
    logger = logging.getLogger()
    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    formatter = logging.Formatter(
        '%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    kwargs = vars(args)
    main(kwargs)