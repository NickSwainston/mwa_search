#!/usr/bin/env python3
import logging
import argparse

from vcstools.progress_bar import progress_bar
from dpp.helper_obs_info import find_pulsars_in_fov, reformat_psrs_pointings
from dpp.helper_files import setup_cfg_dirs, clean_cfg, find_config_files, create_dpp_dir
from dpp.helper_config import create_cfgs_main, dump_to_yaml, from_yaml
from dpp.helper_relaunch import relaunch_ppp
import pulsar_processing_pipeline as ppp

logger = logging.getLogger(__name__)

class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    # This is used for agparse
    pass


def main(kwargs):
    """Initialises the pipeline and begins the run"""
    fresh_run = False
    if kwargs["run_type"] == "rerun_existing" or kwargs["run_type"] == "rerun_broken":
        cfg_names = find_config_files(kwargs["obsid"], kwargs["label"])

        # Separate broken from not broken finishes
        if kwargs["run_type"] == "rerun_broken":
            broken_cfgs = []
            for name in cfg_names:
                cfg = from_yaml(name)
                if cfg["run_ops"]["exit_status"] not in ("100", "101", "200"):
                    broken_cfgs.append(name)

        # Set all 'expected_finishes' to false for the rerun
        elif kwargs["run_type"] == "rerun_existing":
            for name in cfg_names:
                cfg = from_yaml(name)
                cfg["run_ops"]["expected_finish"] = False
                dump_to_yaml(cfg)

    # Create dpp dirs
    elif kwargs["run_type"] == "fresh_run":
        create_dpp_dir(kwargs)

        logger.info("Calculating pulsars in beam...")
        psrs_pointing_list = find_pulsars_in_fov(kwargs["obsid"], kwargs["beg"], kwargs["end"], search_radius=kwargs["search_radius"], offset=kwargs["offset"], angle_offset=kwargs["angle_offset"])
        psrs_pointing_dict = reformat_psrs_pointings(psrs_pointing_list)

        logger.info("Initialising config files...")
        cfgs = create_cfgs_main(kwargs, psrs_pointing_dict)
        cfg_names = []
        for cfg in progress_bar(cfgs, "Setting up config directories: "):
            setup_cfg_dirs(cfg)
            clean_cfg(cfg)
            if cfg: # If there are valid pointing directories
                dump_to_yaml(cfg)
                cfg_names.append(cfg["files"]["my_name"])
        fresh_run = True

    # Launch ppp for each pulsar
    for name in progress_bar(cfg_names, "Launching processing for pulsars: "):
        cfg = from_yaml(name)
        relaunch_ppp(cfg, fresh_run=fresh_run, reset_logs=bool(not kwargs["keep_logs"]))


if __name__ == '__main__':

    loglevels = dict(DEBUG=logging.DEBUG,
                     INFO=logging.INFO,
                     WARNING=logging.WARNING,
                     ERROR=logging.ERROR)

    parser = argparse.ArgumentParser(description="""Initialises and launches a pipeline to process beamformed VCS pulsar data for an entire observation""",
                                     formatter_class=CustomFormatter)
    run_type_help = """
                    fresh_run: Delete all current results in pulsar directories, create new config files and launch ppp for each pulsar.
                    rerun_existing: Load the existing config files for each pulsar and launch the ppp for each of them.
                    rerun_broken: Load the existing config files but only launch the ppp for the files that did not finish properly (error codes other than 100 or 101).
                    """
    required = parser.add_argument_group("Required Options")
    required.add_argument("-o", "--obsid", type=str,help="The obs ID of the data")
    required.add_argument("-O", "--calid", type=str, help="The ID of the calibrator used to calibrate the data")
    required.add_argument("--beg", type=int, help="The beginning of the observation")
    required.add_argument("--end", type=int, help="The end of the observation")
    required.add_argument("--run_type", type=str, choices=["fresh_run", "rerun_existing", "rerun_broken"], help=run_type_help)

    otherop = parser.add_argument_group("Other Options")
    otherop.add_argument("-s", "--search_radius", type=float, default=0.00001,
                         help="The radius to search (create beams within) in degrees to account for ionosphere. Default: 0.00001 degrees (doesn't make a grid)")
    otherop.add_argument("--offset", type=float, default=0, help="The offset to apply to all pointings in arcseconds")
    otherop.add_argument("--angle_offset", type=float, default=0, help="The angle of the offset to apply to all pointings in degrees where zero is north")
    otherop.add_argument("--relaunch", action="store_true", help="Use this tag to relaunch a partially completed opp run.") #TODO: sort out cases for only half present files
    otherop.add_argument("--keep_logs", action="store_true", help="Use this tag to keep the old logs of previous ppp runs")
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