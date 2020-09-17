#!/usr/bin/env python

import logging
import argparse
from dpp.helper_yaml import create_yaml_main

logger = logging.getLogger(__name__)


def main(kwargs):
    create_yaml_main(kwargs)


if __name__ == '__main__':
    loglevels = dict(DEBUG=logging.DEBUG,
                     INFO=logging.INFO,
                     WARNING=logging.WARNING,
                     ERROR=logging.ERROR)
    parser = argparse.ArgumentParser(description="""Initialises a .yaml file with all pulsar info required for a DPP run""",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    obsop = parser.add_argument_group("Observation Options")
    obsop.add_argument("-o", "--obsid", type=int, required=True,
                       help="The obs ID of the data")
    obsop.add_argument("-O", "--cal_id", type=int, required=True,
                       help="The ID of the calibrator used to calibrate the data")
    obsop.add_argument("-p", "--psr", type=str, nargs='*',
                       help="The J name of the pulsar(s). e.g. J2241-5236")
    obsop.add_argument("--obs_beg", type=int, required=True,
                       help="The beginning of the observation")
    obsop.add_argument("--obs_end", type=int, required=True, help="The end of the observation")
    obsop.add_argument("--pointing", type=str, required=True, nargs='*',
                       help="The pointing location of the source in the format HH:MM:SS_+DD:MM:SS. e.g. '19:23:48.53_-20:31:52.95'")

    foldop = parser.add_argument_group("Folding/processing Options")
    foldop.add_argument("--sn_min_thresh", type=float, default=8.0, help="The presto sigma value\
                             above which is deemed a detection.")
    foldop.add_argument("--sn_max_thresh", type=float, default=20.0, help="The presto sigma value\
                             above which is deemed a GOOD detection.")
    foldop.add_argument("--chi_thresh", type=float, default=4.0, help="The presto 'chi' value\
                             above which is deemed a detection.")
    foldop.add_argument("--rvmres", type=int, default=90,
                        help="The number of degree samples to try for alpha and beta.")
    foldop.add_argument("--mask", type=str,
                        help="The pathname of the mask to use for folding")

    otherop = parser.add_argument_group("Other Options")
    otherop.add_argument("--no_ephem", action="store_true",
                         help="Use this to override the use of an ephemeris for foldign the pulsar")
    otherop.add_argument("-L", "--loglvl", type=str, default="INFO",
                         help="Logger verbosity level", choices=loglevels.keys())
    otherop.add_argument("--mwa_search", type=str, default="master",
                         help="The version of mwa_search to use")
    otherop.add_argument("--vcstools", type=str, default="master",
                         help="The version of vcs_tools to use")
    otherop.add_argument("--cand", action="store_true",
                         help="use this tag if this is not a kown pulsar")
    otherop.add_argument("--label", type=str, default="make_pulsar_yaml", help="A label to apply to the .yaml file")
    args = parser.parse_args()
    logger = logging.getLogger()
    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    formatter = logging.Formatter(
        '%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    kwargs = vars(args)
    main(kwargs)