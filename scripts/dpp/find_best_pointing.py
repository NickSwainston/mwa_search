#!/usr/bin/env python

import logging
import argparse

from dpp.helper_bestprof import find_best_pointing_main

logger = logging.getLogger(__name__)


def main(kwargs):
    find_best_pointing_main(kwargs)


if __name__ == '__main__':
    loglevels = dict(DEBUG=logging.DEBUG,
                     INFO=logging.INFO,
                     WARNING=logging.WARNING,
                     ERROR=logging.ERROR)
    parser = argparse.ArgumentParser(description="""Finds the best pointing out of a list of .yaml and .pfd files and deletes the rest""",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--pfds", type=str, nargs="+", required=True, help="The .pfd outputs from a PRESTO fold")
    parser.add_argument("--yamls", type=str, nargs="+", required=True, help="The yaml files corresponding to each pointing")
    parser.add_argument("--label", type=str, default="find_best_pointing", help="The label of the new .yaml files")
    parser.add_argument("-L", "--loglvl", type=str, default="INFO", help="Logger verbosity level", choices=loglevels.keys())
    args = parser.parse_args()
    logger = logging.getLogger()
    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    formatter = logging.Formatter(
        '%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    kwargs = vars(args)
    main(kwargs)
