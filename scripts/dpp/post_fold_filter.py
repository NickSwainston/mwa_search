#!/usr/bin/env python3
import logging
import argparse

from dpp.helper_bestprof import post_fold_filter_main

logger = logging.getLogger(__name__)


def main(kwargs):
    post_fold_filter_main(kwargs)


if __name__ == '__main__':
    loglevels = dict(DEBUG=logging.DEBUG,
                     INFO=logging.INFO,
                     WARNING=logging.WARNING,
                     ERROR=logging.ERROR)
    parser = argparse.ArgumentParser(description="""Creates a prepfold command for each pulsar and writes to file""",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--yamls", type=str, nargs='*', required=True, help="The pathname of the yaml files")
    parser.add_argument("--pfds", type=str, nargs='*', required=True, help="The pathname of the pfd files")
    parser.add_argument("--label", type=str, default="post_fold_filter", help="The label of the new .yaml files")
    parser.add_argument("--loglvl", type=str, default="INFO", help="Logger verbosity level", choices=loglevels.keys())
    args = parser.parse_args()
    logger = logging.getLogger()
    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    formatter = logging.Formatter(
        '%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    kwargs = vars(args)
    main(kwargs)