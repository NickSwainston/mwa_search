#!/usr/bin/env python3
import logging
import argparse

from dpp.helper_polarimetry import pulsar_polarimetry_main

logger = logging.getLogger(__name__)


def main(kwargs):
    pulsar_polarimetry_main(kwargs)


if __name__ == '__main__':
    loglevels = dict(DEBUG=logging.DEBUG,
                     INFO=logging.INFO,
                     WARNING=logging.WARNING,
                     ERROR=logging.ERROR)
    parser = argparse.ArgumentParser(description="""Creates a prepfold command for each pulsar and writes to file""",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--yaml", type=str, required=True, help="The pathname of the yaml file")
    parser.add_argument("--fits", type=str, nargs="*", help="The name(s) of the fits file(s)")
    parser.add_argument("--loglvl", type=str, default="INFO", help="Logger verbosity level", choices=loglevels.keys())
    parser.add_argument("--label", type=str, default="pulsar_polarimetry", help="A label to apply to the .yaml file")
    args = parser.parse_args()
    logger = logging.getLogger()
    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    formatter = logging.Formatter(
        '%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    kwargs = vars(args)
    main(kwargs)