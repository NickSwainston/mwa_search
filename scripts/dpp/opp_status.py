#!/usr/bin/env python3
import logging
import argparse

from dpp.helper_status import opp_status

logger = logging.getLogger(__name__)


def main(kwargs):
    logger.info(f"Reading config files from obsid {kwargs['obsid']}")
    status = opp_status(kwargs["obsid"])
    total_psrs = sum([len(status[i]) for i in status.keys()])
    if status["2"]:
        logger.info("The following pulsars completed their run successfully:")
        logger.info(status["2"])
    if status["1"]:
        logger.info("The following pulsars did were not successfully detected:")
        logger.info(status["1"])
    if status["0"]:
        logger.info("The following pulsars are missing a .cfg file:")
        logger.info(status["1"])
    if status["3"]:
        logger.info("Something went wrong while trying to process the following pulsars:")
        logger.info(status["3"])
    logger.info(f"Total pulsars run through opp:    {total_psrs}")
    logger.info(f"Total successful detections:      {len(status['2'])}")


if __name__ == '__main__':

    loglevels = dict(DEBUG=logging.DEBUG,
                     INFO=logging.INFO,
                     WARNING=logging.WARNING,
                     ERROR=logging.ERROR)

    parser = argparse.ArgumentParser(description="""Initialises and launches a pipeline to process beamformed VCS pulsar data for an entire observation""",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group("Required Options")
    required.add_argument("-o", "--obsid", type=str,help="The obs ID of the opp run")

    otherop = parser.add_argument_group("Other Options")
    otherop.add_argument("-L", "--loglvl", type=str, default="INFO", help="Logger verbosity level", choices=loglevels.keys())

    args = parser.parse_args()
    logger = logging.getLogger()
    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    formatter = logging.Formatter(
        '%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    kwargs = vars(args)
    main(kwargs)