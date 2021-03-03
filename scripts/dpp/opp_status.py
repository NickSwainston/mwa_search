#!/usr/bin/env python3
import logging
import argparse

from dpp.helper_status import opp_status

logger = logging.getLogger(__name__)


def main(kwargs):
    logger.info(f"Reading config files from obsid {kwargs['obsid']}")
    status = opp_status(kwargs["obsid"])
    total_psrs = sum([len(status[i]) for i in status.keys()])
    status_break = "----------------------------------------------"
    if status["100"]:
        logger.info(status_break)
        logger.info("The following pulsars completed their run successfully:")
        logger.info("\n".join(status["100"]))
    if status["101"]:
        logger.info(status_break)
        logger.info("The following pulsars were detected but could not proceed with an RVM fit:")
        logger.info("\n".join(status["101"]))
    if status["200"]:
        logger.info(status_break)
        logger.info("The following pulsars were not successfully detected:")
        logger.info("\n".join(status["200"]))
    if status["201"]:
        logger.info(status_break)
        logger.info("The following pulsars are missing a .cfg file:")
        logger.info("\n".join(status["201"]))
    if status["300"]:
        logger.info(status_break)
        logger.info("Something unkonwn went wrong while trying to process the following pulsars:")
        logger.info("\n".join(status["300"]))
    if status["400"]:
        logger.info(status_break)
        logger.info("I don't know what happened with these pulsars:")
        logger.info("\n".join(status["400"]))
    logger.info(status_break)
    logger.info(f"Total pulsars run through opp:    {total_psrs}")
    logger.info(f"Total successful detections:      {len(status['100']) + len(status['101'])}")
    logger.info(status_break)


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