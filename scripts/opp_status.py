#!/usr/bin/env python3
import logging
import argparse

from dpp.helper_status import opp_status

logger = logging.getLogger(__name__)


def main(kwargs):
    logger.info(f"Reading config files from obsid {kwargs['obsid']}")
    status, detections = opp_status(kwargs["obsid"])
    total_psrs = sum([len(status[i]) for i in status.keys()])
    status_break = "\n########################################################################\n"
    if detections:
        p_string = '\n'.join(detections)
        print("status_break")
        logger.info(f"The following pulsars were detected: \n{p_string}")
    if status["100"]:
        p_string = '\n'.join(status['100'])
        print(status_break)
        logger.info(f"100 || The following pulsars completed their run successfully: \n{p_string}")
    if status["101"]:
        p_string = '\n'.join(status['101'])
        print(status_break)
        logger.info(f"101 || The following pulsars could not proceed with an RVM fit: \n{p_string}")
    if status["102"]:
        p_string = '\n'.join(status['102'])
        print(status_break)
        logger.info(f"102 || The following pulsars did not fold with enough bins for a Gaussian fit: \n{p_string}")
    if status["103"]:
        p_string = '\n'.join(status['103'])
        print(status_break)
        logger.info(f"103 || The following pulsars were too noisy to proceed with a Gaussian fit: \n{p_string}")
    if status["200"]:
        p_string = '\n'.join(status['200'])
        print(status_break)
        logger.info(f"200 || The following pulsars were not classified as a detection: \n{p_string}")
    if status["201"]:
        p_string = '\n'.join(status['201'])
        print(status_break)
        logger.info(f"201 || The following pulsars are missing one or more .pfd files: \n{p_string}")
    if status["202"]:
        p_string = '\n'.join(status['202'])
        print(status_break)
        logger.info(f"202 || The following pulsars were missing classifier files: \n{p_string}")
    if status["203"]:
        p_string = '\n'.join(status['203'])
        print(status_break)
        logger.info(f"203 || The following pulsars attempted to be classified with an incompatible bin count: \n{p_string}")
    if status["204"]:
        p_string = '\n'.join(status['204'])
        print(status_break)
        logger.info(f"204 || The following pulsars were missing one or more .fits files: \n{p_string}")
    if status["205"]:
        p_string = '\n'.join(status['205'])
        print(status_break)
        logger.info(f"205 || The following pulsars were missing one or more pointing directories: \n{p_string}")
    if status["300"]:
        p_string = '\n'.join(status['300'])
        print(status_break)
        logger.info(f"300 || The pipeline was initiated but not started for the following pulsars: \n{p_string}")
    if status["400"]:
        p_string = '\n'.join(status['400'])
        print(status_break)
        logger.info(f"400 || I don't know what happened with these pulsars: \n{p_string}")

    print(status_break)
    logger.info(f"Total pulsars run through opp:    {total_psrs}")
    logger.info(f"Total successful detections:      {len(detections)}")
    print(status_break)


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