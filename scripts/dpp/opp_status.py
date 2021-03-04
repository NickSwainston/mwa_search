#!/usr/bin/env python3
import logging
import argparse

from dpp.helper_status import opp_status

logger = logging.getLogger(__name__)

"""
        "100": "Run completed with detection",
        "101": "Run completed but RVM could not be fit due to unsifficient paswing file",
        "102": "Profile has insufficient bins to continue with Gaussian fit",
        "103": "Profile is too noisy to continue with Gaussian fit",
        "200": "No detections in the run",
        "201": "A PFD file is missing",
        "202": "Classifier files are missing",
        "203": "Bin count for classifer inputs is incorrect",
        "300": "Pipeline has not started",
        "400": "Something unknown went wrong"
"""

def main(kwargs):
    logger.info(f"Reading config files from obsid {kwargs['obsid']}")
    status = opp_status(kwargs["obsid"])
    total_psrs = sum([len(status[i]) for i in status.keys()])
    status_break = "----------------------------------------------"
    if status["100"]:
        p_string = '\n'.join(status['100'])
        logger.info(status_break)
        logger.info(f"100 || The following pulsars completed their run successfully: \n{p_string}")
        logger.info("\n".join(status["100"]))
    if status["101"]:
        p_string = '\n'.join(status['101'])
        logger.info(status_break)
        logger.info(f"101 || The following pulsars were detected but could not proceed with an RVM fit: \n{p_string}")
        logger.info("\n".join(status["101"]))
    if status["102"]:
        p_string = '\n'.join(status['102'])
        logger.info(status_break)
        logger.info(f"102 || The following pulsars were detected but did not fold with enough bins for a Gaussian fit: \n{p_string}")
        logger.info("\n".join(status["102"]))
    if status["103"]:
        p_string = '\n'.join(status['103'])
        logger.info(status_break)
        logger.info(f"103 || The following pulsars were detected but weree too noisy to procees with a Gaussian fit: \n{p_string}")
        logger.info("\n".join(status["103"]))
    if status["200"]:
        p_string = '\n'.join(status['200'])
        logger.info(status_break)
        logger.info(f"200 || The following pulsars were not detected: \n{p_string}")
        logger.info("\n".join(status["200"]))
    if status["201"]:
        p_string = '\n'.join(status['201'])
        logger.info(status_break)
        logger.info(f"201 || The following pulsars are missing a .pfd file: \n{p_string}")
        logger.info("\n".join(status["201"]))
    if status["202"]:
        p_string = '\n'.join(status['202'])
        logger.info(status_break)
        logger.info(f"202 || The following pulsars were missing classifier files: \n{p_string}")
        logger.info("\n".join(status["202"]))
    if status["203"]:
        p_string = '\n'.join(status['203'])
        logger.info(status_break)
        logger.info(f"203 || The following pulsars attempted to be classified with an incompatible bin count: \n{p_string}")
        logger.info("\n".join(status["203"]))
    if status["300"]:
        p_string = '\n'.join(status['300'])
        logger.info(status_break)
        logger.info(f"300 || Something unkonwn went wrong while trying to process the following pulsars: \n{p_string}")
        logger.info("\n".join(status["300"]))
    if status["400"]:
        p_string = '\n'.join(status['400'])
        logger.info(status_break)
        logger.info(f"400 || I don't know what happened with these pulsars: \n{p_string}")

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