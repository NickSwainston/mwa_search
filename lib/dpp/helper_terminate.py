import logging
import sys

logger = logging.getLogger(__name__)

def finish_unsuccessful(cfg, e):
    logger.info("\n")
    logger.info("-------------------------------------------------------------------")
    logger.info("\n")
    logger.info(f"Pipeline has completed its run on pulsar {cfg['source']['name']}")
    logger.info(f"The following tasks were completed")
    for key in cfg["completed"].keys():
        if cfg["completed"][key]:
            logger.info(key)
    logger.info(f"Pipeline was terminated early: {e}")
    sys.exit(0)