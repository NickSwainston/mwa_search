import logging
import sys

logger = logging.getLogger(__name__)

def finish_unsuccessful(cfg, e):
    logger.info("\n")
    logger.info("-------------------------------------------------------------------")
    logger.info(f"Pipeline has completed its run on pulsar {cfg['source']['name']}")
    logger.info(f"The following tasks were completed")
    for key in cfg["completed"].keys():
        if cfg["completed"][key]:
            logger.info(key)
    logger.info(f"Pipeline was terminated early: {e}")
    sys.exit(0)


def finish_successful(cfg):
    logger.info("\n")
    logger.info("-------------------------------------------------------------------")
    logger.info(f"Pipeline has completed its run on pulsar {cfg['source']['name']}")
    logger.info(f"Pipeline successfully completed a full run")
    logger.info("Run Results:")
    logger.info(f"Pointing used:            {cfg['source']['my_pointing']}")
    logger.info(f"Bin count used:           {cfg['source']['my_bins']}")
    logger.info(f"Fold DM:                  {cfg['source']['my_DM']}")
    logger.info(f"Fold Period:              {cfg['source']['my_P']}")
    logger.info(f"Fold Period Derivative:   {cfg['source']['my_Pdot']}")
    logger.info(f"Rotation Measure:         {cfg['pol']['RM']} +/- {cfg['pol']['RM_e']}")
    logger.info(f"l0:                       {cfg['pol']['l0']}")
    logger.info(f"pa0:                      {cfg['pol']['pa0']}")
    logger.info(f"beta:                     {cfg['pol']['beta']}")
    logger.info(f"alpha:                    {cfg['pol']['alpha']}")
    logger.info(f"chi:                      {cfg['pol']['chi']}")
    sys.exit(0)