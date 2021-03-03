import logging
import sys

from dpp.helper_config import dump_to_yaml
from dpp.helper_status import status_from_error, message_from_status

logger = logging.getLogger(__name__)


def finish_unsuccessful(cfg, e):
    status = status_from_error(e)
    logger.info("\n")
    logger.info("-------------------------------------------------------------------")
    logger.info(f"Pipeline has completed its run on pulsar {cfg['source']['name']}")
    logger.info(f"Pipeline was terminated early with code {status}:")
    logger.info(f"Message: {e}")
    cfg["run_ops"]["exit_status"] = status
    dump_to_yaml(cfg)
    sys.exit(0)


def finish_successful(cfg):
    logger.info("\n")
    logger.info("-------------------------------------------------------------------")
    logger.info(f"Pipeline has completed its run on pulsar {cfg['source']['name']}")
    logger.info(f"Pipeline successfully completed a full run")
    my_bins = cfg['source']['my_bins']
    my_pointing = cfg['source']['my_pointing']
    logger.info("Run Results:")
    logger.info(f"Pointing used:            {my_pointing}")
    logger.info(f"Bin count used:           {my_bins}")
    logger.info(f"Presto SN for bin count:  {cfg['folds'][my_pointing]['post'][str(my_bins)]['sn']}")
    logger.info(f"Fold DM:                  {cfg['source']['my_DM']}")
    logger.info(f"Fold Period:              {cfg['source']['my_P']}")
    logger.info(f"Fold Period Derivative:   {cfg['source']['my_Pdot']}")
    logger.info(f"Rotation Measure:         {cfg['pol']['RM']} +/- {cfg['pol']['RM_e']}")
    logger.info(f"l0:                       {cfg['pol']['l0']}")
    logger.info(f"pa0:                      {cfg['pol']['pa0']}")
    logger.info(f"beta:                     {cfg['pol']['beta']}")
    logger.info(f"alpha:                    {cfg['pol']['alpha']}")
    logger.info(f"chi:                      {cfg['pol']['chi']}")
    cfg["run_ops"]["exit_status"] = "100" # Successful finish
    dump_to_yaml(cfg)
    sys.exit(0)