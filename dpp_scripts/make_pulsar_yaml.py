#!/usr/bin/env python 

import psrqpy
import logging
import yaml
import argparse

from mwa_metadb_utils import get_common_obs_metadata
from vcstools import data_load

from dpp.yaml_helper import initiate_pipe, dump_to_yaml

logger = logging.getLogger(__name__)

def main(kwargs):
    metadata, full_meta = get_common_obs_metadata(kwargs["obsid"], return_all=True)
    query = psrqpy.QueryATNF(loadfromdb=data_load.ATNF_LOC).pandas
    pulsars_pointings_dict = {}
    for psrlist, pointing in zip(kwargs["psrs"], kwargs["pointings"]):
        for psr in psrlist.split(":"):
            if psr not in pulsars_pointings_dict.keys():
                pulsars_pointings_dict[psr] = []
            pulsars_pointings_dict[psr].append(pointing)
    for psr in pulsars_pointings_dict.keys():
        logger.info("Processing yaml for PSR: {}".format(psr))
        pipe = initiate_pipe(kwargs, psr, metadata=metadata, full_meta=full_meta, query=query[query['PSRJ'] == psr].reset_index())
        for pointing in pulsars_pointings_dict[psr]:
            # Update the pipe with the pointing specific parameters
            pipe["run_ops"]["pointing"] = pointing
            if pipe["source"]["cand"] == False:
                pipe["run_ops"]["file_precursor"] = f"{pipe['obs']['id']}_{pipe['run_ops']['pointing']}_{pipe['source']['name']}"
                if pipe["source"]["binary"]:
                    from prepfold_cmd_make import create_edited_eph
                    pipe["source"]["edited_eph_name"] = f"{pipe['run_ops']['file_precursor']}.eph"
                    pipe["source"]["edited_eph"] = create_edited_eph(pipe["source"]["name"], pipe["source"]["edited_eph_name"])
            dump_to_yaml(pipe, label=kwargs["label"])


if __name__ == '__main__':
    loglevels = dict(DEBUG=logging.DEBUG,
                     INFO=logging.INFO,
                     WARNING=logging.WARNING,
                     ERROR=logging.ERROR)
    parser = argparse.ArgumentParser(description="""Initialises a .yaml file with all pulsar info required for a DPP run""",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    obsop = parser.add_argument_group("Observation Options")
    obsop.add_argument("-d", "--run_dir", nargs='+',
                       type=str, help="The location of the pointing directory/s")
    obsop.add_argument("-o", "--obsid", type=str,
                       help="The obs ID of the data")
    obsop.add_argument("-O", "--cal_id", type=str,
                       help="The ID of the calibrator used to calibrate the data")
    obsop.add_argument("-p", "--psrs", nargs="+", type=str,
                       help="The J name of the pulsar(s). e.g. J2241-5236")
    obsop.add_argument("--obs_beg", type=int,
                       help="The beginning of the observation")
    obsop.add_argument("--obs_end", type=int, help="The end of the observation")
    obsop.add_argument("--pointings", type=str, nargs="+", help="The pointing(s) location of the source")

    foldop = parser.add_argument_group("Folding/processing Options")
    foldop.add_argument("--sn_min_thresh", type=float, default=8.0, help="The presto sigma value\
                             above which is deemed a detection.")
    foldop.add_argument("--sn_max_thresh", type=float, default=20.0, help="The presto sigma value\
                             above which is deemed a GOOD detection.")
    foldop.add_argument("--chi_thresh", type=float, default=4.0, help="The presto 'chi' value\
                             above which is deemed a detection.")
    foldop.add_argument("--rvmres", type=int, default=90,
                        help="The number of degree samples to try for alpha and beta.")
    foldop.add_argument("--mask", type=str,
                        help="The pathname of the mask to use for folding")

    otherop = parser.add_argument_group("Other Options")
    otherop.add_argument("--no_ephem", action="store_true",
                         help="Use this to override the use of an ephemeris for foldign the pulsar")
    otherop.add_argument("-L", "--loglvl", type=str, default="INFO",
                         help="Logger verbosity level", choices=loglevels.keys())
    otherop.add_argument("--mwa_search", type=str, default="master",
                         help="The version of mwa_search to use")
    otherop.add_argument("--vcstools", type=str, default="master",
                         help="The version of vcs_tools to use")
    otherop.add_argument("--cand", action="store_true",
                         help="use this tag if this is not a kown pulsar")
    otherop.add_argument("--label", type=str, default="", help="A label to apply to the .yaml file")
    args = parser.parse_args()
    logger = logging.getLogger()
    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    formatter = logging.Formatter(
        '%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    kwargs = vars(args)
    main(kwargs)