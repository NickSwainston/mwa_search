from vcstools import data_load
import psrqpy

import logging
import math

logger = logging.getLogger(__name__)


def bin_sampling_limit(pulsar, sampling_rate=1e-4, query=None):
    """Finds the sampling limit of the input pulsar in units of number of bins"""
    if query is None:
        query = psrqpy.QueryATNF(params=["P0"], psrs=[
                                 pulsar], loadfromdb=data_load.ATNF_LOC).pandas
    query_index = list(query["JNAME"]).index(pulsar)
    period = query["P0"][query_index]
    bin_lim = math.ceil(period/sampling_rate) #round up the limit
    return bin_lim


def is_binary(pulsar, query=None):
    """Checks the ATNF database to see if a pulsar is part of a binary system"""
    if query is None:
        query = psrqpy.QueryATNF(params=["BINARY"], psrs=[
                                 pulsar], loadfromdb=data_load.ATNF_LOC).pandas
    query_index = list(query["JNAME"]).index(pulsar)
    if isinstance(query["BINARY"][query_index], str):
        return True
    else:
        return False


def required_bin_folds(pulsar, sampling_rate=1e-4, query=None):
    """Generates a list of integers that are the folding bins required to complete for dpp"""
    sam_lim = bin_sampling_limit(pulsar, sampling_rate=sampling_rate, query=query)
    if sam_lim >= 1024:  # regular period pulsar
        init_folds = [100]
        post_folds = [1024, 512, 256, 128]
    elif sam_lim > 100 and sam_lim < 1024:  # moderate period pulsar
        init_folds = [100]
        f = sam_lim
        post_folds = []
        while f > 100:
            post_folds.append(f)
            f = int(f/2)
    elif sam_lim <= 100:  # if msp
        init_folds = [50]
        post_folds = [sam_lim]
    return init_folds, post_folds