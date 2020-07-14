#!/usr/bin/env python3

import os
import logging
import argparse
from config_vcs import load_config_file
import glob
import sys

from job_submit import submit_slurm
from mwa_metadb_utils import get_common_obs_metadata
import stokes_fold
import binfinder

logger = logging.getLogger(__name__)

class pipe_class:

    def __init__(self, fits_dir=None, obsid=None, cal_id=None, pulsar=None,
                obs_beg=None, obs_end=None, dm=None, period=None,
                obs_ra=None, obs_dec=None, obs_dur=None,
                source_beg=None, source_end=None, central_freq=None,
                loglvl="INFO", mwa_search="master", vcs_tools="master",
                RM=None, RM_err=None, stokes_bins=None, stokes_dep=None,
                no_ephem=False, dspsr_ops="", prep_ops="", rvmres=90, ipfb=False):

        #Obs inormation
        self._pointing_dir   = pointing_dir
        self._cal_id         = cal_id
        self._obsid          = obsid
        self._pulsar         = pulsar
        self._beg            = beg
        self._end            = end
        self._freq           = freq

        #Versions
        self._mwa_search     = mwa_search
        self._vcs_tools      = vcs_tools

        #Run Options
        self._stop           = stop
        self._loglvl         = loglvl
        self._no_ephem       = no_ephem
        self._dm             = dm
        self._period         = period
        self._dspsr_ops      = dspsr_ops
        self._prep_ops       = prep_ops
        self._ipfb           = ipfb

        #Run Utilities
        self._stokes_dep     = stokes_dep

        #Other Parameters
        self._rvmres         = rvmres
        self._threshold      = threshold
        self._subint         = subint
        self._RM             = RM
        self._RM_e              = RM_e
        self._stokes_bins    = stokes_bins
        
    @property
    def obs_beg(self):
        return self._obs_beg
    @obs_beg.setter
    def obs_beg(self, val)
        self._obs_beg = val

    @property
    def obs_end(self):
        return self._obs_end
    @obs_end.setter
    def obs_end(self, val)
        self._obs_end = val

    @property
    def obs_ra(self):
        return self._obs_ra
    @obs_ra.setter
    def obs_ra(self, val)
        self._obs_ra = val

    @property
    def obs_dec(self):
        return self._obs_dec
    @obs_dec.setter
    def obs_dec(self, val)
        self._obs_dec = val

    @property
    def obs_dur(self):
        return self._obs_dur
    @obs_dur.setter
    def obs_dur(self, val)
        self._obs_dur = val

    @property
    def RM(self):
        return self._RM
    @RM.setter
    def RM(self, val):
        self._RM = val

    @property
    def RM_e(self):
        return self._RM_e
    @RM_e.setter
    def RM_e(self, val):
        self._RM_e = val

    @property
    def centre_freq(self):
        return self._freq
    @centre_freq.setter
    def centre_freq(self, val):
        self._centre_freq = val

    def initiate_from_metadata(self):
        md = get_common_obs_metadata(obsid)
        self._obs_ra = md[1]
        self._obs_dec = md[2]
        self._obs_dur = md[3]
        self._obs_freq = md[5]