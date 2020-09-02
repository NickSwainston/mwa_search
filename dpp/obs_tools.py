from astropy.coordinates import SkyCoord
import astropy.units as u

# vcstools imports
import checks
from config_vcs import load_config_file

# mwa_search imports
from grid import get_grid

import logging
logger = logging.getLogger(__name__)

def check_data(obsid, beg=None, dur=None, base_dir=None):
    """
    Checks to see if all of the recombined files exist on disk

    Parameters:
    -----------
    obsid: int
        The observation ID to check
    beg: int
        OPTIONAL - The beginning time of the files to check. If none, will use entire obs. Default: None
    dur: int
        OPTIONAL - The duration in seconds to check since the beginning time. If none, will use entire obs. Default: None
    base_dir: string
        OPTIONAL - The base directory to use. If none, will load from config. Default: None

    Returns:
    ---------
    check: boolean
        True - all files are on disk. False - not all files are on disk
    """
    if base_dir is None:
        comp_config = load_config_file()
        base_dir = comp_config['base_data_dir']
    comb_dir = "{0}{1}/combined".format(base_dir, obsid)

    if not isinstance(beg, int):
        beg = int(beg)
    if not isinstance(dur, int):
        dur = int(dur)

    #Check to see if the files are combined properly
    if beg is not None and dur is not None:
        logger.info("Checking recombined files beginning at {0} and ending at {1}. Duration: {2} seconds"\
                    .format(beg, (beg+dur), dur))
        error = checks.check_recombine(obsid, startsec=beg, n_secs=dur, directory=comb_dir)
    else:
        logger.warn("No start time information supplied. Comparing files with full obs")
        error = checks.check_recombine(obsid, directory=comb_dir)

    if error == True:
        check = False
    else:
        check = True
    return check


def calc_ta_fwhm(freq, array_phase='P2C'):
    """
    Calculates the approximate FWHM of the tied array beam in degrees.

    Parameters:
    -----------
    freq: float
        Frequency in MHz
    array_phase: string
        OPTIONAL - The different array phase (from P1, P2C, P2E) to work out the maximum baseline length. Default = 'P2C'

    Returns:
    --------
    fwhm: float
        FWHM in degrees
    """
    from scipy.constants import c
    from math import degrees

    # Work out baseline in meters
    if array_phase == 'P1':
        # True max_baseline is 2800 but due to the minimal amount of long baselines
        # the following is more realisitic
        max_baseline = 2200.
    elif array_phase == 'P2C':
        # True max_baseline is 700.
        max_baseline = 360.
    elif array_phase == 'P2E':
        max_baseline = 5300.
    else:
        logger.warn("Array phase is OTHER so assumeing the array is in phase 2 extended mode.")
        max_baseline = 5300.

    wavelength = c / (freq * 1e6)
    fwhm = degrees(wavelength / max_baseline)

    return fwhm