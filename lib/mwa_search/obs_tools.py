from mwa_pb.mwa_tile import h2e

from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import astropy.units as u

# vcstools imports
from vcstools.check_files import check_recombine
from vcstools.config import load_config_file

# mwa_search imports
from mwa_search.grid_tools import get_grid

import logging
logger = logging.getLogger(__name__)


def getTargetAZZA(ra,dec,time,lat=-26.7033,lon=116.671,height=377.827):
    """
    Function to get the target position in alt/az at a given EarthLocation and Time.

    Default lat,lon,height is the centre of  MWA.

    Input:
      ra - target right ascension in astropy-readable format
      dec - target declination in astropy-readable format
      time - time of observation in UTC (i.e. a string on form: yyyy-mm-dd hh:mm:ss.ssss)
      lat - observatory latitude in degrees
      lon - observatory longitude in degrees

    Returns:
      a list containing four elements in the following order:
        list[0] = target azimuth in radians
        list[1] = target zenith angle in radians
        list[2] = target azimuth in degrees
        list[3] = target zenith angle in degrees
    """
    #print "Creating EarthLocation for: lat = {0} deg, lon = {1} deg".format(lat,lon)
    location = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=height*u.m)

    #print "Creating SkyCoord for target at (RA,DEC) = ({0},{1})".format(ra,dec)
    coord = SkyCoord(ra,dec,unit=(u.hourangle,u.deg))
    #print "Target at: ({0},{1}) deg".format(coord.ra.deg,coord.dec.deg)

    obstime = Time(time)
    #print "Observation time: {0}".format(obstime.iso)

    #print "Converting to Alt,Az for given time and observatory location..."
    altaz = coord.transform_to(AltAz(obstime=obstime,location=location))
    #print "Target (Alt,Az) = ({0},{1}) deg".format(altaz.alt.deg,altaz.az.deg)

    #print "Converting to (Az,ZA)"
    az = altaz.az.rad
    azdeg = altaz.az.deg

    za = np.pi/2 - altaz.alt.rad
    zadeg = 90 - altaz.alt.deg

    #print "Target (Az,ZA) = ({0},{1}) deg".format(azdeg,zadeg)

    return [az,za,azdeg,zadeg]


def getTargetradec(az,za,time,lst,lat=-26.7033,lon=116.671,height=377.827):
    """
    Function to get the target position in ra dec at a given EarthLocation and Time.

    Default lat,lon,height is the centre of  MWA.

    Input:
      az - target aximuth in radians
      za - target zenith in radians
      time - time of observation in UTC (i.e. a string on form: yyyy-mm-dd hh:mm:ss.ssss)
      lat - observatory latitude in degrees
      lon - observatory longitude in degrees

    Returns:
      a list containing four elements in the following order:
        list[0] = target ra in degrees
        list[1] = target dec in degrees
    """

    ha,dec = h2e(az,za,lat) #hour angle and dec in degrees
    ra = lst-ha


    return [ra,dec]

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
        error = check_recombine(obsid, startsec=beg, n_secs=dur, directory=comb_dir)
    else:
        logger.warn("No start time information supplied. Comparing files with full obs")
        error = check_recombine(obsid, directory=comb_dir)

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