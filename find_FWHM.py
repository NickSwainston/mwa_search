#!/usr/bin/env python

from mwapy.pb import primary_beam as pb
import argparse
from scipy.interpolate import UnivariateSpline
import numpy as np

import mwa_metadb_utils as meta
#import urllib
#import urllib2
#import json

from astropy.coordinates import EarthLocation, SkyCoord, AltAz
from astropy import units as u
from astropy.time import Time
from decimal import Decimal

def get_obs_metadata(obs):
    from mwapy.pb import mwa_db_query as mwa_dbQ

    beam_meta_data = meta.getmeta(service='obs', params={'obs_id':obs})
    channels = beam_meta_data[u'rfstreams'][u"0"][u'frequencies']
    freqs = [float(c)*1.28 for c in channels]
    xdelays = beam_meta_data[u'rfstreams'][u"0"][u'xdelays']
    ydelays = beam_meta_data[u'rfstreams'][u"0"][u'ydelays']
    pointing_AZ, pointing_EL, pointing_ZA = mwa_dbQ.get_beam_pointing(obs)

    return {"channels":channels,
            "frequencies":freqs,
            "xdelays":xdelays,
            "ydelays":ydelays,
            "az":pointing_AZ,
            "za":pointing_ZA
            }

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

    # Create Earth location for MWA
    location = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=height*u.m)

    # Create sky coordinates for target
    coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))

    # Create a time object for desired observing time
    obstime = Time(time)

    # Convert from RA/Dec to Alt/Az
    altaz = coord.transform_to(AltAz(obstime=obstime, location=location))

    az = altaz.az.rad
    azdeg = altaz.az.deg

    za = np.pi/2 - altaz.alt.rad
    zadeg = 90 - altaz.alt.deg


    return az, za, azdeg, zadeg


parser = argparse.ArgumentParser(description="""
       Finds the FWHM from the output file from pabeam.py
       """)
parser.add_argument('-f', '--file_loc',type=str,help='Input file location')
args=parser.parse_args()

temp = args.file_loc.split("/")
obsid, time, freq, ra, dec = temp[-1].split("_")
freq = float(freq[:-3])
#dec,x = dec.split(".")
dec = dec[:-4]
if "." in time:
   time, x = time.split(".")
time = Time(int(time), format="gps")
az, za, azdeg, zadeg = getTargetAZZA(ra,dec,time)

azdeg = round(azdeg / 0.05) * 0.05
print "Search AZ: " + str(azdeg)

print "Getting observation metadata"
metadata = get_obs_metadata(int(obsid))

print "Loading file"
theta, phi, amp = np.genfromtxt(args.file_loc, usecols=(0,1,8), skip_header=14, unpack=True)

#print "Sorting"
#az, za = np.meshgrid(np.radians(sorted(set(phi))), np.radians(sorted(set(theta))))
#delays = [metadata['xdelays'], metadata['ydelays']]
"""
print "Beam sim"
gx, gy = pb.MWA_Tile_full_EE(za, az, freq=freq*1e6, delays=delays, power=True, zenithnorm=True)
tile_beam = (gx + gy) / 2.0

amp /= amp.max()
amp *= np.ravel(tile_beam)
fill_min = 7e-3
fill_max = 0.95 * amp.max()
amp[amp <= fill_min] = 0
amp[amp >= fill_max] = fill_max
"""
max_value = np.max(amp)
local_max_index = np.where(amp==max_value)
max_phi = float(phi[local_max_index[0]])
print "Max phi in plot: " + str(max_phi)

print "Beam sort"
beam_line = []
beam_za = []
print "Amplitude_line: " + str(amp)

for i in range(len(amp)):
    if max_phi == float(phi[i]):
        beam_line.append(amp[i])
        beam_za.append(theta[i])
#print "Beam line to be plotted: " + str(beam_line)
spline = UnivariateSpline(beam_za, beam_line-np.max(beam_line)/2., s=0)
print spline.roots()
r1, r2 = spline.roots()
FWHM = abs(r1-r2)
print ra,dec,FWHM
print "Plot beam"

import matplotlib.pyplot as pl
pl.plot(beam_za, beam_line)
pl.savefig("zenith_line_plot_"+ra+"_"+dec+"_fwhm_"+str(FWHM)+".png")
pl.show()


