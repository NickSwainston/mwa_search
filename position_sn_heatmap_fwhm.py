#! /usr/bin/env python

import argparse
import os

from astropy.coordinates import SkyCoord
import astropy.units as u

from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from matplotlib import cm as CM
from matplotlib import mlab as ML

import numpy as np
import math
from scipy.interpolate import UnivariateSpline



def find_fwhm_and_plot(obsid, pointing):
    pointing_list = []
    sn = []
    for d in os.listdir("/group/mwaops/vcs/{0}/pointings/".format(obsid)):
        if os.path.exists("/group/mwaops/vcs/{0}/pointings/{1}/1221832280_fwhm_test_{1}_PSR_J2145-0750.pfd.bestprof".format(obsid, d)):
            with open("/group/mwaops/vcs/{0}/pointings/{1}/1221832280_fwhm_test_{1}_PSR_J2145-0750.pfd.bestprof".format(obsid, d)) as bestfile:
                for i in bestfile.readlines():
                    if i.startswith("# Prob(Noise)"):
                        pointing_list.append(d)
                        sn.append(float(i.split("(")[-1][1:-7]))

    #find max for a FWHM test
    #max_index = sn.index(max(sn))
    ra_hex = pointing.split("_")[0]
    dec_hex = pointing.split("_")[1]

    coord = SkyCoord(ra_hex,dec_hex,unit=(u.hourangle,u.deg))
    dec_centre = coord.dec.degree
    ra_centre = coord.ra.degree

    ras = []; decs = []
    ra_line = []; ra_sn_line = []
    dec_line = []; dec_sn_line = []
    for i in range(len(sn)):
        rah, dech = pointing_list[i].split("_")
        coord = SkyCoord(rah,dech,unit=(u.hourangle,u.deg))
        ras.append(coord.ra.degree)
        decs.append(coord.dec.degree)
        #print decs,dec_centre
        if decs[i] == dec_centre:
            ra_line.append(ras[i])
            ra_sn_line.append(sn[i])
        if ras[i] == ra_centre:
            dec_line.append(decs[i])
            dec_sn_line.append(sn[i])

    #sort and calc FWHM
    ra_sn_line = [x for _,x in sorted(zip(ra_line,ra_sn_line))]
    ra_line = sorted(ra_line)
    #print ra_sn_line,ra_line

    dec_sn_line = [x for _,x in sorted(zip(dec_line,dec_sn_line))]
    dec_line = sorted(dec_line)
    #print dec_sn_line,dec_line 

    ra_sn_line = np.array(ra_sn_line) ; dec_sn_line = np.array(dec_sn_line)

    spline = UnivariateSpline(ra_line, ra_sn_line-np.max(sn)/2., s=0)
    #print spline.roots()
    r1, r2 = spline.roots()
    ra_FWHM = abs(r1-r2)
    print "raw ra FHWM: " + str(ra_FWHM)
    cor_ra_FWHM = float(ra_FWHM)*math.cos(np.radians(dec_centre))
    print "corrected ra FWHM: {0}".format(cor_ra_FWHM)
    max_ra_i = np.argmax(ra_sn_line)
    print "max sn at ra: {0}".format(ra_line[max_ra_i])

    spline = UnivariateSpline(dec_line, dec_sn_line-np.max(sn)/2., s=0)
    #print spline.roots()
    r1, r2 = spline.roots()
    dec_FWHM = abs(r1-r2)
    print "raw dec FHWM: " + str(dec_FWHM)
    cor_dec_FWHM = float(dec_FWHM)*math.cos(np.radians(dec_centre) + np.radians(26.7))**2
    print "corrected dec FWHM: {0}".format(cor_dec_FWHM)
    max_dec_i = np.argmax(dec_sn_line)
    print "max sn at dec: {0}".format(dec_line[max_dec_i])

    max_coord = SkyCoord(ra_line[max_ra_i], dec_line[max_dec_i],unit=(u.deg,u.deg))
    ra_max_hex = max_coord.ra.to_string(unit=u.hour, sep=':')
    dec_max_hex = max_coord.dec.to_string(unit=u.degree, sep=':')

    print "sn max coord: {0}_{1}".format(ra_max_hex, dec_max_hex)




    ras = np.array(ras); decs = np.array(decs)
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111)

    plt.grid(True)

    cm = plt.cm.get_cmap('cubehelix',20)
    ax.grid(color='r', linestyle='-', linewidth=2)
    sp = plt.scatter(ras, decs, c=sn, s=500, edgecolors='b',cmap = cm)
    #plt.gray()
    plt.colorbar(sp)
    plt.savefig("{0}_position_heatmap.png".format(obsid))
    plt.show()
    #print sorted(ras)
    #print sorted(decs)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""
             For a cross pattern of pointings (probably done by grid.py -m cross)
             that have been prepfolded, uses their signal to noise to work out the fwhm.
             """)
    parser.add_argument('-o','--obsid',type=str,help='The observation ID to be searched')
    parser.add_argument('-p','--pointing',type=str,help='The centre pointing. eg 21:45:50.46_-07:50:18.48.')
    args=parser.parse_args()
    find_fwhm_and_plot(args.obsid, args.pointing) 