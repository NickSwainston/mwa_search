#!/usr/bin/env python 

import numpy as np
import argparse
import math
from math import cos,sin
import glob
import sys

from check_known_pulsars import calc_ta_fwhm
from mwa_metadb_utils import get_common_obs_metadata, get_obs_array_phase
from submit_to_database import sex2deg
from find_pulsar_in_obs import deg2sex, format_ra_dec

import matplotlib.pyplot as plt
from matplotlib import patches

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    Calculate the best position of a source from the singal to noise of several detections.
    """)
    parser.add_argument('-o', '--obsid', type=str,
            help='The observation ID of the fits file to be searched')
    parser.add_argument('-b', '--bestprof_dir', type=str,
            help='The directory of bestprof files of detections.')
    parser.add_argument('-p', '--pdmp_dir', type=str,
            help='The directory of pdmp files of detections.')
    parser.add_argument('-r', '--res', type=float, default=None,
            help='The resolution of the search in degrees.')
    parser.add_argument('-w', '--write', action='store_true',
            help='Write out a file with the predicted poistion.')
    args=parser.parse_args()

    # Set up plots
    fig = plt.figure()
    plt.rc("font", size=8)
    fig.add_subplot(111)
    ax = plt.axes()
    ax.axis('equal')

    # Get the fwhm of the observation
    meta_data = get_common_obs_metadata(args.obsid)
    channels = meta_data[-1]
    oap = get_obs_array_phase(args.obsid)
    if oap == "OTH":
        #Assume it's phase 2 extended array
        oap = "P2E"
    centrefreq = 1.28 * float(min(channels) + max(channels)) / 2.
    fwhm = calc_ta_fwhm(centrefreq, array_phase=oap)
    print("Observation ID: {}".format(args.obsid))
    print("FWHM: {} deg".format(fwhm))

    detections = []
    if args.bestprof_dir:
        for bestprof_file in glob.glob("{}/*bestprof".format(args.bestprof_dir)):
            with open(bestprof_file,"r") as bestprof:
                lines = bestprof.readlines()
                ra, dec = lines[0].split("=")[-1].split("_")[1:3]
                sn = float(lines[13].split("~")[-1].split(" ")[0])
                if sn < 3.:
                    print("skipping RA: {}   Dec: {}  SN: {}".format(ra, dec, sn))
                else:
                    print("RA: {}   Dec: {}  SN: {}".format(ra, dec, sn))
                    rad, decd = sex2deg(ra, dec)
                    detections.append([rad, decd, sn])
    elif args.pdmp_dir:
        for pdmp_file in glob.glob("{}/*posn".format(args.pdmp_dir)):
            with open(pdmp_file,"r") as pdmp:
                lines = pdmp.readlines()
                sn = float(lines[0].split()[3])
                ra, dec = lines[0].split()[9].split("_")[1:3]
                dec = dec[:-3]
                if sn < 3.:
                    print("skipping RA: {}   Dec: {}  SN: {}".format(ra, dec, sn))
                else:
                    print("RA: {}   Dec: {}  SN: {}".format(ra, dec, sn))
                    rad, decd = sex2deg(ra, dec)
                    detections.append([rad, decd, sn])
    else:
        print("Please either use --bestprof_dir or --pdmp_dir. Exiting.")
        sys.exit(1)

    # Find data max mins
    ra_max, dec_max, sn_max = np.max(detections, axis=0)
    ra_min, dec_min, sn_min = np.min(detections, axis=0)

    # Make search area
    if args.res is None:
        res = fwhm / 60
    else:
        res = args.res  #in degrees

    ra_search_range  = np.arange(ra_min  - fwhm, ra_max  + fwhm, res)
    dec_search_range = np.arange(dec_min - fwhm, dec_max + fwhm, res)
    ax.axis([min(ra_search_range), max(ra_search_range), min(dec_search_range), max(dec_search_range)])

    # Looping over search sky area
    RA = []; DEC = []; residual = []
    for dec in dec_search_range:
        for ra in ra_search_range:
            # Loop over sources
            det_gauss = []
            for det in detections:
                # Adjust the FWHM for the projected change
                fwhm_ra  = np.degrees(np.radians(fwhm)/cos(np.radians(dec + 26.7))**2)
                fwhm_dec = np.degrees(np.radians(fwhm)/cos(np.radians(dec)))
                # Calculate the guassian response
                ra_guass =  (ra  - det[0]) / (0.6006*fwhm_ra)
                dec_guass = (dec - det[1]) / (0.6006*fwhm_dec)
                det_gauss.append(math.exp(-(math.sqrt(ra_guass**2 + dec_guass**2))))
                #det_gauss.append(math.exp(-((ra  - det[0])**2 + (dec - det[1])**2) / (0.6006*fwhm)**2))
            # compare SN ratios to guassian ratios for each beam pair
            res_sum = 0.
            for i in range(len(det_gauss)):
                for j in range(len(det_gauss)):
                    if i < j:
                        # This if makes sure we don't do repeat pairs
                        res_sum += (detections[i][2]/detections[j][2]-det_gauss[i]/det_gauss[j])**2
            residual.append(math.sqrt(res_sum))
            RA.append(ra)
            DEC.append(dec)

    ramax = RA[residual.index(min(residual))]
    decmax = DEC[residual.index(min(residual))]
    plt.scatter(ramax,decmax,s=3,c='red', zorder=10)

    rah, dech = deg2sex(ramax, decmax)
    rah, dech = format_ra_dec([[rah, dech]], ra_col = 0, dec_col = 1)[0]
    print("Predicted RA:  {} deg  Dec: {} deg".format(round(ramax, 4), round(decmax, 4)))
    print("Predicted pos: {}_{} ".format(rah, dech))

    if args.write:
        with open("predicted_pos.txt","w") as write_file:
            write_file.write("{}_{} ".format(rah, dech))

    # Calculated predicted SN
    for det in detections:
        if det[2] == sn_max:
            ra, dec, sn = det
            fwhm_ra  = np.degrees(np.radians(fwhm)/cos(np.radians(dec + 26.7))**2)
            fwhm_dec = np.degrees(np.radians(fwhm)/cos(np.radians(dec)))
            # Calculate the guassian response
            ra_guass =  (ra  - ramax)  / (0.6006*fwhm_ra)
            dec_guass = (dec - decmax) / (0.6006*fwhm_dec)
            gauss = math.exp(-(math.sqrt(ra_guass**2 + dec_guass**2)))
            predicted_sn = sn / gauss
    print("Predicted SN:  {}".format(round(predicted_sn, 1)))




    nx = np.array(RA)
    ny = np.array(DEC)
    nz = np.array(residual)

    for det in detections:
        ra, dec, sn = det
        fwhm_ra  = np.degrees(np.radians(fwhm)/cos(np.radians(dec + 26.7))**2)
        fwhm_dec = np.degrees(np.radians(fwhm)/cos(np.radians(dec)) )

        ellipse = patches.Ellipse((ra, dec), fwhm_ra, fwhm_dec,
                                   linewidth=0.3, fill=False, edgecolor='green')
        ax.add_patch(ellipse)
        plt.scatter(ra,dec,s=0.5,c='white', zorder=10)
        ax.text(ra, dec, str(sn), fontsize=12, ha='center', va='center')

    #Start plotting
    colour_map = 'plasma_r'
    nx.shape = (len(dec_search_range),len(ra_search_range))
    ny.shape = (len(dec_search_range),len(ra_search_range))
    nz.shape = (len(dec_search_range),len(ra_search_range))
    if max(residual) > 10*min(residual):
        max_plot = 10*min(residual)
    else:
        max_plot = max(residual)
    plt.pcolor(nx, ny, nz, cmap='plasma_r', vmin=min(residual), vmax=max_plot, zorder=0.5, snap=True)
    plt.colorbar(spacing='uniform', label=r"Residual")
    plt.xlabel("Right Ascension")
    plt.ylabel("Declination")
    fig.savefig("residual.png", dpi=1000, bbox_inches='tight')