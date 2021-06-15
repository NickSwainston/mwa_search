#!/usr/bin/env python 

import numpy as np
import argparse
import math
from math import cos,sin
import glob
import sys
import csv

from mwa_search.obs_tools import calc_ta_fwhm
from vcstools.metadb_utils import get_common_obs_metadata, get_obs_array_phase
from vcstools.pointing_utils import sex2deg, deg2sex, format_ra_dec
from vcstools import prof_utils
from vcstools.gfit import gfit
from vcstools.prof_utils import NoFitError

import matplotlib.pyplot as plt
from matplotlib import patches


def find_pos(dec_search_range, ra_search_range, detections, fwhm, given_fwhm_ra=None, given_fwhm_dec=None, initial_pos=None):
    print("Localising with: (ra, dec, sn)")
    print(detections)
    RA = []; DEC = []; residual = []; psf_guass_grid = []
    for dec in dec_search_range:
        for ra in ra_search_range:
            # Loop over sources
            det_gauss = []
            for det in detections:
                if given_fwhm_ra is None:
                    # Adjust the FWHM for the projected change
                    fwhm_ra  = np.degrees(np.radians(fwhm)/cos(np.radians(dec + 26.7))**2)
                else:
                    fwhm_ra = given_fwhm_ra
                if given_fwhm_dec is None:
                    fwhm_dec = np.degrees(np.radians(fwhm)/cos(np.radians(dec)))
                else:
                    fwhm_dec = given_fwhm_dec
                # Calculate the guassian response
                ra_guass =  (ra  - det[0]) / (0.6006*fwhm_ra)
                dec_guass = (dec - det[1]) / (0.6006*fwhm_dec)
                det_gauss.append(math.exp(-(math.sqrt(ra_guass**2 + dec_guass**2))))
                #det_gauss.append(math.exp(-((ra  - det[0])**2 + (dec - det[1])**2) / (0.6006*fwhm)**2))

            """
            # Calculate a psf for normalisation
            ra_guass =  (ra  - ra_centre)  / (1.5*ra_search_diameter)
            dec_guass = (dec - dec_centre) / (1.5*dec_search_diameter)
            psf_gauss = math.exp(-(math.sqrt(ra_guass**2 + dec_guass**2)))
            psf_guass_grid.append(psf_gauss)
            """

            # compare SN ratios to guassian ratios for each beam pair
            res_sum = 0.
            for i in range(len(det_gauss)):
                for j in range(len(det_gauss)):
                    # This if makes sure we don't do repeat pairs
                    if i < j:
                        ras = [detections[i][0], detections[j][0]]
                        ras.sort()
                        decs = [detections[i][1], detections[j][1]]
                        decs.sort()
                        # if there is an initial pos make sure it is between the two beams
                        if initial_pos is None:
                            res_sum += ( (detections[i][2]/detections[j][2] - det_gauss[i]/det_gauss[j]) )**2
                        else:
                            if (ras[0]  < initial_pos[0] < ras[1]) or\
                               (decs[0] < initial_pos[1] < decs[1]):
                                res_sum += ( (detections[i][2]/detections[j][2] - det_gauss[i]/det_gauss[j]) )**2
            """
            for i in range(len(det_gauss)):
                res_sum += ( (detections[i_sn_max][2]/detections[i][2] - det_gauss[i_sn_max]/det_gauss[i])  )**2#* ((detections[i][2] + detections[j][2]) / 2*sn_max)**2
            """
            residual.append(math.sqrt(res_sum))
            RA.append(ra)
            DEC.append(dec)
    return RA, DEC, residual


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    Calculate the best position of a source from the singal to noise of several detections.
    """)
    parser.add_argument('-o', '--obsid', type=str,
            help='The observation ID of the fits file to be searched')
    parser.add_argument('-O', '--calid', type=str,
            help='The calibration ID of the fits file to be searched')
    parser.add_argument('-b', '--bestprof_dir', type=str,
            help='The directory of bestprof files of detections.')
    parser.add_argument('-p', '--pdmp_dir', type=str,
            help='The directory of pdmp files of detections.')
    parser.add_argument('-r', '--res', type=float, default=None,
            help='The resolution of the search in degrees.')
    parser.add_argument('-w', '--write', action='store_true',
            help='Write out a file with the predicted poistion.')
    parser.add_argument('-fr', '--fwhm_ra', type=float,
            help='Manualy give the RA FWHM in degrees instead of it estimating it from the array phase and frequency')
    parser.add_argument('-fd', '--fwhm_dec', type=float,
            help='Manualy give the declination FWHM in degrees instead of it estimating it from the array phase and frequency')
    parser.add_argument('--orig_pointing', nargs='+', type=str,
            help='The original pointing. If used will output a file with the original pointing and best pointing SNs.')
    parser.add_argument('--label', type=str,
            help='Label the output predicted position.')
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
    rads  = []
    decds = []
    sns   = []
    if args.bestprof_dir:
        for bestprof_file in glob.glob("{}/*bestprof".format(args.bestprof_dir)):
            with open(bestprof_file,"r") as bestprof:
                lines = bestprof.readlines()
                ra, dec = lines[0].split("=")[-1].split("_")[1:3]
                sn = float(lines[13].split("~")[-1].split(" ")[0])
                """
                if sn > 3.:
                    bestprof_data = prof_utils.get_from_bestprof(bestprof_file)
                    _, _, _, period, _, _, _, profile, _ = bestprof_data
                    g_fitter = gfit(profile)
                    try:
                        g_fitter.auto_gfit()
                        sn = g_fitter.fit_dict["sn"]
                    except NoFitError:
                        sn = 1
                    #try:
                    #    sn, sn_e, _ = prof_utils.est_sn_from_prof(profile, period, alpha=2.5)
                    #except:
                    #    sn = 1
                    #if sn is None:
                    #    sn = 1
                """
                rad, decd = sex2deg(ra, dec)
                detections.append([rad, decd, sn, ra, dec])
    elif args.pdmp_dir:
        for pdmp_file in glob.glob("{}/*posn".format(args.pdmp_dir)):
            with open(pdmp_file,"r") as pdmp:
                lines = pdmp.readlines()
                sn = float(lines[0].split()[3])
                ra, dec = lines[0].split()[9].split("_")[1:3]
                #dec = dec[:-3]
                if sn < 3.:
                    print("skipping RA: {}   Dec: {}  SN: {}".format(ra, dec, sn))
                else:
                    rad, decd = sex2deg(ra, dec)
                    detections.append([rad, decd, sn, ra, dec])
    else:
        print("Please either use --bestprof_dir or --pdmp_dir. Exiting.")
        sys.exit(1)

    # sort by SN
    rads  = []
    decds = []
    sns   = []
    detections.sort(key=lambda x: x[2], reverse=True)
    best_sn = detections[0][2]
    best_pointing = "{}_{}".format(detections[0][3], detections[0][4])
    print("Input detections:")
    for ra, dec, sn, raj, decj in detections:
        print("RA: {}  Dec: {}  SN: {}".format(raj, decj, sn))
        rads.append(ra)
        decds.append(dec)
        sns.append(sn)

        if args.orig_pointing:
            print("{}_{}".format(raj, decj))
            if "{}_{}".format(raj, decj) in args.orig_pointing:
                orig_sn = sn
                orig_pointing = "{}_{}".format(raj, decj)

    if args.orig_pointing:
        # output csv of orig and best SN
        with open("{}_{}_{}_orig_best_SN.txt".format(args.label, args.obsid, args.calid), "w") as outfile:
            spamwriter = csv.writer(outfile, delimiter=',')
            spamwriter.writerow([orig_pointing, orig_sn])
            spamwriter.writerow([best_pointing, best_sn])

    # remove raj and decj
    new_detections = []
    for d in detections:
        new_detections.append(d[:3])
    detections = new_detections
    # Find data max mins
    ra_max, dec_max, sn_max = np.max(detections, axis=0)
    #print("sn_max: {}".format(sn_max))
    ra_min, dec_min, sn_min = np.min(detections, axis=0)
    ra_search_diameter  = (ra_max  - ra_min)
    dec_search_diameter = (dec_max - dec_min)
    ra_centre  = (ra_max  + ra_min)  / 2
    dec_centre = (dec_max + dec_min) / 2
    for i in range(len(detections)):
        if detections[i][2] == sn_max:
            ra_sn_max, dec_sn_max, sn = detections[i]
            i_sn_max = i

    # Find the smallest distance between pointings
    dists = []
    for det1 in detections:
        ra1, dec1, _ = det1
        for det2 in detections:
            ra2, dec2, _ = det2
            if ra1 != ra2 and dec1 != dec2:
                dists.append(math.sqrt(math.pow(ra1 - ra2, 2) + math.pow(dec1 - dec2, 2)))
    min_dist = min(dists)

    
    # Only use detections closest to max
    mask   = np.ones(len(detections), dtype=np.bool)
    #mask[flux_data > 10.0] = 0
    radinside = 1.1*min_dist
    print("fit radius: {}".format(radinside))
    for i in range(len(detections)):
        #print("distance: {}".format(math.sqrt((rads[i] - ra_sn_max)**2. + (decds[i] - dec_sn_max)**2.)))
        if math.sqrt((rads[i] - ra_sn_max)**2. + (decds[i] - dec_sn_max)**2.) > radinside:
            mask[i] = 0
    
    #print("{}".format(np.array(detections)[mask,:]))
    #print(mask)
    # Mask manually because there are bugs
    centre_detections = []
    for i in range(len(detections)):
        if mask[i]:
            centre_detections.append(detections[i])
    centre_detections = np.array(centre_detections)

    # Make search area
    if args.res is None:
        res = fwhm / 60
    else:
        res = args.res  #in degrees

    ra_search_range  = np.arange(ra_min  - fwhm, ra_max  + fwhm, res)
    dec_search_range = np.arange(dec_min - fwhm, dec_max + fwhm, res)
    ax.axis([min(ra_search_range), max(ra_search_range), min(dec_search_range), max(dec_search_range)])

    """
    # Find initial estimate using top 3 SN
    RA, DEC, residual = find_pos(dec_search_range, ra_search_range, detections[:3], fwhm)
    ra_initial = RA[residual.index(min(residual))]
    dec_initial = DEC[residual.index(min(residual))]

    # Run again with all detections, only process if the initial detction is between the two beams
    RA, DEC, residual = find_pos(dec_search_range, ra_search_range, detections[:3], fwhm,
                                 given_fwhm_ra=args.fwhm_ra, given_fwhm_dec=args.fwhm_dec,
                                 initial_pos=[ra_initial, dec_initial])
    """

    RA, DEC, residual = find_pos(dec_search_range, ra_search_range, centre_detections, fwhm,
                                 given_fwhm_ra=args.fwhm_ra, given_fwhm_dec=args.fwhm_dec)

    ramax = RA[residual.index(min(residual))]
    decmax = DEC[residual.index(min(residual))]
    plt.scatter(ramax,decmax,s=3,c='red', zorder=10)

    rah, dech = deg2sex(ramax, decmax)
    rah, dech = format_ra_dec([[rah, dech]], ra_col = 0, dec_col = 1)[0]
    print("Predicted RA:  {} deg  Dec: {} deg".format(round(ramax, 4), round(decmax, 4)))
    print("Predicted pos: {}_{} ".format(rah, dech))

    if args.write:
        with open("predicted_pos.txt","w") as write_file:
            if args.label:
                write_file.write("{},{}_{} ".format(args.label, rah, dech))
            else:
                write_file.write("{}_{} ".format(rah, dech))

    # Calculated predicted SN
    fwhm_ra  = np.degrees(np.radians(fwhm)/cos(np.radians(dec_sn_max + 26.7))**2)
    fwhm_dec = np.degrees(np.radians(fwhm)/cos(np.radians(dec_sn_max)))
    # Calculate the guassian response
    ra_guass =  (ra_sn_max  - ramax)  / (0.6006*fwhm_ra)
    dec_guass = (dec_sn_max - decmax) / (0.6006*fwhm_dec)
    gauss = math.exp(-(math.sqrt(ra_guass**2 + dec_guass**2)))
    predicted_sn = sn_max / gauss
    print("Predicted SN:  {}".format(round(predicted_sn, 1)))


    nx = np.array(RA)
    ny = np.array(DEC)
    nz = np.array(residual)
    #nz = np.array(psf_guass_grid)

    for det in np.array(detections):
        ra, dec, sn = det
        fwhm_ra  = np.degrees(np.radians(fwhm)/cos(np.radians(dec + 26.7))**2)
        fwhm_dec = np.degrees(np.radians(fwhm)/cos(np.radians(dec)) )

        ellipse = patches.Ellipse((ra, dec), fwhm_ra, fwhm_dec,
                                   linewidth=0.3, fill=False, edgecolor='green')
        ax.add_patch(ellipse)
        plt.scatter(ra,dec,s=0.5,c='white', zorder=10)
        #print(det, np.array(detections)[~mask])
        if det[2] in (np.array(sns)[~mask]):
            ax.text(ra, dec, "{:.2f}".format(sn), fontsize=8, ha='center', va='center', color='0.5')
        else:
            ax.text(ra, dec, "{:.2f}".format(sn), fontsize=8, ha='center', va='center', color='0')

    #Start plotting
    colour_map = 'plasma_r'
    nx.shape = (len(dec_search_range),len(ra_search_range))
    ny.shape = (len(dec_search_range),len(ra_search_range))
    nz.shape = (len(dec_search_range),len(ra_search_range))
    if np.amax(nz) > 10*np.amin(nz):
        max_plot = 10*np.amin(nz)
    else:
        max_plot = np.amax(nz)
    max_plot = np.amax(nz)
    plt.pcolor(nx, ny, nz, cmap='plasma_r', vmin=np.amin(nz), vmax=max_plot, zorder=0.5, snap=True)
    plt.colorbar(spacing='uniform', label=r"Residual")
    plt.xlabel("Right Ascension")
    plt.ylabel("Declination")
    fig.savefig("{}_{}_{}_residual.png".format(args.label, args.obsid, args.calid), dpi=1000, bbox_inches='tight')
