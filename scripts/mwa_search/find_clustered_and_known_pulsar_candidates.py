import csv
import numpy as np
from scipy.stats import gaussian_kde
import pandas as pd
import argparse
import glob
import os
import shutil
import sys
from psrqpy import QueryATNF

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import ScalarFormatter, LogFormatter

import astropy.units as u
from astropy.coordinates import SkyCoord, search_around_sky


def find_clustered_cands(cand_data,
                         mdist=18.56,
                         mperiod=5.,
                         mdm=5.,
                         no_plot=False,
                         sn_min=5.):
    """
    Find candidate that are cluster in period, DM and position

    Parameters
    ----------
    cand_data: list
        List of lists in the format
        [RA, Dec, Period(ms), DM, SN, file_loc]
    mdist: float
        Distance between candidates to consider them clustered.
        Default 18.56 arcminutes (1 beamwidth)
    mperiod: float
        Maximum period (in ms) difference between candidates to consider them clustered.
        Default 1 ms
    mdm: float
        Maximum DM difference between candidates to consider them clustered.
        Default 5

    Returns
    -------
    clustered_cands: list
        Lost of the file locations of the clustered cands
    """
    # Transpose list for easier use
    tlist = list(zip(*cand_data))
    ra  = list(tlist[0])
    dec = list(tlist[1])

    clustered_cands = []
    #print("Making coords")
    coords = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle,u.deg))
    #print("Searching")
    # Check if they're within a beam width
    idx1, idx2, sep2d, _ = search_around_sky(coords, coords, mdist*u.arcminute)

    # Check if they have the similar period and DM
    period = []
    dm = []
    sn = []
    for i in range(len(sep2d)):
        cand_1 = cand_data[idx1[i]]
        cand_2 = cand_data[idx2[i]]
        if ( idx1[i] != idx2[i]) and \
            ( abs(float(cand_1[2]) - float(cand_2[2])) < mperiod ) and \
            ( abs(float(cand_1[3]) - float(cand_2[3])) < mdm ) and \
            ( cand_1[4] > sn_min or cand_2[4] > sn_min ):
            # and \
            #not ( 1642. < float(cand_1[2]) < 1644. ):
            #print("")
            #print(cand_1)
            #print(cand_2)
            period.append(float(cand_1[2]))
            period.append(float(cand_2[2]))
            dm.append(float(cand_1[3]))
            dm.append(float(cand_2[3]))
            sn.append(float(cand_1[4]))
            sn.append(float(cand_2[4]))
            clustered_cands.append(cand_1[5][:-8]+"png")
            clustered_cands.append(cand_2[5][:-8]+"png")
    if len(clustered_cands) == 0:
        print("No clustered candidates found")
        sys.exit(0)
    clustered_cands = pd.unique(clustered_cands).tolist()

    if not no_plot:
        #sort by SN
        period = np.array(period)
        dm = np.array(dm)
        sn = np.array(sn)
        idx = sn.argsort()
        dm, period, sn = dm[idx], period[idx], sn[idx]

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        plt.scatter(dm, period, s=10, c=sn, norm=LogNorm())
        plt.title("{} Clustered Candidates".format(len(clustered_cands)))
        ax.set_yscale('log')
        ax.get_yaxis().set_major_formatter(ScalarFormatter())
        ax.get_yaxis().set_minor_formatter(ScalarFormatter())
        ax.set_xlabel("DM")
        ax.set_ylabel("Period (ms)")
        formatter = LogFormatter(10, labelOnlyBase=False)
        cb = plt.colorbar(label=r"$\sim$SN (PRESTO sigma)", format=formatter)
        cb.ax.minorticks_on()
        plt.show()
        #plt.savefig("candidate_clustering_all.png")
    return clustered_cands


def find_known_pulsars(cand_data,
                       mdist=30,
                       mperiod=1.,
                       mdm=1.):
    """
    Find candidate that are cluster in period, DM and position

    Parameters
    ----------
    cand_data: list
        List of lists in the format
        [RA, Dec, Period(ms), DM, SN, file_loc]
    mdist: float
        Distance between candidates to consider them clustered.
        Default 18.56 arcminutes (1 beamwidth)
    mperiod: float
        Maximum period (in ms) difference between candidates to consider them clustered.
        Default 1 ms
    mdm: float
        Maximum DM difference between candidates to consider them clustered.
        Default 5

    Returns
    -------
    clustered_cands: list
        Lost of the file locations of the clustered cands
    """
    pminf = (100 - mperiod) / 100
    pmaxf = (100 + mperiod) / 100

    query = QueryATNF().pandas
    #query.save('full_query.pkl')
    coord_kp = SkyCoord(ra=[query["RAJ"]], dec=[query["DECJ"]], unit=(u.hourangle,u.deg))

    # possible harmonics
    harms = list(range(2,32))
    # possible number of peaks in profile
    possible_peaks = list(range(1,12))
    possible_fractions = []
    for h in harms:
        for p in possible_peaks:
            possible_fractions.append(h/p)
    #print(possible_fractions)

    known_pulsars = []
    no_known_pulsars = []
    for cand in cand_data:
        RA, Dec, period, DM, SN, file_loc = cand
        period = period / 1000. # convert to s
        kp = False

        # Check if RFI
        if (period*pminf < 1. < period*pmaxf) or (period*pminf < 2. < period*pmaxf):
            print("{:<10}          {}".format("RFI", file_loc))
            continue

        # Do a simple check first (not a harmonic)
        """
        print('(DM > {} && DM < {}) && (P0 > {} && P0 < {})'.format(DM-mdm, DM+mdm, period*pminf, period*pmaxf))
        kpquery = QueryATNF(condition='(DM > {} && DM < {}) && (P0 > {} && P0 < {})'.format(DM-mdm, DM+mdm, period*pminf, period*pmaxf)
                            , params=['PSRJ', 'P0', 'DM'], cache=True
                            , loadfromdb=data_load.ATNF_LOC)
                            #, loadquery='full_query.pkl')
        print(kpquery.pandas)
        """
        for i in range(len(query)):
            jname = query["PSRJ"][i]
            kp0 = query["P0"][i]
            kdm = query["DM"][i]
            if (DM-mdm < kdm < DM+mdm) and (period*pminf < kp0 < period*pmaxf):
                known_pulsars.append([jname, file_loc, None])
                kp = True
                print("{:<10}  {:6.3f}  {}".format(jname, 1.0, file_loc))
                break
        if kp:
            continue

        # Check if it's a harmoinc
        coord_cand = SkyCoord(ra=RA, dec=Dec, unit=(u.hourangle,u.deg))
        for i in range(len(query)):
            jname = query["PSRJ"][i]
            kp0 = query["P0"][i]
            kdm = query["DM"][i]

            # Find pulsars within distance
            coord_kp = SkyCoord(ra=[query["RAJ"][i]], dec=[query["DECJ"][i]], unit=(u.hourangle,u.deg))
            if coord_kp.separation(coord_cand) < mdist*u.deg:
                cand_harm = kp0 / period
                # Check the pulsars that are within the distance if they have a fraction of the candidate period
                for pf in possible_fractions:
                    if abs(cand_harm - pf) < 0.01 and (DM-mdm < kdm < DM+mdm):
                        known_pulsars.append([jname, file_loc, cand_harm])
                        print("{:<10}  {:6.3f}  {}".format(jname, cand_harm, file_loc))
                        kp = True
                        break
            if kp:
                break
        if kp:
            continue

        no_known_pulsars.append(file_loc)
        print(file_loc)
    return known_pulsars, no_known_pulsars


def get_from_bestprof(file_loc):
    """
    Get info from a bestprof file

    Parameters
    ----------
    file_loc: string
        The path to the bestprof file

    Returns
    -------
    [ra, dec, period, dm, sigma, file_loc]: list
    """

    with open(file_loc,"r") as bestprof:
        lines = bestprof.readlines()

        # Get a pointing from the input fits file name
        ra, dec = lines[0].split("_ch")[0].split("_")[-2:]

        # Get period and DM
        period = lines[15][22:-1]
        period = float(period.split(' +/- ')[0])
        dm = float(lines[14][22:-1])

        # Get sigma
        sigma = float(lines[13].split("~")[-1].split(" sigma")[0])


    return [ra, dec, period, dm, sigma, file_loc]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    A script designed to organise the SMART candidates and look for clustering in period and DM in adjacent beams
    """)
    cand_input = parser.add_argument_group('Candidate input')
    cand_input.add_argument('-o', '--obsid', type=str,
                           help='Obsid to be used on Prometheus. Will find the directory, do a bit of sorting and read the bestprofs')
    cand_input.add_argument('-g', '--glob', type=str,# nargs='+',
                           help='A wilcard expression used to find all the bestprofs files. Should be surround by double quotes.')
    cand_input.add_argument('-c', '--cand_file', type=str,
                           help='Candidate input file (generated by Sammy).')

    clust_par = parser.add_argument_group('Clustering Parameters')
    clust_par.add_argument('-cd', '--cdist', type=float, default=18.56,
                           help='Distance between candidates to consider them clustered. Default 18.56 arcminutes (1 beamwidth)')
    clust_par.add_argument('-cp', '--cperiod', type=float, default=5.,
                           help='Maximum period (in ms) difference between candidates to consider them clustered. Default 1 ms')
    clust_par.add_argument('--cdm', type=float, default=5.,
                           help='Maximum DM difference between candidates to consider them clustered. Default 5')
    clust_par.add_argument('--csn', type=float, default=5.,
                           help='SN cut off. Default 5 (basially no cut off)')

    clust_par = parser.add_argument_group('Search params')
    clust_par.add_argument('-sd', '--sdist', type=float, default=30,
                           help='Radius to search in degrees. Default 30 degrees')
    clust_par.add_argument('-sp', '--speriod', type=float, default=1.,
                           help='Percentage difference in period to search. Default 1')
    clust_par.add_argument('--sdm', type=float, default=1.5,
                           help='Maximum DM difference to search. Default 1.5')

    out_arg = parser.add_argument_group('Outputs')
    out_arg.add_argument('--no_plot', action='store_true',
                         help="Don't plot the candidate clustering P and DM space")
    out_arg.add_argument('--out_file', type=str,
                         help="An output file containing the bestprof file location of each clustered candidates.")
    out_arg.add_argument('--out_dir', type=str,
                         help="An output directory containing all of the clustered candidates.")
    args=parser.parse_args()

    # Parse candidate inputs
    if sum(map(bool, [args.obsid, args.glob, args.cand_file])) > 1:
        print("Please only use either --obsid, --glob or --candfile")
        sys.exit(0)
    elif args.obsid:
        # Move all positive candidates into one file if needed
        pos_sub_dirs = glob.glob('/data/nswainston/SMART_cand_sorting/{}/S*/positive_detections'.format(args.obsid))
        if len(pos_sub_dirs) > 0:
            print("Moving all cands into one directory")
            if not os.path.isdir('/data/nswainston/SMART_cand_sorting/{}/positive_detections'.format(args.obsid)):
                os.mkdir('/data/nswainston/SMART_cand_sorting/{}/positive_detections'.format(args.obsid))
            pos_files = glob.glob('/data/nswainston/SMART_cand_sorting/{}/S*/positive_detections/*'.format(args.obsid))
            for pf in pos_files:
                shutil.move(pf, '/data/nswainston/SMART_cand_sorting/{}/positive_detections'.format(args.obsid))

            # Delete old directories
            pos_dirs = glob.glob('/data/nswainston/SMART_cand_sorting/{}/S*/positive_detections/'.format(args.obsid))
            for pd in pos_dirs:
                # Checking whether the folder is empty or not
                if len(os.listdir(pd)) == 0:
                    # Removing the file using the os.remove() method
                    os.rmdir(pd)

        # Read in all ML positively classified candidates
        bestprof_files = glob.glob('/data/nswainston/SMART_cand_sorting/{}/positive_detections/*bestprof'.format(args.obsid))
        if len(bestprof_files) == 0:
            print("No bestprof files found in /data/nswainston/SMART_cand_sorting/{}/positive_detections/*bestprof".format(args.obsid))
            sys.exit(0)
        cand_data = []
        for bf in bestprof_files:
            cand_data.append(get_from_bestprof(bf))
    elif args.glob:
        # Read in all wildcard candidates
        cand_data = []
        #for bf in args.glob:
        #print(glob.glob(args.glob))
        for bf in glob.glob(args.glob):
            cand_data.append(get_from_bestprof(bf))
    elif args.cand_file:
        with open(args.cand_file) as csvfile:
            spamreader = csv.reader(csvfile, delimiter=' ')
            cand_data = []
            for row in spamreader:
                ra, dec, period, dm, sigma, file_loc = row
                cand_data.append([ra, dec, float(period), float(dm), float(sigma), file_loc])
    else:
        print("Please use either --obsid, --glob or --candfile")
        sys.exit(0)


    # Find clustered canidates
    print("Finding clusted candidates")
    clustered_cands = find_clustered_cands(cand_data,
                                              mdist=args.cdist,
                                              mperiod=args.cperiod,
                                              mdm=args.cdm,
                                              no_plot=args.no_plot,
                                              sn_min=args.csn)
    print("Found {} clustered candidates".format(len(clustered_cands)))

    clustered_cands_data = []
    for cd in cand_data:
        if cd[-1] in clustered_cands:
            clustered_cands_data.append(cd)

    print("Finding known pulsars")
    clustered_known_pulsars, clustered_not_known_pulsars = find_known_pulsars(clustered_cands_data,
                                              mdist=args.sdist,
                                              mperiod=args.speriod,
                                              mdm=args.sdm)
    files_to_follow_up = clustered_not_known_pulsars
    print("Found {} clustered not known pulsar candidates".format(len(files_to_follow_up)))

    # Output files if requested
    if args.out_file:
        print("Clustered candidates are listed in {}".format(args.out_file))
        with open(args.out_file, 'w') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=',')
            for file_name in files_to_follow_up:
                spamwriter.writerow([file_name])
    if args.out_dir:
        print("Moving clustered candidates to {}".format(args.out_dir))
        if not os.path.isdir(args.out_dir):
            # Create path if it doesn't exist
            os.mkdir(args.out_dir)
        for file_name in files_to_follow_up:
            # grab all presto files
            for pf in glob.glob("{}*".format(file_name[:-4])):
                shutil.move(pf, args.out_dir)
