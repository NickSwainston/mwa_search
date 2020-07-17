#! /usr/bin/env python

import sys
import argparse
import os
import glob
import math
import numpy as np
import psrqpy
from astropy.coordinates import SkyCoord
import astropy.units as u
import csv

# vcstools and mwa_search imports
import find_pulsar_in_obs as fpio
import mwa_search_pipeline as search_pipe
from mwa_metadb_utils import get_common_obs_metadata as get_meta
from mwa_metadb_utils import obs_max_min, get_obs_array_phase
from check_known_pulsars import calc_ta_fwhm, get_pointings_required, get_sources_in_fov
from pulsar_obs_helper import find_pulsars_power
from config_vcs import load_config_file
from grid import get_grid
import checks
import sn_flux_est as snfe

import logging
logger = logging.getLogger(__name__)

#get ATNF db location
try:
    ATNF_LOC = os.environ['PSRCAT_FILE']
except KeyError:
    logger.warn("ATNF database could not be loaded on disk. This may lead to a connection failure")
    ATNF_LOC = None


    

def find_pulsars_in_fov(obsid, psrbeg, psrend):
    """
    Find all pulsars in the field of view and return all the pointings sorted into vdif and normal lists:
    -----------
    obsid: int
        The observation ID
    psrbeg: int
        The begining of the observation you are processing in GPS time
    psrend: int
        The end of the observation you are processing in GPS time

    Returns:
    --------
    list of lists:
           [pulsar_name_list,
            pulsar_pointing_list,
            vdif_name_list,
            vdif_pointing_list,
            sp_name_list,
            sp_pointing_list]
    """
    
    #Find all pulsars in beam at at least 0.3 of zenith normlaized power
    names_ra_dec = np.array(fpio.grab_source_alog(max_dm=250))
    pow_dict, meta_data = find_pulsars_power(obsid, powers=[0.3, 0.1], names_ra_dec=names_ra_dec)
    channels = meta_data[-1][-1]
    obs_psrs = pow_dict[0.3][obsid]
    psrs_list_03 = [x[0] for x in obs_psrs]
    #Include all bright pulsars in beam at at least 0.1 of zenith normalized power
    psrs_01 = [x[0] for x in pow_dict[0.1][obsid]]
    sn_dict_01 = snfe.multi_psr_snfe(psrs_01, obsid, beg=psrbeg, end=psrend, min_z_power=0.1)
    for psr in pow_dict[0.1][obsid]:
        if psr[0] not in psrs_list_03:
            sn, sn_err = sn_dict_01[psr[0]]
            if sn is not None and sn_err is not None:
                if sn - sn_err >= 10.:
                    obs_psrs.append(psr)

    #get all the pulsars periods
    pulsar_list = []
    for o in obs_psrs:
        pulsar_list.append(o[0])
    periods = psrqpy.QueryATNF(params=["P0"], psrs=pulsar_list,
                               loadfromdb=ATNF_LOC).pandas["P0"]


    oap = get_obs_array_phase(obsid)
    centrefreq = 1.28 * float(min(channels) + max(channels)) / 2.
    fwhm = calc_ta_fwhm(centrefreq, array_phase=oap)

    # Sort all the sources into 3 categories, pulsars which is for slow pulsars, vdif
    # for fast pulsars that require vdif and sp for singple pulse searches (FRBs,
    # RRATs and pulsars without ATNF periods)
    pulsar_pointing_list = []
    pulsar_name_list = []
    vdif_pointing_list = []
    vdif_name_list = []
    pulsar_search_pointing_list = []
    pulsar_search_name_list = []
    sp_pointing_list = []
    sp_name_list = []
    
    for pi, pulsar_line in enumerate(obs_psrs):
        vdif_check = False
        sp_check = False

        PSRJ = pulsar_line[0]
        #See if pulsar is in beam for times
        #enter, exit = snfe.pulsar_beam_coverage(obsid, PSRJ, beg=psrbeg, end=psrend)
        #if enter is None or exit is None:
        #    print("{0} is not in beam for time range. Not beamforming".format(PSRJ))
        #    continue
        #TODO: uncomment this section when sn_flux_est is pulled to vcstools master

        if not (len(PSRJ) < 11 or PSRJ[-1] == 'A' or PSRJ[-2:] == 'aa'):
            continue

        for line in names_ra_dec:
            if PSRJ == line[0]:
                temp = [line]

        #temp = fpio.get_psrcat_ra_dec(pulsar_list=[PSRJ])
        temp = fpio.format_ra_dec(temp, ra_col = 1, dec_col = 2)
        jname, raj, decj = temp[0]
        #get pulsar period
        period = periods[pi]
        if math.isnan(period):
            print("WARNING: Period not found in ephermeris for {0} so assuming "
                  "it's an RRAT".format(jname))
            sp_check = True
            period = 0.
        elif float(period) < .05:
            vdif_check = True
        period = float(period)*1000.
        print("{0:12} RA: {1} Dec: {2} Period: {3:8.2f} (ms) Begin {4} End {5}".format(PSRJ,
               raj, decj, period, psrbeg, psrend))

        jname_temp_list = []
        if PSRJ[-1] == 'A' or PSRJ[-2:] == 'aa':
            #Got to find all the pulsar J names with other letters
            vdif_check = True
            for pulsar_l in obs_psrs:
                pulsar_name = pulsar_l[0]
                if pulsar_name.startswith(PSRJ[:-2]):
                    jname_temp_list.append(pulsar_name)
        else:
            jname_temp_list.append(jname)

        # grid the pointings to fill 2 arcminute raduis to account for ionosphere shift
        pointing_list_list = get_pointings_required(raj, decj, fwhm, 2./60.)

        # sort the pointings into the right groups
        for prd in pointing_list_list:
            if vdif_check:
                vdif_name_list.append(jname_temp_list)
                vdif_pointing_list.append("{0}_{1}".format(prd[0], prd[1]))
            elif sp_check:
                sp_name_list.append(jname_temp_list)
                sp_pointing_list.append("{0}_{1}".format(prd[0], prd[1]))
            else:
                pulsar_name_list.append(jname_temp_list)
                pulsar_pointing_list.append("{0}_{1}".format(prd[0], prd[1]))

    
    #Get the rest of the singple pulse search canidates
    #-----------------------------------------------------------------------------------------------------------
    temp = get_sources_in_fov(obsid, 'RRATs', fwhm)
    sp_name_list = sp_name_list + temp[0]
    sp_pointing_list = sp_pointing_list + temp[1]

    # Find all of the FRB candidates
    #-----------------------------------------------------------------------------------------------------------
    temp = get_sources_in_fov(obsid, 'FRB', fwhm)
    sp_name_list = sp_name_list + temp[0]
    sp_pointing_list = sp_pointing_list + temp[1]

    # Find all of the Fermi candidates
    #-----------------------------------------------------------------------------------------------------------
    fermi_list = get_sources_in_fov(obsid, 'Fermi', fwhm)
    pulsar_search_name_list = pulsar_search_name_list + fermi_list[0]
    pulsar_search_pointing_list = pulsar_search_pointing_list + fermi_list[1]

    # Find all of the points of interest candidates
    #-----------------------------------------------------------------------------------------------
    poi_list = get_sources_in_fov(obsid, 'POI', fwhm)
    pulsar_search_name_list = pulsar_search_name_list + poi_list[0]
    pulsar_search_pointing_list = pulsar_search_pointing_list + poi_list[1]

    # Sometimes we get redundant RRATs that are found in RRAT and ANTF catalogues so they need to be removed
    sp_name_list     = list(dict.fromkeys([ ";".join(s) for s in sp_name_list]))
    sp_pointing_list = list(dict.fromkeys(sp_pointing_list))

    # Changing the format of the names list to make it easier to format
    return [[ ";".join(s) for s in pulsar_name_list],
            pulsar_pointing_list,
            [ ";".join(s) for s in vdif_name_list],
            vdif_pointing_list,
            [ ";".join(s) for s in pulsar_search_name_list],
            pulsar_search_pointing_list,
            sp_name_list,
            sp_pointing_list]

if __name__ == "__main__":
    loglevels = dict(DEBUG=logging.DEBUG,
                     INFO=logging.INFO,
                     WARNING=logging.WARNING)
    parser = argparse.ArgumentParser(description="""
    Writes a file full of all the names and pointings of sources in the FoV.
    """)
    parser.add_argument('-o', '--obsid', type=str,
            help='The observation ID of the fits file to be searched')
    parser.add_argument("-b", "--begin", type=int, help="First GPS time to process [no default]")
    parser.add_argument("-e", "--end", type=int, help="Last GPS time to process [no default]")
    parser.add_argument("-a", "--all", action="store_true", default=False,
            help="Perform on entire observation span. Use instead of -b & -e")
    parser.add_argument("-L", "--loglvl", type=str, help="Logger verbosity level. Default: INFO",
                        default="INFO")
    args=parser.parse_args()

    # set up the logger for stand-alone execution
    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    ch.setLevel(loglevels[args.loglvl])
    formatter = logging.Formatter('%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.propagate = False

    #option parsing
    if not args.obsid:
        print("Please input observation id by setting -o or --obsid. Exiting")
        quit()

    if args.begin and args.end:
        beg = args.begin
        end = args.end
    elif args.all:
        beg, end = obs_max_min(args.obsid)

    output_list = find_pulsars_in_fov(args.obsid, beg, end)
    with open('{}_fov_sources_temp.csv'.format(args.obsid), 'w', newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',')
        for ol in output_list:
            if len(ol) == 0:
                # Print a space to the line which prevents Nextflow formatting erorrs
                spamwriter.writerow([" "])
            else:
                spamwriter.writerow(ol)

    with open('{}_fov_sources_temp.csv'.format(args.obsid), 'r') as readfile:
        csv_read = readfile.readlines()

    with open('{}_fov_sources.csv'.format(args.obsid), 'w') as csvfile:
        for line in csv_read:
            if len(line) == 0:
                csvfile.write(" \n")
            else:
                csvfile.write(line)