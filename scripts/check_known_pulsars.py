#! /usr/bin/env python3

import sys
import argparse
import os
import glob
import math
import numpy as np
import psrqpy
from astropy.coordinates import SkyCoord
import astropy.units as u
import config

# vcstools and mwa_search imports
import find_pulsar_in_obs as fpio
import mwa_search_pipeline as search_pipe
from mwa_metadb_utils import get_common_obs_metadata as get_meta
from mwa_metadb_utils import obs_max_min, get_obs_array_phase
from grid import get_grid
import checks
import sn_flux_est as snfe
from mwa_pulsar_client import client
import submit_to_database as std

import logging
logger = logging.getLogger(__name__)

#get ATNF db location
try:
    ATNF_LOC = os.environ['PSRCAT_FILE']
except KeyError:
    logger.warn("ATNF database could not be loaded on disk. This may lead to a connection failure")
    ATNF_LOC = None

def search_for_cal_srclist(obsid, cal_id, all_cal_returns=False, all_srclist_returns=False):
    """
    Given an obsid, searches common locations for the rts folder(s) as well as the sourcelist(s)

    Parameters:
    -----------
    obsid: int
        The observation ID
    cal_id: int
        The calibrator ID
    all_cal_returns: boolean
        OPTIONAL - If true, will return all RTS directories found. Default: False
    all_srclist_returns: boolean
        OPTIONAL - If true, will return all sourcelist files found. Default: False

    Returns:
    --------
    cal_path: string
        The location of the calibrator path
    srclist: string
        The pathname of the sourcelist
    """
    comp_config = config.load_config_file()
    base_dir = comp_config['base_product_dir']
    cal_dir = os.path.join(base_dir, str(obsid), "cal", str(cal_id))
    cal_dirs=[]
    srclists=[]
    #search all subdirectories
    for root, dirs, files in os.walk(cal_dir):
        if "rts" in dirs:
            cal_dirs.append(os.path.join(root, "rts"))
        for f in files:
            if f.endswith(".txt") and "srclist" in f:
                srclists.append(os.path.join(root, f))

    #handle multiple rts folders with user input
    if not all_cal_returns and len(cal_dirs) > 1:
        valid = False
        while not valid:
            print("Multiple RTS files found. Please choose one")
            for i, a_dir in enumerate(cal_dirs):
                print("{0}: {1}".format(i+1, a_dir))
            choice = int(input("Choose a number between 1 and {0}: ".format(len(cal_dirs))))
            if choice >= 1 and choice <= len(cal_dirs):
                valid=True
                my_cal_dir = cal_dirs[choice-1]
                print("Using RTS directory: {}".format(my_cal_dir))
                cal_dirs = [my_cal_dir]
            else:
                print("## Not a valid choice! ##")

    #handle multiple sourcelist files with user input
    if not all_srclist_returns and len(srclists) > 1:
        valid = False
        print("Multiple sourcelist files found. Please choose one")
        while not valid:
            for i, a_file in enumerate(srclists):
                print("{0}: {1}".format(i+1, a_file))
            choice = int(input("Choose a number between 1 and {0}: ".format(len(srclists))))
            if choice >= 1 and choice <= len(srclists):
                valid=True
                my_srclist = srclists[choice-1]
                print("Using sourcelist directory: {}".format(my_srclist))
                srclists = [my_srclist]
            else:
                print("## Not a valid choice! ##")

    return cal_dirs, srclists

def find_beg_end(obsid, base_path="/group/mwaops/vcs/"):
    """
    looks through the comined files of the obsid to find the beginning and end gps times

    Parameters:
    -----------
    obsid: int
        The observation ID
    base_path: string
        OPTIONAL - The system's base working path. Default = '/group/mwaops/vcs/'

    Returns:
    --------
    beg: int
        The beginning time for on-disk files
    end: int
        The end time for on-disk files
    """
    #TODO have some sort of check to look for gaps
    if glob.glob("{0}/{1}/combined/{1}*_ics.dat".format(base_path, obsid)):
        combined_files = glob.glob("{0}/{1}/combined/{1}*_ics.dat".format(base_path, obsid))
    else:
        meta_data = get_meta(obsid)
        channels = meta_data[-1]
        combined_files = glob.glob("{0}/{1}/combined/{1}*_ch{2}.dat".\
                                   format(base_path, obsid, channels[-1]))
    comb_times = []
    for comb in combined_files:
        comb_times.append(int(comb.split("_")[1]))
    beg = min(comb_times)
    end = max(comb_times)

    return beg, end

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
        comp_config = config.load_config_file()
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
    if array_phase == 'P2C':
        # True max_baseline is 700.
        max_baseline = 360.
    elif array_phase == 'P2E':
        max_baseline = 5300.

    wavelength = c / (freq * 1e6)
    fwhm = degrees(wavelength / max_baseline)

    return fwhm

def get_pointings_required(source_ra, source_dec, fwhm, search_radius):
    """
    Gets the number of grid pointings required to cover the search radius

    Parameters
    ----------
    source_ra, source_dec: string
        A string separated representing the RA and dec respectively.
        Expected format is 'hh:mm[:ss.s]'
    fwhm: float
        FWHM of the tied-array beam in degrees.
        Can be calculated in the calc_ta_fwhm function
    search_radius: float
        The radius of the circle that you would like to search

    Returns
    -------
    pointing_list_list: list of lists
        A list of pointings where each pointing contains an RA and a Dec in the format 'hh:mm:ss.ss'
        [[RA, Dec]]
    """
    #convert to radians
    coord = SkyCoord(source_ra, source_dec, unit=(u.hourangle,u.deg))
    rar = coord.ra.radian #in radians
    decr = coord.dec.radian

    #make a grid around each pulsar
    grid_sep = fwhm * 0.6
    #work out how many loops are required
    loops = int( (search_radius - fwhm/2.) / grid_sep )
    if loops < 0:
        loops = 0
    logger.debug("loops: {}".format(loops))
    rads, decds = get_grid(rar, decr, np.radians(grid_sep), loops)

    #convert back to sexidecimals
    coord = SkyCoord(rads,decds,unit=(u.deg,u.deg))
    rajs = coord.ra.to_string(unit=u.hour, sep=':')
    decjs = coord.dec.to_string(unit=u.degree, sep=':')
    temp = []
    for raj, decj in zip(rajs, decjs):
        temp.append([raj, decj])
    pointing_list_list = fpio.format_ra_dec(temp, ra_col = 0, dec_col = 1)
    return pointing_list_list

def find_pulsars_power(obsid, powers=None, names_ra_dec=None):
    """
    Finds the beam power information for pulsars in a specific obsid

    Parameters:
    -----------
    obsid: int
        The observation ID
    powers: list/tuple
        OPTIONAL - A list of minimum beam powers to evaluate the pulsar coverage at. If none, will use [0.3, 0.1]. Default: None
    names_ra_dec: list
        OPTIONAL - A list of puslars and their RA and Dec values to evaluate (generated from fpio.get_source_alog).
                   If none, will look for all pulsars. Default: None

    Returns:
    --------
    pulsar_power_dict: dictionary
        Contains keys - power
            Contains key - obsid
                Contains one list for each pulsar found in that power
                    Each list is constructed as [jname, enter, exit, max_power]
    meta_data: list
        A list of the output of get_common_obs_metadata for the input obsid
    """
    if not powers:
        powers = [0.3, 0.1]
    elif not (isinstance(powers, list) or isinstance(powers, tuple)):
        #try this if powers isn't iterable
        powers=list(powers)

    if names_ra_dec is None:
        names_ra_dec = np.array(fpio.grab_source_alog(max_dm=250))

    pulsar_power_dict = {}
    for pwr in powers:
        obs_data, meta_data = fpio.find_sources_in_obs([obsid], names_ra_dec, dt_input=100, min_power=pwr)
        pulsar_power_dict[pwr] = obs_data

    return pulsar_power_dict, meta_data

def get_sources_in_fov(obsid, source_type, fwhm):
    """
    Find all sources of the input type in the observations field-of-view

    Parameters:
    -----------
    obsid: str
        observation ID to search in
    source_type: str
        the source type input to fpio.grab_source_alog
    fwhm: float
        FWHM of the tied-array beam in degrees.
        Can be calculated in the calc_ta_fwhm function

    Returns:
    --------
    list:
        name_list: list
            A list of pulsars in the FOV
        pointing_list: list
            A list of pointings corresponding to the pulsars in name_list
    """
    names_ra_dec = fpio.grab_source_alog(source_type=source_type)
    obs_data, meta_data = fpio.find_sources_in_obs([obsid], names_ra_dec, dt_input=100)

    name_list = []
    pointing_list = []
    for pulsar_line in obs_data[obsid]:
        jname = pulsar_line[0]
        for line in names_ra_dec:
            if jname == line[0]:
                jname, raj, decj = line
        jname_temp_list = [jname]

        # grid the pointings to fill the position uncertaint (given in arcminutes)
        pointing_list_list = get_pointings_required(raj, decj, fwhm, 1./60.)

        # sort the pointings into the right groups
        for prd in pointing_list_list:
            name_list.append(jname_temp_list)
            pointing_list.append("{0}_{1}".format(prd[0], prd[1]))
    return [name_list, pointing_list]


def submit_folds(obsid, DI_dir, cal_obs, args, psrbeg, psrend,
                      product_dir='/group/mwaops/vcs',
                      mwa_search_version='master',
                      vcstools_version='master',
                      relaunch=False):
    """
    Beamforms on all pulsar locations in the obsid field between some time range. Launches search pipeline when done

    Parameters:
    -----------
    obsid: int
        The observation ID
    DI_dir: str
        The directory containing the Jones matrices solutions
    cal_obs:
        The calibration observation ID
    args: object
        The args object generated from argparse. Used for documentation purposes
    psrbeg: int
        The beginning of the observing time
    psrend: int
        The end of the observing time
    product_dir: string
        OPTIONAL - The base directory to store data products. Default = '/group/mwaops/vcs'
    mwa_search_version: string
        OPTIONAL - The version of mwas_search to use. Default = 'master'
    vcstools_version: string
        OPTIONAL - The version of vcstools to use. Default = 'master'
    """
    base_dir = os.path.join(product_dir, obsid, "dpp_pointings")
    nfiles = ( psrend - psrbeg + 1 ) // 200
    if ( ( psrend - psrbeg + 1 )%200 != 0 ):
        nfiles += 1

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
    sp_pointing_list = []
    sp_name_list = []
    for pi, pulsar_line in enumerate(obs_psrs):
        vdif_check = False
        sp_check = False

        PSRJ = pulsar_line[0]
        #See if pulsar is in beam for times
        enter, leave = snfe.pulsar_beam_coverage(obsid, PSRJ, beg=psrbeg, end=psrend)
        if enter is None or leave is None:
            print("{0} is not in beam for time range. Not folding".format(PSRJ))
            continue

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

    print('\nSENDING OFF PULSAR PROCESSING')
    print('----------------------------------------------------------------------------------------')
    # Send off pulsar search
    relaunch_script = 'mwa_search_pipeline.py -o {0} -O {1} --DI_dir {2} -b {3} -e {4} --cand_type Pulsar --vcstools_version {5} --mwa_search_version {6} --channels'.format(obsid, cal_obs, DI_dir, psrbeg, psrend, vcstools_version, mwa_search_version)
    for ch in channels:
        relaunch_script = "{0} {1}".format(relaunch_script, ch)
    search_opts = search_pipe.search_options_class(obsid, cal_id=cal_obs,
                              begin=psrbeg, end=psrend, channels=channels,
                              args=args, DI_dir=DI_dir, relaunch_script=relaunch_script,
                              search_ver=mwa_search_version,
                              vcstools_ver=vcstools_version,
                              data_process=True)
    pulsar_pointing_dirs = [os.path.join(base_dir, s) for s in pulsar_pointing_list]
    for pdir in pulsar_pointing_dirs:
        # Check if fits files are there
        if len(glob.glob("{0}/{1}_*fits".format(pdir, obsid))) < nfiles:
            logger.error("Can not find the {0} expected files in {1}. Exiting".format(nfiles, pdir))
            sys.exit(1)
    for ppd, pnl in zip(pulsar_pointing_dirs, pulsar_name_list):
        for pulsar_name in pnl:
            # Not sure if this works with extended array obs with gridded pointings
            search_pipe.multibeam_binfind(search_opts, [ppd], None, pulsar_name)

    print('\nSENDING OFF VDIF PULSAR PROCESSING')
    print('----------------------------------------------------------------------------------------')
    #Send off vdif pulsar search
    relaunch_script = "{0} --vdif".format(relaunch_script)
    search_opts = search_pipe.search_options_class(obsid, cal_id=cal_obs,
                              begin=psrbeg, end=psrend, channels=channels,
                              args=args, DI_dir=DI_dir, relaunch_script=relaunch_script,
                              search_ver=mwa_search_version,
                              vcstools_ver=vcstools_version,
                              vdif=True, data_process=True)
    vdif_pointing_dirs = [os.path.join(base_dir, s) for s in vdif_pointing_list]
    for pdir in vdif_pointing_dirs:
        # Check if fits files are there
        if len(glob.glob("{0}/{1}_*fits".format(pdir, obsid))) < nfiles:
            print("Can not find the {0} expected files in {1}. Exiting".format(nfiles, pdir))
            sys.exit(1)
    for vpd, vnl in zip(vdif_pointing_dirs, vdif_name_list):
        for vdif_name in vnl:
            # Not sure if this works with extended array obs with gridded pointings
            search_pipe.multibeam_binfind(search_opts, [vpd], None, vdif_name)

    return

if __name__ == "__main__":
    loglevels = dict(DEBUG=logging.DEBUG,
                     INFO=logging.INFO,
                     WARNING=logging.WARNING)
    parser = argparse.ArgumentParser(description="""
    Beamforms, folds on all known pulsars for an observation. If a pulsar is detected it's uploaded to the MWA Pulsar Database.

    Based on a script written by Mengyao Xue.
    """)
    parser.add_argument('-o', '--obsid', type=str,
            help='The observation ID of the fits file to be searched')
    parser.add_argument("--DI_dir", default=None,
            help="Directory containing either Direction Independent Jones Matrices (as created by the RTS) or calibration_solution.bin as created by Andre Offringa's tools.[no default]")
    parser.add_argument('-O','--cal_obs', type=int,
            help="Observation ID of calibrator you want to process.", default=None)
    parser.add_argument("-b", "--begin", type=int, help="First GPS time to process [no default]")
    parser.add_argument("-e", "--end", type=int, help="Last GPS time to process [no default]")
    parser.add_argument("-a", "--all", action="store_true", default=False,
            help="Perform on entire observation span. Use instead of -b & -e")
    parser.add_argument("-n", "--no_comb_check", action="store_true",
            help="Don't do a recombined files check")
    parser.add_argument("--cal_dir_to_tar", type=str, default=None,
            help="The directory (usually rts) of the calibration solutions to upload. If unsupplied, will search for it.")
    parser.add_argument("--srclist", type=str, default=None,
            help="The pathname of the sourcelist to upload as part of the calibration solutions. If unsupplied, will search for it.")
    parser.add_argument("--no_upload", action="store_true",
            help="Use this tag if you don't want to upload calibration solutions to the pulsar database")
    parser.add_argument("-r", "--relaunch", action="store_true", default=False,
            help="Will relaunch searches is they have already completed")
    parser.add_argument('--mwa_search_version', type=str, default='master',
            help="The module version of mwa_search to use. Default: master")
    parser.add_argument('--vcstools_version', type=str, default='master',
            help="The module version of vcstools to use. Default: master")
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
        sys.exit(1)

    if not args.cal_obs:
        print("Please input calibration observation id by setting -O or --cal_obs. Exiting")
        sys.exit(1)

    comp_config = config.load_config_file()
    if not args.DI_dir:
        args.DI_dir = "{0}/{1}/cal/{2}/rts/".format(comp_config['base_product_dir'],
                                                    args.obsid, args.cal_obs)
        print("No DI_dir given so assuming {0} is the directory".format(args.DI_dir))

    if not args.no_upload:
        cal_path, srclist = search_for_cal_srclist(obsid, cal_id, all_cal_returns=False, all_srclist_returns=False)
        std.upload_cal_files(args.obsid, args.cal_obs, cal_path, srclist)

    if args.begin and args.end:
        beg = args.begin
        end = args.end
    elif args.all:
        beg, end = obs_max_min(args.obsid)
    else:
        find_beg_end(args.obsid, base_path=comp_config['base_product_dir'])

    #Perform data checks
    dur = end - beg + 1
    if not args.no_comb_check:
        check = check_data(args.obsid, beg=beg, dur=dur)
        if not check:
            logger.error("Recombined check has failed. Cannot continue.")
            sys.exit(1)
        else:
            logger.info("Recombined check passed, all files present.")

    submit_folds(args.obsid, args.DI_dir, args.cal_obs, args, beg, end,
                      product_dir=comp_config['base_data_dir'],
                      mwa_search_version=args.mwa_search_version,
                      vcstools_version=args.vcstools_version,
                      relaunch=args.relaunch)


