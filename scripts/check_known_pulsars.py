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

# vcstools and mwa_search imports
import find_pulsar_in_obs as fpio
import mwa_search_pipeline as search_pipe
from mwa_metadb_utils import get_common_obs_metadata as get_meta
from mwa_metadb_utils import obs_max_min, get_obs_array_phase
import config
from grid import get_grid
import checks
#import sn_flux_est as snfe

import logging
logger = logging.getLogger(__name__)

#get ATNF db location
try:
    ATNF_LOC = os.environ['PSRCAT_FILE']
except KeyError:
    logger.warn("ATNF database could not be loaded on disk. This may lead to a connection failure")
    ATNF_LOC = None


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
        logger.error("Recombined files check has failed. Cannot continue")
        sys.exit(1)
    else:
        logger.info("Recombined check passed")

    return


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
        A list of pointings were each pointing contains an RA and a Dec in the format 'hh:mm:ss.ss'
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


def beamform_and_fold(obsid, DI_dir, cal_obs, args, psrbeg, psrend,
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


    #obsbeg, obsend, obsdur = file_maxmin.print_minmax(obsid)

    #wrapping for find_pulsar_in_obs.py
    names_ra_dec = np.array(fpio.grab_source_alog(max_dm=250))
    obs_data, meta_data = fpio.find_sources_in_obs([obsid], names_ra_dec, dt_input=100)
    channels = meta_data[-1][-1]

    #get all the pulsars periods
    pulsar_list = []
    for o in obs_data[obsid]:
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
    for pi, pulsar_line in enumerate(obs_data[obsid]):
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
            for pulsar_l in obs_data[obsid]:
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
                vdif_pointing_list.append("{0} {1}".format(prd[0], prd[1]))
            elif sp_check:
                sp_name_list.append(jname_temp_list)
                sp_pointing_list.append("{0} {1}".format(prd[0], prd[1]))
            else:
                pulsar_name_list.append(jname_temp_list)
                pulsar_pointing_list.append("{0} {1}".format(prd[0], prd[1]))
    
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
    search_pipe.beamform(search_opts, pulsar_pointing_list,
                         pulsar_list_list=pulsar_name_list)

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
    search_pipe.beamform(search_opts, vdif_pointing_list,
                         pulsar_list_list=vdif_name_list)


    #Get the rest of the singple pulse search canidates
    #-----------------------------------------------------------------------------------------------------------
    orig_names_ra_dec = fpio.grab_source_alog(source_type='RRATs',
                                              max_dm=250, include_dm=True)
    # remove any RRATs without at least arc minute accuracy
    names_ra_dec = np.array([s for s in orig_names_ra_dec if (len(s[1]) > 4 and len(s[2]) > 4)])
    obs_data, meta_data = fpio.find_sources_in_obs([obsid], names_ra_dec, dt_input=100)
    
    for pulsar_line in obs_data[obsid]:
        jname = pulsar_line[0]
        for line in names_ra_dec:
            if jname == line[0]:
                jname, raj, decj, _ = line
        jname_temp_list = [jname]

        # grid the pointings to fill 2 arcminute raduis to account for ionosphere shift
        pointing_list_list = get_pointings_required(raj, decj, fwhm, 2./60.)
               
        # sort the pointings into the right groups
        for prd in pointing_list_list:
            sp_name_list.append(jname_temp_list)
            sp_pointing_list.append("{0} {1}".format(prd[0], prd[1]))
    
    print('\nSENDING OFF RRAT SINGLE PULSE SEARCHS')
    print('----------------------------------------------------------------------------------------')
    # Send off pulsar search
    relaunch_script = 'mwa_search_pipeline.py -o {0} -O {1} --cand_type RRATs --DI_dir {2} -b {3} -e {4} --single_pulse --vcstools_version {5} --mwa_search_version {6} --channels'.format(obsid, cal_obs, DI_dir, psrbeg, psrend, vcstools_version, mwa_search_version)
    for ch in channels:
        relaunch_script = "{0} {1}".format(relaunch_script, ch)
    search_opts = search_pipe.search_options_class(obsid, cal_id=cal_obs,
                              begin=psrbeg, end=psrend, channels=channels,
                              args=args, DI_dir=DI_dir, relaunch_script=relaunch_script,
                              search_ver=mwa_search_version,
                              vcstools_ver=vcstools_version,
                              single_pulse=True, cand_type='RRATs')
    search_pipe.beamform(search_opts, sp_pointing_list,
                         pulsar_list_list=sp_name_list,
                         code_comment="RRATs single pulse search",
                         relaunch=relaunch)


    # Find all of the FRB candidates
    #-----------------------------------------------------------------------------------------------------------
    orig_names_ra_dec = fpio.grab_source_alog(source_type='FRB',
                                              max_dm=10000, include_dm=True)
    # remove any RRATs without at least arc minute accuracy
    names_ra_dec = np.array([s for s in orig_names_ra_dec if (len(s[1]) > 4 and len(s[2]) > 4)])
    obs_data, meta_data = fpio.find_sources_in_obs([obsid], names_ra_dec, dt_input=100)
    
    sp_name_list = []
    sp_pointing_list = []
    for pulsar_line in obs_data[obsid]:
        jname = pulsar_line[0]
        for line in names_ra_dec:
            if jname == line[0]:
                jname, raj, decj, dm = line
        jname_temp_list = [jname]

        # grid the pointings to fill 2 arcminute raduis to account for ionosphere shift
        pointing_list_list = get_pointings_required(raj, decj, fwhm, 2./60.)
               
        # sort the pointings into the right groups
        for prd in pointing_list_list:
            sp_name_list.append(jname_temp_list)
            sp_pointing_list.append("{0} {1}".format(prd[0], prd[1]))
    
    print('\nSENDING OFF FRB SINGLE PULSE SEARCHS')
    print('----------------------------------------------------------------------------------------')
    # Send off pulsar search
    relaunch_script = 'mwa_search_pipeline.py -o {0} -O {1} --cand_type FRB --DI_dir {2} -b {3} -e {4} --single_pulse --vcstools_version {5} --mwa_search_version {6} --channels'.format(obsid, cal_obs, DI_dir, psrbeg, psrend, vcstools_version, mwa_search_version)
    for ch in channels:
        relaunch_script = "{0} {1}".format(relaunch_script, ch)
    search_opts = search_pipe.search_options_class(obsid, cal_id=cal_obs,
                              begin=psrbeg, end=psrend, channels=channels,
                              args=args, DI_dir=DI_dir, relaunch_script=relaunch_script,
                              search_ver=mwa_search_version,
                              vcstools_ver=vcstools_version,
                              single_pulse=True, cand_type='FRB')
    search_pipe.beamform(search_opts, sp_pointing_list,
                         pulsar_list_list=sp_name_list,
                         code_comment="FRB single pulse search",
                         relaunch=relaunch)


    # Find all of the Fermi candidates
    #-----------------------------------------------------------------------------------------------------------
    orig_names_ra_dec = fpio.grab_source_alog(source_type='Fermi',
                                              max_dm=10000, include_dm=True)
    # remove any RRATs without at least arc minute accuracy
    names_ra_dec = np.array([s for s in orig_names_ra_dec if (len(s[1]) > 4 and len(s[2]) > 4)])
    obs_data, meta_data = fpio.find_sources_in_obs([obsid], names_ra_dec, dt_input=100)
    
    sp_name_list = []
    sp_pointing_list = []
    for pulsar_line in obs_data[obsid]:
        jname = pulsar_line[0]
        for line in names_ra_dec:
            if jname == line[0]:
                jname, raj, decj, pos_u = line
        jname_temp_list = [jname]

        # grid the pointings to fill the position uncertaint (given in arcminutes)
        pointing_list_list = get_pointings_required(raj, decj, fwhm, float(pos_u)/60.)
               
        # sort the pointings into the right groups
        for prd in pointing_list_list:
            sp_name_list.append(jname_temp_list)
            sp_pointing_list.append("{0} {1}".format(prd[0], prd[1]))
    
    print('\nSENDING OFF FERMI CANDIDATE SEARCHS')
    print('----------------------------------------------------------------------------------------')
    # Send off pulsar search
    relaunch_script = 'mwa_search_pipeline.py -o {0} -O {1} --cand_type Fermi --DI_dir {2} -b {3} -e {4} --search --vcstools_version {5} --mwa_search_version {6} --channels'.format(obsid, cal_obs, DI_dir, psrbeg, psrend, vcstools_version, mwa_search_version)
    for ch in channels:
        relaunch_script = "{0} {1}".format(relaunch_script, ch)
    search_opts = search_pipe.search_options_class(obsid, cal_id=cal_obs,
                              begin=psrbeg, end=psrend, channels=channels,
                              args=args, DI_dir=DI_dir, relaunch_script=relaunch_script,
                              search_ver=mwa_search_version,
                              vcstools_ver=vcstools_version,
                              search=True, cand_type='Fermi')
    search_pipe.beamform(search_opts, sp_pointing_list,
                         pulsar_list_list=sp_name_list,
                         code_comment="Fermi candidate pulsar search",
                         relaunch=relaunch)

    # Find all of the points of interest candidates
    #-----------------------------------------------------------------------------------------------
    names_ra_dec = fpio.grab_source_alog(source_type='POI')
    obs_data, meta_data = fpio.find_sources_in_obs([obsid], names_ra_dec, dt_input=100)
    
    poi_name_list = []
    poi_pointing_list = []
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
            poi_name_list.append(jname_temp_list)
            poi_pointing_list.append("{0} {1}".format(prd[0], prd[1]))
    
    print('\nSENDING OFF POINTS OF INTEREST SEARCHS')
    print('----------------------------------------------------------------------------------------')
    # Send off pulsar search
    relaunch_script = "mwa_search_pipeline.py -o {0} -O {1} --cand_type POI " \
                      "--DI_dir {2} -b {3} -e {4} --search --vcstools_version {5} " \
                      "--mwa_search_version {6} --channels".format(obsid, cal_obs,
                        DI_dir, psrbeg, psrend, vcstools_version, mwa_search_version)
    for ch in channels:
        relaunch_script = "{0} {1}".format(relaunch_script, ch)
    search_opts = search_pipe.search_options_class(obsid, cal_id=cal_obs,
                              begin=psrbeg, end=psrend, channels=channels,
                              args=args, DI_dir=DI_dir, relaunch_script=relaunch_script,
                              search_ver=mwa_search_version,
                              vcstools_ver=vcstools_version,
                              search=True, cand_type='POI')
    search_pipe.beamform(search_opts, poi_pointing_list,
                         pulsar_list_list=poi_name_list,
                         code_comment="Points of interest candidate pulsar search",
                         relaunch=relaunch)

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
        quit()

    if not args.cal_obs:
        print("Please input calibration observation id by setting -O or --cal_obs. Exiting")
        quit()

    comp_config = config.load_config_file()
    if not args.DI_dir:
        args.DI_dir = "{0}/{1}/cal/{2}/rts/".format(comp_config['base_product_dir'],
                                                    args.obsid, args.cal_obs)
        print("No DI_dir given so assuming {0} is the directory".format(args.DI_dir))

    if args.begin and args.end:
        beg = args.begin
        end = args.end
    elif args.all:
        beg, end = obs_max_min(args.obsid)
    else:
        find_beg_end(args.obsid, base_path=comp_config['base_product_dir'])

    #Perform data checks
    dur = end-beg
    if not args.no_comb_check:
        check_data(args.obsid, beg=beg, dur=dur)

    beamform_and_fold(args.obsid, args.DI_dir, args.cal_obs, args, beg, end,
                      product_dir=comp_config['base_product_dir'],
                      mwa_search_version=args.mwa_search_version,
                      vcstools_version=args.vcstools_version,
                      relaunch=args.relaunch)



