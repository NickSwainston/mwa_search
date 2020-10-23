import logging
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import os
import glob
import math
import psrqpy
import sys

# vcstools imports
import find_pulsar_in_obs as fpio
import sn_flux_est as snfe
from config_vcs import load_config_file
comp_config = load_config_file()
from mwa_metadb_utils import get_common_obs_metadata
from mwa_metadb_utils import obs_max_min, get_obs_array_phase
import checks
from vcstools import data_load

# mwa_search imports
from mwa_search.grid_tools import get_grid
from mwa_search.obs_tools import calc_ta_fwhm

logger = logging.getLogger(__name__)

def _argcheck_find_fold_times(pulsar, obsid, beg, end, min_z_power):
    """Checks that the arguments for find_fold_times are valid"""
    # Pulsar
    if not isinstance(pulsar, str):
        raise TypeError(f"Pulsar is not a string value: {pulsar}")
    # Obsid
    if not isinstance(obsid, int):
        try:
            obsid = int(obsid)
            logger.warn("Obsid had to be converted to int. This may be evidence of a bug")
        except (ValueError, TypeError) as e:
            e.message = f"Invalid Observation ID: {obsid}. Cannot be converted to int"
            raise
    # Beg and end
    if not isinstance(beg, int):
        try:
            beg = int(beg)
            logger.warn("Begin time had to be converted to int. This may be evidence of a bug")
        except (ValueError, TypeError) as e:
            e.message = f"Invalid begin time: {beg}. Cannot be converted to int"
            raise
    if not isinstance(end, int):
        try:
            end = int(end)
            logger.warn("End time had to be converted to int. This may be evidence of a bug")
        except (ValueError, TypeError) as e:
            e.message = f"Invalid end time: {end}. Cannot be converted to int"
            raise
    if beg>end:
        raise ValueError(f"Begining time {beg} greater than end time {end}")
    #convert min_z_power
    try:
        if not isinstance(min_z_power, list):
            min_z_power = list(min_z_power)
    except TypeError as e:
        e.message = f"Invalid min_z_power: {min_z_power}"
        raise
    return pulsar, obsid, beg, end, min_z_power

def find_fold_times(pulsar, obsid, beg, end, min_z_power=(0.3, 0.1), metadata=None, full_meta=None):
    """
    Finds the fractional time the pulsar is in the beam at some zenith normalized power

    Parameters:
    -----------
    pulsar: string
        Pulsar J name
    obsid: int
        The observation ID
    beg: int
        The beginning of the observation time in gps time
    end: int
        The end of the observation time in gps time
    min_z_power: tuple/list
        OPTIONAL - evaluated the pulsar as 'in the beam' at this normalized zenith power. If None will use [0.3, 0.1] Default: None

    Returns:
    [enter, leave, power]: list
        enter: float
            The time the pulsar enters the beam as a normalized fraction of beg and end. None if pulsar not in beam
        leave: float
            The time the pulsar leaves the beam as a normalized fraction of beg and end. None if pulsar not in beam
        power: float
            The power for which enter and leave are calculated
    """
    pulsar, obsid, beg, end, min_z_power = _argcheck_find_fold_times(pulsar, obsid, beg, end, min_z_power)
    min_z_power = sorted(min_z_power, reverse=True)
    names_ra_dec = fpio.grab_source_alog(pulsar_list=[pulsar])
    pow_dict, _ = find_pulsars_power(obsid, powers=min_z_power,
                                     names_ra_dec=names_ra_dec, metadata_list=[[metadata, full_meta]])
    for power in pow_dict.keys():
        psr_list = pow_dict[power][obsid]
        enter = None
        leave = None
        if psr_list:  # if pulsar is in beam for this power coverage
            this_enter, this_leave = snfe.pulsar_beam_coverage(
                obsid, pulsar, beg=beg, end=end, min_z_power=power,
                metadata=metadata, full_meta=full_meta)
            if this_enter is not None and this_leave is not None:
                enter = this_enter
                leave = this_leave
                break

    return [enter, leave, power]


def find_pulsars_power(obsid, powers=None, names_ra_dec=None, metadata_list=None):
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
    metadata: list
        A list of the output of get_common_obs_metadata for the input obsid

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
        # try this if powers isn't iterable
        powers = list(powers)

    if names_ra_dec is None:
        names_ra_dec = np.array(fpio.grab_source_alog(max_dm=250))

    pulsar_power_dict = {}
    for pwr in powers:
        obs_data, meta_data = fpio.find_sources_in_obs(
            [obsid], names_ra_dec,
            dt_input=100, min_power=pwr, metadata_list=metadata_list)
        pulsar_power_dict[pwr] = obs_data

    return pulsar_power_dict, meta_data


def find_beg_end(obsid, base_path=comp_config["base_data_dir"]):
    """
    looks through the comined files of the obsid to find the beginning and end gps times

    Parameters:
    -----------
    obsid: int
        The observation ID
    base_path: string
        OPTIONAL - The system's base working pat

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
        meta_data = get_common_obs_metadata(obsid)
        channels = meta_data[-1]
        combined_files = glob.glob("{0}/{1}/combined/{1}*_ch{2}.dat".\
                                   format(base_path, obsid, channels[-1]))
    comb_times = []
    for comb in combined_files:
        comb_times.append(int(comb.split("_")[1]))
    beg = min(comb_times)
    end = max(comb_times)

    return beg, end


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


def find_pulsars_in_fov(obsid, psrbeg, psrend, fwhm=None, search_radius=0.02):
    """
    Find all pulsars in the field of view and return all the pointings sorted into vdif and normal lists:
    -----------
    obsid: int
        The observation ID
    psrbeg: int
        The begining of the observation you are processing in GPS time
    psrend: int
        The end of the observation you are processing in GPS time
    fwhm: float
        The FWHM of the beam in degrees.
        Default None: Value will be estimated
    search_radius: float
        The radius to search (create beams within) in degrees to account for ionosphere.
        Default: 0.02 degrees

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
    period_query = psrqpy.QueryATNF(params=["PSRJ", "P0"], psrs=pulsar_list,
                               loadfromdb=data_load.ATNF_LOC).pandas

    if fwhm is None:
        # Estimate FWHM
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
        period = period_query[period_query['PSRJ'] == PSRJ].reset_index()["P0"][0]

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
        pointing_list_list = get_pointings_required(raj, decj, fwhm, search_radius)

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
    print(fermi_list)
    pulsar_search_name_list = pulsar_search_name_list + fermi_list[0]
    pulsar_search_pointing_list = pulsar_search_pointing_list + fermi_list[1]

    # Find all of the points of interest candidates
    #-----------------------------------------------------------------------------------------------
    poi_list = get_sources_in_fov(obsid, 'POI', fwhm)
    print(poi_list)
    pulsar_search_name_list = pulsar_search_name_list + poi_list[0]
    pulsar_search_pointing_list = pulsar_search_pointing_list + poi_list[1]

    # Sometimes we get redundant RRATs that are found in RRAT and ANTF catalogues so they need to be removed
    sp_name_list     = list(dict.fromkeys([ ":".join(s) for s in sp_name_list]))
    sp_pointing_list = list(dict.fromkeys(sp_pointing_list))

    # Changing the format of the names list to make it easier to format
    return [[ ":".join(s) for s in pulsar_name_list],
            pulsar_pointing_list,
            [ ":".join(s) for s in vdif_name_list],
            vdif_pointing_list,
            [ ":".join(s) for s in pulsar_search_name_list],
            pulsar_search_pointing_list,
            sp_name_list,
            sp_pointing_list]


def find_pulsars_in_fov_main(kwargs):
    import csv
    #option parsing
    if not kwargs["obsid"]:
        raise ValueError("Please input observation id by setting -o or --obsid. Exiting")

    if not (kwargs["begin"] and kwargs["end"]):
        kwargs["beg"], kwargs["end"] = obs_max_min(args.obsid)

    output_list = find_pulsars_in_fov(kwargs["obsid"], kwargs["begin"], kwargs["end"], fwhm=kwargs["fwhm"], search_radius=kwargs["search_radius"])
    with open(f"{kwargs['obsid']}_fov_sources_temp.csv", 'w', newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',')
        for ol in output_list:
            if len(ol) == 0:
                # Print a space to the line which prevents Nextflow formatting erorrs
                spamwriter.writerow([" "])
            else:
                spamwriter.writerow(ol)

    with open(f'{kwargs["obsid"]}_fov_sources_temp.csv', 'r') as readfile:
        csv_read = readfile.readlines()

    with open(f'{kwargs["obsid"]}_fov_sources.csv', 'w') as csvfile:
        for line in csv_read:
            if len(line) == 0:
                csvfile.write(" \n")
            else:
                csvfile.write(line)