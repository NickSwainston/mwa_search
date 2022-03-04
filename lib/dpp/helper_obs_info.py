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
from vcstools.radiometer_equation import multi_psr_snfe
from vcstools.metadb_utils import get_common_obs_metadata, obs_max_min, get_obs_array_phase
from vcstools import data_load
from vcstools.pointing_utils import format_ra_dec
from vcstools.catalogue_utils import grab_source_alog, deg2sex
from vcstools.beam_calc import find_sources_in_obs, source_beam_coverage_and_times
from vcstools.config import load_config_file

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


def reformat_psrs_pointings(psr_pointing_list):
    """Reformats the output from find_pulsars_in_fov for the dpp"""
    _dict = {}
    # Loop over the concatenated lists of regular and vdif pulsars/pointings
    # Multiple pulsars per pointing are listed as psr1:psr2 etc... so we need to fix for that
    for psr_per_point, pointing in zip(psr_pointing_list[0] + psr_pointing_list[2], psr_pointing_list[1] + psr_pointing_list[3]):
        psr_list = psr_per_point.split(":")
        for psr in psr_list:
            if psr not in _dict.keys():
                _dict[psr] = []
            _dict[psr].append(pointing)
    return _dict


def find_fold_times(pulsars, obsid, beg, end, min_z_power=(0.3, 0.1), metadata=None, full_meta=None, query=None):
    """
    Finds the fractional time the pulsar is in the beam at some zenith normalized power

    Parameters:
    -----------
    pulsar: list
        Pulsar J names of pulsars to evaluate
    obsid: int
        The observation ID
    beg: int
        The beginning of the observation time in gps time
    end: int
        The end of the observation time in gps time
    min_z_power: tuple/list
        OPTIONAL - evaluated the pulsar as 'in the beam' at this normalized zenith power. If None will use [0.3, 0.1] Default: None

    Returns:
    fold_times_dict: dict
        keys:
            psr: dict
                keys:
                    power: float
                        The power of the recorded enter and leave times
                    enter: float
                        The normalised time the pulsar enters the beam
                    exit: float
                        The normalisd time the pulsar leaves the beam
    """
    if not metadata or not full_meta:
        metadata, full_meta = get_common_obs_metadata(obsid, return_all=True)
    obs_beg, obs_end = obs_max_min(obsid)
    min_z_power = sorted(min_z_power, reverse=True)
    names_ra_dec = grab_source_alog(pulsar_list=pulsars)
    pow_dict, _ = find_pulsars_power(obsid, powers=min_z_power, names_ra_dec=names_ra_dec, metadata_list=[[metadata, full_meta]])
    fold_time_dict = {}
    for psr in pulsars:
        fold_time_dict[psr] = {}
        for power in pow_dict.keys():
            fold_time_dict[psr]["power"] = power
            fold_time_dict[psr]["enter"] = None
            fold_time_dict[psr]["leave"] = None
            psr_list = pow_dict[power][obsid]
            if psr_list:  # if pulsar is in beam for this power coverage
                files_beg_norm, files_end_norm = source_beam_coverage_and_times(
                    obsid, psr, files_beg=beg, files_end=end,
                    obs_beg=obs_beg, obs_end=obs_end,
                    min_z_power=power, query=query,
                    common_metadata=metadata)[4:-3]
                if files_beg_norm is not None and files_end_norm is not None:
                    fold_time_dict[psr]["enter"] = files_beg_norm
                    fold_time_dict[psr]["leave"] = files_end_norm
                    break

    return fold_time_dict

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
        OPTIONAL - A list of puslars and their RA and Dec values to evaluate (generated from get_source_alog).
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
        names_ra_dec = np.array(grab_source_alog(max_dm=250))

    pulsar_power_dict = {}
    if len(names_ra_dec) > 0:
        for pwr in powers:
            obs_data, meta_data = find_sources_in_obs(
                [obsid], names_ra_dec,
                dt_input=100, min_z_power=pwr, metadata_list=metadata_list)
            pulsar_power_dict[pwr] = obs_data

    return pulsar_power_dict, meta_data


def find_beg_end(obsid, base_path=None):
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
    if base_path is None:
        comp_config = load_config_file()
        base_path = comp_config["base_data_dir"]
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
    rads, decds = get_grid(rar, decr, np.radians(grid_sep), loops, verbose=False)

    #convert back to sexidecimals
    coord = SkyCoord(rads,decds,unit=(u.deg,u.deg))
    rajs = coord.ra.to_string(unit=u.hour, sep=':')
    decjs = coord.dec.to_string(unit=u.degree, sep=':')
    temp = []
    for raj, decj in zip(rajs, decjs):
        temp.append([raj, decj])
    pointing_list_list = format_ra_dec(temp, ra_col = 0, dec_col = 1)
    return pointing_list_list


def get_sources_in_fov(obsid, source_type, fwhm):
    """
    Find all sources of the input type in the observations field-of-view

    Parameters:
    -----------
    obsid: str
        observation ID to search in
    source_type: str
        the source type input to grab_source_alog()
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
    names_ra_dec = grab_source_alog(source_type=source_type)
    if len(names_ra_dec) == 0 :
        return [[],[]]
    obs_data, _ = find_sources_in_obs([obsid], names_ra_dec, dt_input=100)
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


def apply_offset(pointing_list, offset, angle_offset):
    """
    Apply an offset to the input list of pointings

    Parameters:
    -----------
    pointing_list: list
        A list of pointings where each pointing contains an RA and a Dec in the format 'hh:mm:ss.ss_dd:mm:ss.ss'
    offset: float
        The offset to apply to all pointings in arcseconds
    angle_offset: float
        The angle of the offset to apply to all pointings in degrees where zero is north

    Returns:
    --------
    offset_pointing_list: list
        A list of pointings where each pointing contains an RA and a Dec in the format 'hh:mm:ss.ss_dd:mm:ss.ss'
    """
    # Create astropy SkyCoords
    rajs  = []
    decjs = []
    for p in pointing_list:
        rajs.append(p.split("_")[0])
        decjs.append(p.split("_")[1])
    orig_pos = SkyCoord(rajs, decjs, frame='icrs', unit=(u.hourangle,u.deg))

    # Apply offset
    offset_pos = orig_pos.directional_offset_by(angle_offset*u.deg, offset*u.arcsec)

    # Convert back to pointing list
    rajs = offset_pos.ra.to_string(unit=u.hour, sep=':')
    decjs = offset_pos.dec.to_string(unit=u.degree, sep=':')
    temp = []
    for raj, decj in zip(rajs, decjs):
        temp.append([raj, decj])
    temp_formated = format_ra_dec(temp, ra_col = 0, dec_col = 1)
    offset_pointing_list = []
    for raj, decj in temp_formated:
        offset_pointing_list.append("{}_{}".format(raj, decj))
    return offset_pointing_list


def find_pulsars_in_fov(obsid, psrbeg, psrend,
                        fwhm=None, search_radius=0.02,
                        meta_data=None, full_meta=None,
                        no_known_pulsars=False, no_search_cands=False,
                        offset=0, angle_offset=0):
    """
    Find all pulsars in the field of view and return all the pointings sorted into vdif and normal lists:

    Parameters:
    -----------
    obsid: int
        The observation ID
    psrbeg: int
        The begining of the observation you are processing in GPS time
    psrend: int
        The end of the observation you are processing in GPS time

    fwhm: float
        OPTIONAL - The FWHM of the beam in degrees.
        Default None: Value will be estimated
    search_radius: float
        OPTIONAL - The radius to search (create beams within) in degrees to account for ionosphere.
        Default: 0.02 degrees
    meta_data: list
        OPTIONAL - Comon observation metadata in the format from vcstools.metadb_utils.get_common_obs_metadata
        Default None: Will perform the metadata call
    full_meta: dict
        OPTIONAL - Full observation metadata in the format from vcstools.metadb_utils.getmeta
        Default None: Will perform the metadata call
    no_known_pulsars: bool
        OPTIONAL - Will return no known pulsars
        Default: False
    no_search_cands: bool
        OPTIONAL - Will return no search candidates
        Default: False

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
    if not meta_data or not full_meta:
        meta_data, full_meta = get_common_obs_metadata(obsid, return_all=True)
    channels = meta_data[-1]

    if fwhm is None:
        # Estimate FWHM
        oap = get_obs_array_phase(obsid)
        centrefreq = 1.28 * float(min(channels) + max(channels)) / 2.
        fwhm = calc_ta_fwhm(centrefreq, array_phase=oap)

    # Find all pulsars in beam at at least 0.3 and 0.1 of zenith normlaized power
    names_ra_dec = np.array(grab_source_alog(max_dm=250))
    pow_dict, _ = find_pulsars_power(obsid, powers=[0.3, 0.1], names_ra_dec=names_ra_dec)
    obs_psrs = pow_dict[0.3][obsid]
    # Find pulsars with power between 0.3 and 0.1 and calculate their SN
    psrs_list_03 = [x[0] for x in obs_psrs]
    psrs_list_01 = [x[0] for x in pow_dict[0.1][obsid]]
    psrs_03_01 = psrs_list_01
    for psr in psrs_list_03:
        if psr in psrs_list_01:
            psrs_03_01.remove(psr)
    sn_dict_01 = multi_psr_snfe(psrs_03_01, obsid, psrbeg, psrend, min_z_power=0.1, common_metadata=meta_data)
    # Include all bright pulsars in beam at at least 0.1 of zenith normalized power
    for psr in psrs_03_01:
        sn, sn_err, _, _ = sn_dict_01[psr]
        if sn is not None and sn_err is not None:
            if sn - sn_err >= 10.:
                for psr_list in pow_dict[0.1][obsid]:
                    if psr in psr_list:
                        obs_psrs.append(psr_list)

    #get all the pulsars periods
    pulsar_list = []
    for o in obs_psrs:
        pulsar_list.append(o[0])
    period_query = psrqpy.QueryATNF(params=["PSRJ", "P0"], psrs=pulsar_list,
                                    loadfromdb=data_load.ATNF_LOC).pandas

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
        if not (len(PSRJ) < 11 or PSRJ[-1] == 'A' or PSRJ[-2:] == 'aa'):
            continue

        for line in names_ra_dec:
            if PSRJ == line[0]:
                temp = [line]

        temp = format_ra_dec(temp, ra_col = 1, dec_col = 2)
        jname, raj, decj = temp[0]

        #get pulsar period
        period = period_query[period_query['PSRJ'] == PSRJ].reset_index()["P0"][0]

        if math.isnan(period):
            logger.warn("Period not found in ephermeris for {0} so assuming "
                  "it's an RRAT".format(jname))
            sp_check = True
            period = 0.
        elif float(period) < .05:
            vdif_check = True
        period = float(period)*1000.
        logger.debug("{0:12} RA: {1} Dec: {2} Period: {3:8.2f} (ms) Begin {4} End {5}".format(PSRJ,
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
    logger.info(f"{obsid} Fermi candidates: {fermi_list}")
    pulsar_search_name_list = pulsar_search_name_list + fermi_list[0]
    pulsar_search_pointing_list = pulsar_search_pointing_list + fermi_list[1]

    # Find all of the points of interest candidates
    #-----------------------------------------------------------------------------------------------
    poi_list = get_sources_in_fov(obsid, 'POI', fwhm)
    logger.info(f"{obsid} Points of interest: {poi_list}")
    pulsar_search_name_list = pulsar_search_name_list + poi_list[0]
    pulsar_search_pointing_list = pulsar_search_pointing_list + poi_list[1]

    # Sometimes we get redundant RRATs that are found in RRAT and ANTF catalogues so they need to be removed
    sp_name_list     = list(dict.fromkeys([ ":".join(s) for s in sp_name_list]))
    sp_pointing_list = list(dict.fromkeys(sp_pointing_list))

    # Changing the format of the names list to make it easier to format
    pulsar_name_list        = [ ":".join(s) for s in pulsar_name_list]
    vdif_name_list          = [ ":".join(s) for s in vdif_name_list]
    pulsar_search_name_list = [ ":".join(s) for s in pulsar_search_name_list]

    if no_known_pulsars:
        # Return empty list for all known pulsar categories
        pulsar_name_list = []
        pulsar_pointing_list = []
        vdif_name_list = []
        vdif_pointing_list = []

    if no_search_cands:
        # Return empty list for all search candidate categories
        pulsar_search_name_list = []
        pulsar_search_pointing_list = []
        sp_name_list = []
        sp_pointing_list = []

    # Apply offsets
    pulsar_pointing_list        = apply_offset(pulsar_pointing_list,        offset, angle_offset)
    vdif_pointing_list          = apply_offset(vdif_pointing_list,          offset, angle_offset)
    pulsar_search_pointing_list = apply_offset(pulsar_search_pointing_list, offset, angle_offset)
    sp_pointing_list            = apply_offset(sp_pointing_list,            offset, angle_offset)

    return [pulsar_name_list,
            pulsar_pointing_list,
            vdif_name_list,
            vdif_pointing_list,
            pulsar_search_name_list,
            pulsar_search_pointing_list,
            sp_name_list,
            sp_pointing_list]


def find_pulsars_in_fov_main(kwargs):
    import csv
    #option parsing
    if not kwargs["obsid"]:
        raise ValueError("Please input observation id by setting -o or --obsid. Exiting")

    if not (kwargs["begin"] and kwargs["end"]):
        kwargs["beg"], kwargs["end"] = obs_max_min(kwargs["obsid"])

    output_list = find_pulsars_in_fov(kwargs["obsid"], kwargs["begin"], kwargs["end"],
                                      fwhm=kwargs["fwhm"],
                                      search_radius=kwargs["search_radius"],
                                      no_known_pulsars=kwargs["no_known_pulsars"],
                                      no_search_cands=kwargs["no_search_cands"],
                                      offset=kwargs["offset"],
                                      angle_offset=kwargs["angle_offset"])
    if kwargs['n_pointings'] is None:
        with open(f"{kwargs['obsid']}_fov_sources.csv", 'w', newline='') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=',')
            for ol in output_list:
                if len(ol) == 0:
                    # Print a space to the line which prevents Nextflow formatting erorrs
                    #spamwriter.writerow([" "])
                    csvfile.write(" \n")
                else:
                    spamwriter.writerow(ol)
    else:
        name_pointing_pairs = [output_list[i:i + 2] for i in range(0, len(output_list), 2)]
        pair_names = ["pulsar", "vdif", "pulsar_search", "single_pulse"]
        # Loop over each output type
        for i, name in enumerate(pair_names):
            name_list     = output_list[2*i]
            pointing_list = output_list[2*i+1]
            if len(name_list) != 0:
                # split each file into the required length
                pointing_list_chunks = [pointing_list[x:x+kwargs['n_pointings']] for x in range(0, len(pointing_list), kwargs['n_pointings'])]
                for ci in range(len(pointing_list_chunks)):
                    first_id = ci * kwargs['n_pointings'] + 1
                    last_id  = ci * kwargs['n_pointings'] + len(pointing_list_chunks[ci])
                    out_file_name = f"{kwargs['obsid']}_fov_{name}_sources_{first_id}_{last_id}.txt"
                    print(f"Recording {name} sources in {out_file_name}")
                    with open(out_file_name, 'w') as out_file:
                        for out_pointing in pointing_list_chunks[ci]:
                            out_file.write(f"{out_pointing}\n")
