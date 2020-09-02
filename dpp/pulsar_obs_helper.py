#!/usr/bin/env python3
import logging

# vcstools imports
import find_pulsar_in_obs as fpio
import sn_flux_est as snfe
from config_vcs import load_config_file
comp_config = load_config_file()
from mwa_metadb_utils import get_common_obs_metadata


logger = logging.getLogger(__name__)

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
    if min_z_power is None:
        min_z_power = [0.3, 0.1]
    if not isinstance(min_z_power, list):
        min_z_power = list(min_z_power)

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