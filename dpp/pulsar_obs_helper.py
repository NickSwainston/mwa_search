#!/usr/bin/env python3
import logging


import find_pulsar_in_obs as fpio
import sn_flux_est as snfe


logger = logging.getLogger(__name__)

def find_fold_times(pulsar, obsid, beg, end, min_z_power=(0.3, 0.1)):
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
    pow_dict, _ = find_pulsars_power(
        obsid, powers=min_z_power, names_ra_dec=names_ra_dec
    )
    for power in pow_dict.keys():
        psr_list = pow_dict[power][obsid]
        enter = None
        leave = None
        if psr_list:  # if pulsar is in beam for this power coverage
            this_enter, this_leave = snfe.pulsar_beam_coverage(
                obsid, pulsar, beg=beg, end=end, min_z_power=power
            )
            if this_enter is not None and this_leave is not None:
                enter = this_enter
                leave = this_leave
                break

    return [enter, leave, power]


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
        # try this if powers isn't iterable
        powers = list(powers)

    if names_ra_dec is None:
        names_ra_dec = np.array(fpio.grab_source_alog(max_dm=250))

    pulsar_power_dict = {}
    for pwr in powers:
        obs_data, meta_data = fpio.find_sources_in_obs(
            [obsid], names_ra_dec, dt_input=100, min_power=pwr
        )
        pulsar_power_dict[pwr] = obs_data

    return pulsar_power_dict, meta_data
