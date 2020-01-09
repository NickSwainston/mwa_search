#!/usr/bin/env python3

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import logging
import argparse
import json
import os
import glob

import prof_utils
import binfinder

logger = logging.getLogger(__name__)

try:
    ATNF_LOC = os.environ['PSRCAT_FILE']
except KeyError:
    logger.warning("ATNF database could not be found on disk.")
    ATNF_LOC = None
EPNDB_LOC = os.environ["EPNDB_LOC"]

#--------------------------------------------------------------------------
def pulsar_db_search(pulsar=None, obsid=None):
    #from mwa_pulsar_client import client
    #address = "https://mwa-pawsey-volt01.pawsey.org.au"
    #auth = ('mwapulsar','veovys9OUTY=')
    #pul_list_dict = client.pulsar_list(address, auth)
    #detection_list = client.detection_list(address, auth)
    #return_stuff = client.detection_get(address, auth, observationid=obsid)
    #print(client.pulsar_get(address, auth, name="J2241-5236"))

    """
    Pulsar list is a list of dictionaries. Each dictionary has the following keys:
    {'name', 'ra', 'dec'}

    Detection list is a list of dictionaries. Each dictionary is a detection with the following keys:
    {'observationid', 'subband', 'observation_type', 'observation_type_str', 'startcchan', 'stopcchan', 'coherent', 'flux', 'flux_error', 'width', 'width_error', 'scattering', 'scattering_error', 'dm', 'dm_error', 'pulsar', 'calibrator', 'version'}

    Detection get is a list of dictionaries. The list contains a dictionary for each detection in the OBSID. Each of these dictionaries has the folloing keys:
    'observationid', 'subband', 'observation_type', 'observation_type_str', 'startcchan', 'stopcchan', 'coherent', 'flux', 'flux_error', 'width', 'width_error', 'scattering', 'scattering_error', 'dm', 'dm_error', 'pulsar', 'calibrator', 'detection_files'

    Pulsar get returns a dict including pulsar's name, ra, dec and a list of detections. This includes all parameters listed above ^
     """
    pass

#--------------------------------------------------------------------------
def roll_data(data, idx_to_roll=None, roll_to=None):
    """
    Rolls a list of data to about some index to some new index.

    Parameters:
    -----------
    data: list
        The data to roll
    idx_to_roll: int
        OPTIONAL - The index number of the original array to roll. If None, uses max value index. Default: None
    roll_to: int
        OPTIONAL - The index to move the 'idx_to_roll' to. If None, rolls to centre of list. Default: None

    Return:
    -------
    [idx_to_roll, roll_to, rolled_data]: list
        idx_to_roll: int
            The index that the data was rolled about
        roll_to: int
            The index of the array that the roll index was moved to
        rolled_data: list
            The rolled data
    """
    if idx_to_roll is None:
        val = max(data)
        idx_to_roll = list(data).index(val)
    if roll_to is None:
        roll_to = int(len(data)/2-1)

    #This is the amount to cycle the data
    roll_n = int(roll_to - idx_to_roll)
    rolled_data = np.roll(data, roll_n)

    return [idx_to_roll, roll_to, rolled_data]

#--------------------------------------------------------------------------
def get_data_from_epndb(pulsar):
    """
    Searches the EPN database and returns all information of the chosen pulsar

    Parameters:
    -----------
    pulsar: string
        The name of the pulsar to search for

    Returns:
    --------
    pulsar_dict: dictionary
        A dictionary in which each value is a list corresponding to a different entry on the databse. Keys:
        I: list
            An list of stokes I values for the profile
        Q: list
            An list of stokes Q values for the profile
        U: list
            An list of stokes U values for the profile
        V: list
            An list of stokes V values for the profile
        freq: float
            The frequency of the observation in MHz
        dm: float
            The measured dispersion measure
        rm: float
            The measured rotation measure. Returns 0.0 for no rm measurement
        site: string
            The location the pulsar was observed at

    """
    #find all directories in epndb
    all_json_dirs = glob.glob(EPNDB_LOC + "/json/*/*/*.json")
    json_pulsar_dirs = []

    for directory in all_json_dirs:
        name = directory.split("/")[-2]
        if name==pulsar:
            json_pulsar_dirs.append(directory)

    logger.debug("All json directories:            {}".format(len(all_json_dirs)))
    logger.debug("Positive json directories:       {}".format(len(json_pulsar_dirs)))

    pulsar_dict={"Ix":[], "Qx":[], "Ux":[], "Vx":[], "Iy":[], "Qy":[], "Uy":[], "Vy":[],\
                "freq":[], "dm":[], "rm":[], "site":[]}
    for i, file_path in enumerate(json_pulsar_dirs):
        with open(file_path) as json_file:
            data = json.load(json_file)
            header = data["hdr"]
            series = data["series"]

            #Ignore anything without frequency info
            if "freq" in header:
                pulsar_dict["freq"].append(float(header["freq"]))
                pulsar_dict["dm"].append(float(header["dm"]))
                pulsar_dict["rm"].append(float(header["rm"]))
                pulsar_dict["site"].append(header["site"])
                pulsar_dict["Ix"].append([k[0] for k in series["I"]])
                pulsar_dict["Iy"].append([k[1] for k in series["I"]])

                if "Q" in series:
                    pulsar_dict["Qx"].append([k[0] for k in series["Q"]])
                    pulsar_dict["Qy"].append([k[1] for k in series["Q"]])
                else:
                    pulsar_dict["Qx"].append(None)
                    pulsar_dict["Qy"].append(None)
                if "U" in series:
                    pulsar_dict["Ux"].append([k[0] for k in series["U"]])
                    pulsar_dict["Uy"].append([k[1] for k in series["U"]])
                else:
                    pulsar_dict["Ux"].append(None)
                    pulsar_dict["Uy"].append(None)
                if "V" in series:
                    pulsar_dict["Vx"].append([k[0] for k in series["V"]])
                    pulsar_dict["Vy"].append([k[1] for k in series["V"]])
                else:
                    pulsar_dict["Vx"].append(None)
                    pulsar_dict["Vy"].append(None)

    #sort by frequency
    freq_list = pulsar_dict["freq"]
    for key in pulsar_dict.keys():
        _, pulsar_dict[key] = (list(t) for t in zip(*sorted(zip(freq_list, pulsar_dict[key]))))

    return pulsar_dict

#--------------------------------------------------------------------------
def plot_bestprof(bestprof, freq=None, out_dir="./"):
    """
    Plots a .pfd.bestprof file from a presto output and saves as a .png

    Parameters:
    -----------
    bestprof: string
        The path to the bestprof file
    freq: float
        OPTIONAL - The frequency of the observation. Default: None
    out_dir: string
        OPTIONAL - The directory to save the file to. Default: './'

    Returns:
    --------
    fig_path: string
        The path of the .png plot
    """
    #retrieve data from bestprof
    logger.info("Plotting profile from file: {0}".format(bestprof))
    y = prof_utils.get_from_bestprof(bestprof)[-2]

    #normalize and align
    y = np.array(y)/max(y)
    y = roll_data(y)[-1]
    x = np.linspace(-0.5, 0.5, len(y))
    info_dict = binfinder.bestprof_info(filename=bestprof)

    #make the title
    title = "{0} {1} Pulse profile".format(info_dict["pulsar"], info_dict["obsid"])
    save_name = "pulse_profile_{0}_{1}".format(info_dict["obsid"], info_dict["pulsar"])
    if freq is not None:
        title += " - {}MHz".format(freq)
        save_name += "_{}MHz".format(freq)
    save_name += ".png"

    #plot
    _, ax = plt.subplots(figsize=(12, 8))
    plt.plot(x, y, color="black")
    plt.title(title)
    plt.text(0.05, 0.95,  "S/N:             {0}".format(info_dict["sn"]), fontsize=10, color="black", transform=ax.transAxes)
    plt.text(0.05, 0.925, "Chi Sq:          {0}".format(info_dict["chi"]), fontsize=10, color="black", transform=ax.transAxes)
    plt.text(0.05, 0.9,   "DM:              {0}".format(info_dict["dm"]), fontsize=10, color="black", transform=ax.transAxes)
    plt.text(0.05, 0.875, "Period (ms):     {0} +/- {1}".format(info_dict["period"], info_dict["period_error"]), fontsize=10,\
            color="black", transform=ax.transAxes)

    fig_path = os.path.join(out_dir, save_name)
    logger.info("Saving bestprof figure: {0}".format(fig_path))
    plt.savefig("{0}".format(fig_path))

    return fig_path

#--------------------------------------------------------------------------
def plot_ascii(archive, pulsar=None, freq=None, obsid=None, out_dir="./"):
    """
    Plots an ascii text file and saves as a .png

    Parameters:
    -----------
    archive: string
        The path to the ascii text file
    pulsar: string
        OPTIONAL - The name of the pulsar
    freq: float
        OPTIONAL - The frequency of the observation. Default: None
    obsid: int
        OPTIONAL - The ID of the MWA observation. Default: None
    out_dir: string
        OPTIONAL - The directory to save the file to. Default: './'

    Returns:
    --------
    fig_path: string
        The path of the .png plot
    """
    #Read the archive
    sI = prof_utils.get_from_ascii(archive)[0]
    logger.info("Plotting profile from file: {0}".format(archive))

    #normalize and align
    sI = np.array(sI)/max(sI)
    sI = roll_data(sI)[-1]
    x = np.linspace(-0.5, 0.5, len(sI))

    #make the title
    title = ""
    save_name = "pulse_profile"
    if pulsar is not None:
        title += "{}".format(pulsar)
        save_name += "_{}".format(pulsar)
    title += " Pulse Profile"
    if obsid is not None:
        title += " {}".format(obsid)
        save_name += "_{}".format(obsid)
    if freq is not None:
        title += " - {}MHz".format(freq)
        save_name += "_{}MHz".format(freq)
    save_name += ".png"

    #plot -
    plt.figure(figsize=(20, 12))
    plt.title(title, fontsize=36)
    plt.xlabel("Pulse Phase", fontsize=20)
    plt.ylabel("Intensity", fontsize=20)
    plt.plot(x, sI, color="black")

    fig_path = os.path.join(out_dir, save_name)
    logger.info("Saving ascii figure: {0}".format(fig_path))
    plt.savefig(fig_path)

    return fig_path

#--------------------------------------------------------------------------
def calc_lin_pa(stokes_Q, stokes_U):
    """
    Calculates the linear polsarization and position angle from stokes Q and U components

    Parameters:
    -----------
    stokes_Q: list
        A list of the Stokes Q values
    stokes_U: list
        A list of the Stokes U values

    Returns:
    --------
    lin_pol: list
        The linear polarization
    pa: list
        The position angles
    """
    lin_pol = []
    pa = []
    for q, u in zip(stokes_Q, stokes_U):
        lin_pol.append(np.sqrt(q**2+u**2))
        pa.append(0.5*np.arctan(q/u))

    return lin_pol, pa

#--------------------------------------------------------------------------
def plot_archive_stokes(archive, pulsar=None, freq=None, obsid=None, out_dir="./"):
    """
    Plots a polarimetry profile as a .png using a full-stokes ascii text file

    Parameters:
    -----------
    archive: string
        The path to the ascii text file
    pulsar: string
        OPTIONAL - The name of the pulsar. Default: None
    freq: float
        OPTIONAL - The frequency of the observation. Default: None
    obsid: int
        OPTIONAL - The observation ID. Default: None
    out_dir: string
        OPTIONAL - The directory to ouput the .png to

    Returns:
    --------
    fig_path: string
        The path of the output .png file
    """

    #make the title
    title = ""
    save_name = "Polarimetry_profile"
    if pulsar is not None:
        title += "{}".format(pulsar)
        save_name += "_{}".format(pulsar)
    title += " Polarimetry Profile"
    if obsid is not None:
        title += " {}".format(obsid)
        save_name += "_{}".format(obsid)
    if freq is not None:
        title += " - {}MHz".format(freq)
        save_name += "_{}MHz".format(freq)
    save_name += ".png"

    #Read the archive
    sI = []
    sQ = []
    sU = []
    sV = []
    f = open(archive)
    lines = iter(f.readlines())
    next(lines) #skip first line
    for line in lines:
        thisline=line.split()
        sI.append(float(thisline[3]))
        sQ.append(float(thisline[4]))
        sU.append(float(thisline[5]))
        sV.append(float(thisline[6]))

    #normalize, aign and find linear polarization and position angle
    max_I = max(sI)
    sI = np.array(sI)/max_I
    sQ = np.array(sQ)/max_I
    sU = np.array(sU)/max_I
    sV = np.array(sV)/max_I
    roll_idx, roll_to, sI = roll_data(sI)
    sQ = roll_data(sQ, idx_to_roll=roll_idx, roll_to=roll_to)[-1]
    sU = roll_data(sU, idx_to_roll=roll_idx, roll_to=roll_to)[-1]
    sV = roll_data(sV, idx_to_roll=roll_idx, roll_to=roll_to)[-1]
    lin_pol, pa = calc_lin_pa(sQ, sU)
    x = np.linspace(-0.5, 0.5, len(sI))

    #plot
    fig = plt.figure(figsize=(20, 12))
    fig.subplots_adjust(hspace=0)

    ax_1 = plt.subplot2grid((4,1),(0,0), colspan=1, rowspan=3)
    ax_1.tick_params(labelsize=14)
    ax_1.set_xticks([])
    ax_1.set_title(title, fontsize=36)
    ax_1.set_ylabel("Intensity", fontsize=20)
    ax_1.set_xlim(-0.5, 0.5)

    ax_2 = plt.subplot2grid((4,1),(3,0), colspan=1, rowspan=1)
    ax_2.tick_params(labelsize=14)
    ax_2.set_yticks([-0.5, 0.5])
    ax_2.set_xticks(np.linspace(-0.5, 0.5, 11))
    ax_2.set_xlabel("Pulse Phase", fontsize=20)
    ax_2.set_ylabel("Position Angle", fontsize=20)
    ax_2.set_xlim(-0.5, 0.5)

    ax_1.plot(x, sI, color="k", label="Stokes I")
    ax_1.plot(x, lin_pol, color="r", label="Linear Polarization")
    ax_1.plot(x, sV, color="b", label="Circular Polarization")
    ax_2.scatter(x, pa, color="k", label="Position Angle")

    ax_1.legend(loc="upper right", fontsize=18)
    fig_path = os.path.join(out_dir, save_name)
    logger.info("Saving polarimetry figure: {0}".format(fig_path))
    plt.savefig(fig_path)

    return fig_path

#--------------------------------------------------------------------------
def plot_stack(frequencies, profs_x, profs_y, pulsar_name,\
                out_dir="./", mybuffer=0.75, ignore_duplicates=True, special_freqs=None, ignore_freqs=None):
    """
    Plots multiple profiles stacked on top of one anothre in order of frequency. Saves as a .png

    Parameters:
    -----------
    frequencies: list
        The frequencies of the profiles to plot
    profs_x: list
        The phases of the profiles to plot. Where one phsae rotation begins as -0.5 and ends at 0.5
    profs_y: list
        The Intesities of the profiles
    pulsar_name: string
        The name of the pulsar
    out_dir: string
        OPTIONAL - The directory to output the .png to. Default: './'
    mybuffer: float
        OPTIONAL - The separation in intensity between profiles. Default:0.75
    ignore_dupicates: boolean
        OPTIONAL - If True, will not plot duplicate frequencies. Default: True
    special_freqs: list
        OPTIONAL - Any frequencies to be highlighted in the plot. Default: None
    ignore_freqs: list
        OPTIONAL - Any frequencies to not plot. Default: None

    Returns:
    --------
    fig_name: string
        The path of the saved .png
    """
    #initialize nones
    if not special_freqs:
        special_freqs = []
    if not ignore_freqs:
        ignore_freqs = []

    #Find and remove unwanted data
    ignore_idxs = []
    for i, freq in enumerate(frequencies):
        if freq in ignore_freqs:
            ignore_idxs.append(i)
        if ignore_duplicates and frequencies[i-1]==freq:
            ignore_idxs.append(i)
        ignore_idxs = list(set(ignore_idxs))

    for i in sorted(ignore_idxs, reverse=True):
        del frequencies[i]
        del profs_x[i]
        del profs_y[i]

    #roll the profiles to align the first maxima
    rolled_profs = []
    for profile in profs_y:
        _, _, _, _, _, _, maxima, _ = prof_utils.prof_eval_gfit(I)
        rl_prof = roll_data(I, idx_to_roll=maxima[0])[-1]
        rolled_profs.append(rl_prof)

    #Initialize figure
    plt.figure(figsize=(24, 20 + 2*len(frequencies)))
    #Loop over all frequencies
    for i, freq in enumerate(frequencies):
        if freq in special_freqs:
            clr = "red"
        else:
            clr = "black"

        x = profs_x[i]
        y = np.array(rolled_profs[i])/max(rolled_profs[i])
        y = np.array(y) + mybuffer*i
        plt.plot(x, y, color=clr)
        plt.text(0.35, 0.2+mybuffer*i, "{}MHz".format(round(freq, 2)), fontsize = 30, color = "black")

    #Finalizing the plot and saving the figure
    plt.xlim(-0.5, 0.5)
    plt.yticks([])
    plt.xticks(fontsize=30)
    plt.xlabel("Pulse Phase", fontsize=40)
    plt.ylabel("Intensity", fontsize=40)
    plt.title(pulsar_name + " Pulse Profiles", fontsize=60)
    fig_name = os.path.join(out_dir, pulsar_name + "_stacked_profiles.png")
    logger.info("Saving stacked profiles: {}".format(fig_name))
    plt.savefig(fig_name)

    return fig_name

def plot_stack_pol(frequencies, I_x, I_y, lin_y, circ_y, pulsar_name,\
                out_dir="./", mybuffer=1.1, ignore_duplicates=True, ignore_freqs=None):
    """
    Plots multiple profiles stacked on top of one anothre in order of frequency. Saves as a .png

    Parameters:
    -----------

    out_dir: string
        OPTIONAL - The directory to output the .png to. Default: './'
    mybuffer: float
        OPTIONAL - The separation in intensity between profiles. Default:0.75
    ignore_dupicates: boolean
        OPTIONAL - If True, will not plot duplicate frequencies. Default: True
    special_freqs: list
        OPTIONAL - Any frequencies to be highlighted in the plot. Default: None
    ignore_freqs: list
        OPTIONAL - Any frequencies to not plot. Default: None

    Returns:
    --------
    fig_name: string
        The path of the saved .png
    """
    #initialize nones
    if not ignore_freqs:
        ignore_freqs = []

    #Find and remove unwanted data
    ignore_idxs = []
    for i, freq in enumerate(frequencies):
        if freq in ignore_freqs:
            ignore_idxs.append(i)
        if ignore_duplicates and frequencies[i-1]==freq:
            ignore_idxs.append(i)
        ignore_idxs = list(set(ignore_idxs))

    for i in sorted(ignore_idxs, reverse=True):
        del frequencies[i]
        del Ix[i]
        del Iy[i]
        del lin_x[i]
        del lin_y[i]
        del circ_x[i]
        del circ_y[i]

    #roll the profiles to align the first maxima
    rolled_I = []
    rolled_lin = []
    rolled_circ = []
    for I, lin, circ in zip(I_y, lin_y, circ_y):
        _, _, _, _, _, _, maxima, _ = prof_utils.prof_eval_gfit(I)
        idx, roll_to, new_I = roll_data(I, idx_to_roll=maxima[0])
        new_lin = roll_data(lin, idx_to_roll=idx, roll_to=roll_to)[-1]
        new_circ = roll_data(circ, idx_to_roll=idx, roll_to=roll_to)[-1]
        rolled_I.append(new_I)
        rolled_lin.append(new_lin)
        rolled_circ.append(new_circ)

    #Initialize figure
    plt.figure(figsize=(24, 20 + 2*len(frequencies)))
    #Loop over all frequencies
    for i, freq in enumerate(frequencies):
        max_I = max(rolled_I[i])
        I_y = np.array(rolled_I[i])/max_I + mybuffer*i
        lin_y = np.array(rolled_lin[i])/max_I + mybuffer*i
        circ_y = np.array(rolled_circ[i])/max_I + mybuffer*i
        #logger.info("I_x: {}".format(I_x))
        #logger.info("I_y: {}".format(I_y))
        plt.plot(I_x[i], I_y, color="black")
        plt.plot(I_x[i], lin_y, color="red")
        plt.plot(I_x[i], circ_y, color="blue")
        plt.text(0.35, 0.2+mybuffer*i, "{}MHz".format(round(freq, 2)), fontsize = 30, color = "blue")

    #Finalizing the plot and saving the figure
    plt.xlim(-0.5, 0.5)
    plt.yticks([])
    plt.xticks(fontsize=30)
    plt.xlabel("Pulse Phase", fontsize=40)
    plt.ylabel("Intensity", fontsize=40)
    plt.title(pulsar_name + " Pulse Profiles", fontsize=60)
    fig_name = os.path.join(out_dir, pulsar_name + "_stacked_profiles_pol.png")
    logger.info("Saving stacked profiles: {}".format(fig_name))
    plt.savefig(fig_name)

    return fig_name

#--------------------------------------------------------------------------
if __name__ == '__main__':

    #dictionary for choosing log-levels
    loglevels = dict(DEBUG=logging.DEBUG,
                    INFO=logging.INFO,
                    WARNING=logging.WARNING,
                    ERROR = logging.ERROR)

    #Arguments
    parser = argparse.ArgumentParser(description="""Plots pulse profiles for a given pulsar""")

    obsop = parser.add_argument_group("Observation Options")
    obsop.add_argument("-p ", "--pulsar", type=str, help="J name of the pulsar. eg. J2241-5326")
    obsop.add_argument("-o", "--obsid", type=int, help="The MWA observation ID")
    obsop.add_argument("-f", "--freq", type=float, help="The observing frequency in MHz")

    ioop = parser.add_argument_group("Input and Output Opttions")
    ioop.add_argument("-b", "--bestprof", type=str, help="Location of the MWA bestprof file.")
    ioop.add_argument("-a", "--ascii", type=str, help="location of the dspsr RM fixed archive file in ascii format.")
    ioop.add_argument("-d", "--out_dir", type=str, default="./", help="Directory for output figure(s)")

    modeop = parser.add_argument_group("Modes")
    modeop.add_argument("--plt_bestprof", action="store_true", help="Plot a bestprof profile")
    modeop.add_argument("--plt_ascii", action="store_true", help="Plot an ascii profile")
    modeop.add_argument("--plt_pol", action="store_true", help="Plot a polarimetry profile from a supplied ascii archive")
    modeop.add_argument("--plt_stack", action="store_true", help="Plot data from epndb")
    modeop.add_argument("--plt_stack_pol", action="store_true", help="Plot data from epndb with full polarisation information")
    modeop.add_argument("--plt_bp_stack", action="store_true", help="Plot data from epndb and include supplied bestprof")
    modeop.add_argument("--plt_ascii_stack", action="store_true", help="Plot data from epndb and include supplied ascii file")
    modeop.add_argument("--plt_ascii_stack_pol", action="store_true", help="Plot data from epndb with full polarisation information\
                        and include supplied ascii file")

    otherop = parser.add_argument_group("Other Options")
    otherop.add_argument("-L", "--loglvl", type=str, help="Logger verbosity level. Default: INFO", choices=loglevels.keys(), default="INFO")

    args = parser.parse_args()

    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    ch.setLevel(loglevels[args.loglvl])
    formatter = logging.Formatter('%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    #Assertions
    if args.plt_bestprof:
        if not args.bestprof:
            logger.error("Please supply a bestprof profile to plot")

    if args.plt_ascii or args.plt_pol:
        if not args.ascii:
            logger.error("Please supply an ascii profile to plot")

    if args.plt_stack or args.plt_stack_pol:
        if not args.pulsar:
            logger.error("Please supply a pulsar name")

    if args.plt_bp_stack:
        if not args.pulsar or not args.bestprof or not args.freq:
            logger.error("Please ensure you have suppled a pulsar name as well as a bestprof profile")

    if args.plt_ascii_stack or args.plt_ascii_stack_pol:
        if not args.pulsar or not args.ascii or not args.freq:
            logger.error("Please ensure you have suppled a pulsar name, frequency and an ascii profile")

    #Do the things
    if args.plt_bestprof:
        plot_bestprof(args.bestprof, out_dir=args.out_dir)

    if args.plt_ascii:
        plot_ascii(args.ascii, pulsar=args.pulsar, freq=args.freq, obsid=args.obsid, out_dir=args.out_dir)

    if args.plt_pol:
        plot_archive_stokes(args.ascii, pulsar=args.pulsar, freq=args.freq, obsid=args.obsid, out_dir=args.out_dir)

    if args.plt_stack:
        pulsar_dict = get_data_from_epndb(args.pulsar)
        plot_stack(pulsar_dict["freq"], pulsar_dict["Ix"], pulsar_dict["Iy"], args.pulsar,\
            out_dir=args.out_dir, special_freqs=[args.freq])

    if args.plt_bp_stack or args.plt_ascii_stack:
        if args.plt_bp_stack:
            y = prof_utils.get_from_bestprof(args.bestprof)[-2]
        if args.plt_ascii_stack:
            y = prof_utils.get_from_ascii(args.ascii)[0]
        x = np.linspace(-0.5, 0.5, len(y))
        pulsar_dict = get_data_from_epndb(args.pulsar)
        for key in pulsar_dict.keys():
            if key == "Ix":
                pulsar_dict[key].append(x)
            elif key == "Iy":
                pulsar_dict[key].append(y)
            elif key == "freq":
                pulsar_dict[key].append(args.freq)
            else:
                pulsar_dict[key].append(None)

        #sort by frequency
        freq_list = pulsar_dict["freq"]
        for key in pulsar_dict.keys():
            _, pulsar_dict[key] = (list(t) for t in zip(*sorted(zip(freq_list, pulsar_dict[key]))))
        #plot
        plot_stack(pulsar_dict["freq"], pulsar_dict["Ix"], pulsar_dict["Iy"], args.pulsar,\
            out_dir=args.out_dir, special_freqs=[args.freq])

    if args.plt_stack_pol:
        pulsar_dict = get_data_from_epndb(args.pulsar)
        #remove anything without pol. info
        del_idxs=[]
        for i, Q, U, V in zip(range(len(pulsar_dict["Qy"])), pulsar_dict["Qy"], pulsar_dict["Uy"], pulsar_dict["Vy"]):
            if not Q or not U or not V:
                del_idxs.append(i)
        for i in sorted(del_idxs, reverse=True):
            del pulsar_dict["Ix"][i]
            del pulsar_dict["Iy"][i]
            del pulsar_dict["Qy"][i]
            del pulsar_dict["Uy"][i]
            del pulsar_dict["Vy"][i]
            del pulsar_dict["freq"][i]

        #calc lin pol
        lin=[]
        for Qy, Uy in zip(pulsar_dict["Qy"], pulsar_dict["Uy"]):
            lin.append(calc_lin_pa(Qy, Uy)[0])

        #plot
        plot_stack_pol(pulsar_dict["freq"], pulsar_dict["Ix"], pulsar_dict["Iy"], lin, pulsar_dict["Vy"], args.pulsar,\
            out_dir=args.out_dir)

    if args.plt_ascii_stack_pol:
        pulsar_dict = get_data_from_epndb(args.pulsar)
        I, Q, U, V, length = prof_utils.get_stokes_from_ascii(args.ascii)
        x = np.linspace(-0.5, 0.5, length)
        for key in pulsar_dict.keys():
            if key == "Ix":
                pulsar_dict[key].append(x)
            elif key == "Iy":
                pulsar_dict[key].append(I)
            elif key == "Qy":
                pulsar_dict[key].append(Q)
            elif key == "Uy":
                pulsar_dict[key].append(U)
            elif key == "Vy":
                pulsar_dict[key].append(V)
            elif key == "freq":
                pulsar_dict[key].append(args.freq)
            else:
                pulsar_dict[key].append(None)

        #remove anything without pol. info
        del_idxs=[]
        for i, Q, U, V in zip(range(len(pulsar_dict["Qy"])), pulsar_dict["Qy"], pulsar_dict["Uy"], pulsar_dict["Vy"]):
            if not Q or not U or not V:
                del_idxs.append(i)
        for i in sorted(del_idxs, reverse=True):
            del pulsar_dict["Ix"][i]
            del pulsar_dict["Iy"][i]
            del pulsar_dict["Qy"][i]
            del pulsar_dict["Uy"][i]
            del pulsar_dict["Vy"][i]
            del pulsar_dict["freq"][i]

        #sort by frequency
        freq_list = pulsar_dict["freq"]
        for key in pulsar_dict.keys():
            _, pulsar_dict[key] = (list(t) for t in zip(*sorted(zip(freq_list, pulsar_dict[key]))))

        #calc lin pol
        lin=[]
        for Qy, Uy in zip(pulsar_dict["Qy"], pulsar_dict["Uy"]):
            lin.append(calc_lin_pa(Qy, Uy)[0])

        #plot
        plot_stack_pol(pulsar_dict["freq"], pulsar_dict["Ix"], pulsar_dict["Iy"], lin, pulsar_dict["Vy"], args.pulsar,\
            out_dir=args.out_dir)