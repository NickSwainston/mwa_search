#!/usr/bin/env python3

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import logging
import argparse
import json
import time

import search_epndb
import binfinder
import mwa_metadb_utils
import os
from data_process_pipeline import run_params_class

logger = logging.getLogger(__name__)


#--------------------------------------------------------------------------
def align_data_with_peak(stokes_I, stokes_Q=None, stokes_U=None, stokes_V=None):
    #feed this only the y values of each dataset
    maxval = 0
    max_i = 0
    for i, point in enumerate(stokes_I):
        if point > maxval:
            maxval = point
            max_i = i
    if max_i == 0 and maxval == 0:
        logger.warn("data has not been aligned properly. Your data may be faulty")
    #This is the amount to cycle the data to align the max value with the centre
    roll_n = len(stokes_I) + len(stokes_I)//2 - int(max_i)
    stokes_I = np.roll(stokes_I, roll_n)
    if stokes_Q is not None:
        stokes_Q = np.roll(stokes_Q, roll_n)
    if stokes_U is not None:
        stokes_U = np.roll(stokes_U, roll_n)
    if stokes_V is not None:
        stokes_V = np.roll(stokes_V, roll_n)

    if stokes_Q is not None:
        return stokes_I, stokes_Q, stokes_U, stokes_V
    else:
        return stokes_I

#--------------------------------------------------------------------------
def normalize(stokes_I, stokes_Q=None, stokes_U=None, stokes_V=None, maxval=None):

    #allows for normalization wrt a custom max value
    if maxval==None:
        maxval = 0
        for i in stokes_I:
            if abs(i)>maxval:
                maxval = abs(i)

    for i, _ in enumerate(stokes_I):
        if maxval <= 0:
            logger.warn("Division by zero or a negative number")
        stokes_I[i] = stokes_I[i]/maxval
        if stokes_Q:
            stokes_Q[i] = stokes_Q[i]/maxval
        if stokes_U:
            stokes_U[i] = stokes_U[i]/maxval
        if stokes_V:
            stokes_V[i] = stokes_V[i]/maxval

    if stokes_Q:
        return stokes_I, stokes_Q, stokes_U, stokes_V
    else:
        return stokes_I

#--------------------------------------------------------------------------
def add_buffer(data, mybuffer):
    for i, _ in enumerate(data):
        data[i] = data[i] + mybuffer
    return data

def split_list(alist):
    x = []
    y = []
    for i, x_y in enumerate(alist[0]):
        x.append(x_y[i])
        y.append(x_y[i])

    return x, y

#--------------------------------------------------------------------------
def calc_lin_pa(stokes_Q, stokes_U):
    lin_pol = []
    pa = []
    for q, u in zip(stokes_Q, stokes_U):
        lin_pol.append(np.sqrt(q**2+u**2))
        pa.append(np.arctan(q/u))
    return lin_pol, pa

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


    return frequency, stokes_I, stokes_Q, stokes_U, stokes_V


#--------------------------------------------------------------------------
def plot_bestprof(bestprof, out_dir, nocrop=False):

    #retrieve data from bestprof
    x = []
    y = []
    print("Plotting profile from file: {0}".format(bestprof))
    f = open(bestprof)
    lines = f.readlines()
    for line in lines:
        line=line.split()
        if "#" not in line[0]:
            x.append(float(line[0]))
            y.append(float(line[1]))
    f.close()


    y = normalize(y)
    y = align_data_with_peak(y)
    x = normalize(x)

    info_dict = binfinder.bestprof_info(filename=bestprof)

    x_len=12
    y_len=8
    fig, ax = plt.subplots(figsize=(x_len, y_len))
    #Crop the image
    if nocrop==False:
        ax.set_xlim(0.25, 0.75)
    else:
        ax.set_xlim(0., 1.)

    plt.title("{0}_{1}_Profile".format(info_dict["obsid"], info_dict["pulsar"]))
    plt.text(0.05, 0.95,    "S/N:             {0}".format(info_dict["sn"]), fontsize=9, color="black", transform=ax.transAxes)
    plt.text(0.05, 0.925,     "Chi Sq:          {0}".format(info_dict["chi"]), fontsize=9, color="black", transform=ax.transAxes)
    plt.text(0.05, 0.9,    "DM:              {0}".format(info_dict["dm"]), fontsize=9, color="black", transform=ax.transAxes)
    plt.text(0.05, 0.875,     "Period (ms):     {0}+/-{1}".format(info_dict["period"], info_dict["period_error"]), fontsize=9, color="black", transform=ax.transAxes)

    ax.plot(x, y, color="black")
    fig_name = "{0}_{1}_presto_pulse_prof.png".format(info_dict["obsid"], info_dict["pulsar"]) 
    fig_path = os.path.join(out_dir, fig_name)
    print("Saving figure:   {0}".format(fig_path))
    plt.savefig("{0}".format(fig_path))

    return fig_path




#--------------------------------------------------------------------------
def plot_archive(archive, obsid, pulsar, freq, out_dir="./", nocrop=False):
    #Read the archive
    x=[]
    sI=[]
    sQ=[]
    sU=[]
    sV=[]
    f = open(archive)
    lines = iter(f.readlines())
    next(lines) #skip first line
    for line in lines:
        thisline=line.split()
        x.append(int(float(thisline[2])))
        sI.append(float(thisline[3]))
        sQ.append(float(thisline[4]))
        sU.append(float(thisline[5]))
        sV.append(float(thisline[6]))

    #normalize, align and find linear polarization and position angle
    x=normalize(x)
    sI, sQ, sU, sV = normalize(sI, stokes_Q=sQ, stokes_U=sU, stokes_V=sV)
    sI, sQ, sU, sV = align_data_with_peak(sI, stokes_Q=sQ, stokes_U=sU, stokes_V=sV)
    lin_pol, pa = calc_lin_pa(sQ, sU)

    #plot -
    fig = plt.figure(figsize=(20, 12))
    fig.subplots_adjust(hspace=0)

    if nocrop==False:
        crop=(0.25, 0.75)
    else:
        crop=(0,1)

    ax_1 = plt.subplot2grid((4,1),(0,0), colspan=1, rowspan=3)
    ax_1.tick_params(labelsize=14)
    ax_1.set_xticks([])#empty
    ax_1.set_title("{0} Pulse Profile - {1} MHz".format(pulsar, freq), fontsize=36)
    ax_1.set_ylabel("Amplitude", fontsize=20)
    ax_1.set_xlim(crop)

    ax_2 = plt.subplot2grid((4,1),(3,0), colspan=1, rowspan=1)
    ax_2.tick_params(labelsize=14)
    ax_2.set_xlim(crop)
    ax_2.set_yticks([-1.5,-1.0,-0.5,0,0.5,1.0,1.5])
    ax_2.set_xlabel("Pulse Phase", fontsize=20)
    ax_2.set_ylabel("Position Angle", fontsize=20)

    ax_1.plot(x, sI, color="k", label="Stokes I")
    if lin_pol is not None:
        ax_1.plot(x, lin_pol, color="r", label="Linear Polarization")
    if sV is not None:
        ax_1.plot(x, sV, color="b", label="Circular Polarization")
    if pa is not None:
        ax_2.scatter(x, pa, color="k", label="Position Angle")

    ax_1.legend(loc="upper right", fontsize=18)
    plt.savefig("{0}/{1}_polarimetry_profile.png".format(out_dir, pulsar))

    return "{0}/{1}_polarimetry_profile.png".format(out_dir, pulsar)

#--------------------------------------------------------------------------
def plot_stack(frequencies, stokes_I, lin_pol, stokes_V, pulsar_name, out_dir, mybuffer=0.75):

    #Initialize figure
    plt.figure(figsize=(24, 20 + 4*len(data_dict)))
    #Loop over all frequencies
    for i, freq in enumerate(frequencies):
        #roll_n = 0 #safety valve
        #maxval = 0

        #Stokes I:
        I_x = stokes_I[i][0]
        I_y = add_buffer(stokes_I[i][1], mybuffer*i)
        plt.plot(I_x, I_y, ["k"])
        plt.text(0.35, 0.2+buffer*i, freq+" MHz", fontsize = 40, color = "red")

        #linear Polarization:
        if lin_pol[i][0]:
            lin_pol_x = lin_pol[i][0]
            lin_pol_y = add_buffer(lin_pol[i][1], mybuffer*i)
            plt.plot(lin_pol_x, lin_pol_y, ["b"])

        #Circular Polarization:
        if stokes_V[i][0]:
            V_x = stokes_V[i][0]
            V_y = add_buffer(stokes_V[i][1], my_buffer*i)
            plt.plot(V_x, V_y, ["r"])


    #Finalizing the plot and saving the figure
    plt.yticks([])
    custom_lines = [Line2D([0], [0], color="k", lw=2),
                   Line2D([0], [0], color="b", lw=2),
                   Line2D([0], [0], color="r", lw=2)]
    plt.legend(custom_lines, ["Stokes I", "Linear Polarization", "Circular Polarization"], loc="best", fontsize=40)
    plt.title(pulsar_name + " Pulse Profiles", fontsize=60)
    plt.savefig(out_dir + pulsar_name + "_stacked_pulse_priofiles.png")


#--------------------------------------------------------------------------
def sort_data(pulsar_name, epndb_dir, out_dir, prof_path=None, full_stokes=False):

    #Plot MWA profile if requested
    if prof_path is not None:
        logger.info("Plotting MWA single profile")
        #mwa_x, mwa_y = plot_MWA_profile(out_dir, pulsar_name, prof_path)

    #Grab the dictionary of epn puslars
    epn_dict = search_epndb.get_epn_paths(pulsar_name, epndb_dir)
    paths = epn_dict[pulsar_name]
    json_files = []
    for file in paths:
        logger.debug("Reading data from {0}".format(file))
        f = open(file, "r")
        json_files.append(json.load(f))
        f.close()
    logger.info("Found {0} different profiles on the EPN database for this pulsar".format(len(paths)))

     #Each json file is a dictionary that points to 3 dictionaries described as follows:
    #The high level dictionary has keys {"hdr", "files", "series"} - each of these keys point to a dictionary
    #"hdr" contains all of the relevant header information as follows:
        #'name', 'freq', 'site', 'rm', 'dm', 'scale', 'state', 'npol', 'length', 'rcvr:name', 'be:name', 'license', 'ref', 'basename'
    #"files" is pretty much useless, don't worry about it
    #"series" contains all of the raw data split into stokes params:
        #'I', 'Q', 'U', 'V', 'PA', 'PAE', 'L', 'P', 'EL', 'ELE', 'max', 'min'
        #each data point in a stokes series is an array in the format: [x_val, y_val]

    #Some json files don't have frequency data - this makes them useless to us:
    for fileno, data_dict in enumerate(json_files):
        newdict = json_files
        if "freq" not in data_dict['hdr']:
            logger.warn("No frequency information found for file {0} This data will be ignored!".format(paths[fileno]))
            del newdict[fileno]
        json_files = newdict
        del newdict

    #TODO: get MWA profile and append it to the json_files dictionary


    frequencies = []
    stokes_I = []
    stokes_Q = []
    stokes_U = []
    stokes_V = []
    #Now put everything in a less complex (stupid) format
    for i, data_dict in enumerate(json_files):
        f = (data_dict["hdr"]["freq"])
        I_x, I_y = split_list(data_dict["series"]["I"])

        if data_dict["series"]["Q"]:
            Q_x, Q_y = split_list(data_dict["series"]["Q"])
        else:
            Q_x = []
            Q_y = []
        if data_dict["series"]["U"]:
            U_x, U_y = split_list(data_dict["series"]["U"])
        else:
            U_x = []
            U_y = []
        if data_dict["series"]["U"]:
            V_x, V_y = split_list(data_dict["series"]["V"])
        else:
            V_x = []
            V_y = []

        frequencies.append(f)
        stokes_I.append([I_x, I_y])
        stokes_Q.append([Q_x, Q_y])
        stokes_U.append([U_x, U_y])
        stokes_V.append([V_x, V_y])

    #Stokes_I looks like:[[I_x_1, I_y_1][I_x_2, I_y_2]]
    #where [I_x, I_y] looks like: [[x_1, x_2, etc][y_1, y_2, etc]]
    #Meaning Stokes_I[i][0] = [I_x]
    #Stokes_I[i][1] = [I_y]
    #which corresponds to the freqeuncy at frequencies[i]

    #order the datasets based on frequency if there are more than one:
    if len(frequencies)>1:
        f = []
        for i in frequecies:
            f.append(str(i))
        frequencies, stokes_I, stokes_Q, stokes_U, stokes_V = zip(sorted*(zip(frequencies, stokes_I, stokes_Q, stokes_U, stokes_V)))


    #Align the data with the peak value in stokes I
    lin_pol = []
    pa = []
    for i in range(len(frequencies)):
        stokes_I[i][1], stokes_Q[i][1], stokes_U[i][1], stokes_V[i][1] = align_data_with_peak(stokes_I[i][1], stokes_Q[i][1], stokes_U[i][1], stokes_V[i][1])
        #Normalize:
        stokes_I, stokes_Q, stokes_U, stokes_V = normalize(stokes_I[i][1], stokes_Q[i][1], stokes_U[i][1], stokes_V[i][1])
        #Now get linear Pol and Position Angle
        if stokes_Q[i] and stokes_U[i]:
            lin_pol[i], pa[i] = calc_lin_pa(stokes_Q[i], stokes_U[i])
        else:
            lin_pol[i] = [[],[]]
            pa[i] = [[],[]]


    #find the mwa data and plot it
    for i, freq in frequencies:
        if mwa_freq == freq:
            #Normalize range
            stokes_I = normalize(stokes_I[i][0])
            for j in stokes_I[i][0]:
                stokes_I[i][0][j] = stokes_I[i][0][j]-0.5
            #Set all the other stuff to the same axis
            lin_pol[i][0] = stokes_I[i][0]
            pa[i][0] = stokes_I[i][0]
            stokes_V[i][0] = stokes_I[i][0]

            #Now it's ready to plot
            plot_mwa_profile(stokes_I[i], lin_pol[i], pa[i], stokes_V[i], freq, pulsar_name, out_dir)
            break

    #Plot the stacked profiles
    plot_stacked(frequencies, stokes_I, lin_pol, stokes_V, pulsar_name, out_dir)


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
    obsop.add_argument("-p ", "--pulsar", type=str, help="J name of the pulsar. eg. J2241-5326. Default: ''")
    obsop.add_argument("-o", "--obsid", type=str, default="", help="The observation ID. Default: ''")
    obsop.add_argument("-O", "--cal_id", type=str, default="", help="The ID of the calibrator. Default: ''")

    ioop = parser.add_argument_group("Input and Output Opttions")
    ioop.add_argument("-b", "--bestprof", type=str, help="Location of the MWA bestprof file.")
    ioop.add_argument("-a", "--archive", type=str, help="location of the dspsr RM fixed archive file in ascii format.")
    ioop.add_argument("-d", "--out_dir", type=str, help="Directory for output figure(s)")
    ioop.add_argument("--epndb_dir", type=str, default="/group/mwaops/k_smith/www.epta.eu.org/epndb/json", help="location of the epn database json folder. Default: /group/mwaops/k_smith/www.epta.eu.org/epndb/json")

    otherop = parser.add_argument_group("Other Options")
    otherop.add_argument("-L", "--loglvl", type=str, help="Logger verbosity level. Default: INFO", choices=loglevels.keys(), default="INFO")
    otherop.add_argument("--nocrop", action="store_true", help="Use this tag to plot the full pulse profile, otherwise the image will be cropped by 25% on each side")

    modeop = parser.add_argument_group("Modes")
    parser.add_argument("-m", "--mode", type=str, help="The desired plotting mode. Options:\n\
                                                    'b' - plot the bestprof profile\n\
                                                    's' - plot a dspsr ascii file with full stokes.")

    args = parser.parse_args()


    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    ch.setLevel(loglevels[args.loglvl])
    formatter = logging.Formatter('%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    if args.out_dir is None:
        logger.error("Please supply an output directory and rerun.")
        sys.exit(1)

    run_params = run_params_class(pulsar=args.pulsar, obsid=args.obsid, cal_id=args.cal_id,\
                    bestprof=args.bestprof, archive=args.archive, out_dir=args.out_dir,\
                    epndb_dir=args.epndb_dir, loglvl=args.loglvl, mode=args.mode, nocrop=args.nocrop)


    #sort_data(args.pulsar_name, args.epndb_dir, args.out_dir, args.bestprof, args.full_stokes)
    if run_params.mode=="b":
        plot_bestprof(args.bestprof, args.out_dir, args.nocrop)
    elif run_params.mode=="s":
        plot_archive(run_params=run_params)
