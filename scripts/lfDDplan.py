#! /usr/bin/env python3

import argparse
import math
import numpy as np
import matplotlib.pyplot as plt

def plot_sensitivity(DD_plan_array, time, centrefreq, freqres, bandwidth):
    base_sensitivity = 3 #mJy. This could be done properly but this will for now
    # adjust for time
    base_sensitivity = base_sensitivity * math.sqrt(4800) / math.sqrt(time)
    # adjust for unscattered pulsar
    #base_sensitivity = base_sensitivity / math.sqrt( ( 1. - 0.05) / 0.05 )

    # period to work with
    periods = np.array([ 1., 0.1, 0.01, 0.001 ])*1000.
    widths = periods*0.05

    #fig, ax = plt.subplots(1, 1)
    plt.subplots(1, 1)
    #ax.set_ylim(0, 100)
    #plt.ylim(0, 100)

    for period, width in zip(periods, widths):
        sensitivties = []
        DMs = []
        for dm_row in DD_plan_array:
            DM_start, D_DM, DM_step, _, timeres = dm_row
            for DM in np.arange(DM_start, D_DM, DM_step):
                # For each DM you're going to search do a sensitivy cal based on how
                # a pulse will get

                #Dm smear over a frequency channel
                dm_smear = DM * freqres * 8.3 * 10.**6 / centrefreq**3
                #Dm smear due to maximum incorrect DM
                dm_step_smear = 8.3 * 10.**6 * DM_step / 2. * bandwidth / centrefreq**3

                effective_width = math.sqrt(width**2 + dm_smear**2 + dm_step_smear**2 + timeres**2)
                #print(effective_width)
                #sensitivity given new effectiv width
                if effective_width >= period:
                    sensitivties.append(1000.)
                else:
                    sensitivties.append(base_sensitivity /
                                    math.sqrt( ( period - effective_width) / effective_width ) *
                                    math.sqrt( ( period - width) / width ))
                DMs.append(DM)
        plt.plot(DMs, sensitivties, label="P={0} ms".format(period))
    #plt.yscale('log')
    plt.legend()
    plt.yscale('log')
    plt.xscale('log')
    plt.ylabel(r"Detection Sensitivity, 10$\sigma$ (mJy)")
    plt.xlabel(r"Dispersion measure (pc cm$^{-3}$ ")
    plt.title("Sensitivy using a minimum DM step size of {0}".format(DD_plan_array[0][2]))
    #plt.show()
    plt.savefig("DM_step_sens_mDMs_{0}.png".format(DD_plan_array[0][2]))


def calc_nsub(centrefreq, dm):
    #work out how many subbands to use based on the dm smear over a subband
    nsub = 2
    #time_res and dm_smear are in ms
    time_res = 0.1
    dm_smear = dm * 0.01 / nsub * 8.3 * 10.**6 / centrefreq**3
    while dm_smear > time_res:
        nsub *= 2.
        dm_smear = dm * 0.01 / nsub * 8.3 * 10.**6 / centrefreq**3
    return int(nsub)


def dd_plan(centrefreq, bandwidth, nfreqchan, timeres, lowDM, highDM, min_DM_step=0.02):
    """
    Work out the dedisperion plan

    Parameters
    ----------
    centrefreq: float
        The center frequency of the observation in MHz
    bandwidth: float
        The bandwidth of the observation in MHz
    nfreqchan: int
        The number of frequency channels
    timeres: float
        The time resolution of the observation in ms
    lowDM: float
        The lowest dispersion measure
    highDM: float
        The highest dispersion measure
    min_DM_step: float
        Will overwrite the minimum DM step with this value

    Returns
    -------
    DD_plan_array: list list
        dedispersion plan format:
        [[low_DM, high_DM, DM_step, nDM_step, timeres, downsample, nsub ]]
    """

    DD_plan_array = []
    freqres = bandwidth / float(nfreqchan)
    previous_DM = lowDM

    #number of time samples smeared over before moving to next D_dm
    smear_fact = 3.

    #Loop until you've made a hit your range max
    D_DM = 0.
    downsample = 1
    while D_DM < round(highDM, 2):
        #calculate the DM where the current time resolution equals the
        #dispersion in a frequency channel (a bit of an overkill)

        #Dm smear over a frequency channel
        dm_smear = previous_DM * freqres * 8.3 * 10.**6 / centrefreq**3
        total_smear = math.sqrt(timeres**2 +
                                dm_smear**2)


        D_DM = smear_fact * timeres * centrefreq**3 /\
               (8.3 * 10.**6 * freqres)

        #difference in DM that will double the effective width (eq 6.4 of pulsar handbook)
        #TODO make this more robust
        #DM_step = math.sqrt( (2.*timeres)**2 - timeres**2 )/\
        #          (8.3 * 10**6 * bandwidth / centrefreq**3)
        DM_step = smear_fact * total_smear * centrefreq**3 /\
                  (8.3 * 10.**6 * 0.5 * bandwidth)


        #round to nearest 0.01
        DM_step = round(DM_step, 2)
        if DM_step < min_DM_step:
            #set DM to 0.01 as a zero DM doesn't make sense
            DM_step = min_DM_step


        if D_DM > highDM:
            #last one so range from to max
            D_DM = highDM
        #range from last to new
        D_DM = round(D_DM, 2)
        nDM_step = int((D_DM - previous_DM) / DM_step)
        if D_DM > lowDM:
            nsub = calc_nsub(centrefreq, D_DM)
            DD_plan_array.append([ previous_DM, D_DM, DM_step, nDM_step, timeres, downsample, nsub ])
            previous_DM = D_DM

        #Double time res to account for incoherent dedispersion
        timeres *= 2.
        downsample *= 2

    return DD_plan_array



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    Used to calculate the Dedispersion plan for low-frequency telescopes. Inspired by PRESTO's DDplan.py.
    """)
    parser.add_argument('-f', '--centrefreq', type=float, default=150.,
                        help='Centre frequency of the observation in MHz.')
    parser.add_argument('-b', '--bandwidth', type=float, default=30.72,
                        help='Bandwidth of the observation in MHz.')
    parser.add_argument('-nf', '--nfreqchan', type=int, default=3072,
                        help='Number of frequency channels.')
    parser.add_argument('-t', '--timeres', type=float, default=0.1,
                        help='Time resolution in ms.')
    parser.add_argument('-ld', '--lowDM', type=float, default=1.,
                        help='Lowest DM of the required range.')
    parser.add_argument('-hd', '--highDM', type=float, default=250.,
                        help='Highest DM of the required range.')
    parser.add_argument('-o', '--obsid', type=int,
                        help='The MWA observation ID of an observation. Using this command will get the require observation parameters.')
    parser.add_argument('-m', '--min_DM_step', type=float, default=0.02,
                        help='The  minimun DM step size, default 0.02')
    parser.add_argument('-p', '--plot', action='store_true',
                        help='Plot the sensitivty of the DM plan')
    parser.add_argument('--time', type=int, default=4800,
                        help='Time in seconds to calculate the sensitivity')
    #parser.add_argument()
    args=parser.parse_args()

    if args.obsid:
        #get the centrefreq from the obsid metadata
        beam_meta_data = get_common_obs_metadata(args.obsid)
        obs, ra, dec, dura, [xdelays, ydelays], centrefreq, channels = beam_meta_data

        args.centrefreq = channels


    DD_plan_array = dd_plan( args.centrefreq, args.bandwidth, args.nfreqchan, args.timeres, args.lowDM, args.highDM, min_DM_step=args.min_DM_step)
    print(" low DM | high DM | DeltaDM | Nsteps | Downsamp | nsub | Effective time resolution (ms) ")
    total_steps = 0
    for d in DD_plan_array:
        print("{0:7.1f} | {1:7.1f} | {2:7.2f} | {3:6d} | {4:8d} | {5:4d} | {6:7.3f}".\
               format(d[0], d[1], d[2], d[3], d[5], d[6], d[4]))
        total_steps += d[3]
    print("Total DM steps required: {}".format(total_steps))

    if args.plot:
        #work out time to use
        if args.time:
            time = args.time
        elif args.obsid:
            time = dura
        else:
            #using default
            time = args.time
        freq_res = args.bandwidth / args.nfreqchan
        plot_sensitivity(DD_plan_array, time, args.centrefreq, freq_res, args.bandwidth)
