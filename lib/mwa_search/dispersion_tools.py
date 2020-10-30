import math
import numpy as np
from matplotlib import use
use('Agg')
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
            DM_start, D_DM, DM_step, _, timeres, _, _ = dm_row
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


def dd_plan(centrefreq, bandwidth, nfreqchan, timeres, lowDM, highDM,
            min_DM_step=0.02, max_DM_step=500.0, max_dms_per_job=5000):
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
    max_dms_per_job: int
        If Nsteps is greater than this value split it into multiple lines

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
        total_smear = math.sqrt(timeres**2 + dm_smear**2)


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
        if DM_step > max_DM_step:
            DM_step = max_DM_step


        if D_DM > highDM:
            #last one so range from to max
            D_DM = highDM
        #range from last to new
        D_DM = round(D_DM, 2)
        nDM_step = int((D_DM - previous_DM) / DM_step)
        if D_DM > lowDM:
            nsub = calc_nsub(centrefreq, D_DM)
            if downsample > 16:
                DD_plan_array.append([ previous_DM, D_DM, DM_step, nDM_step, timeres, 16, nsub ])
            else:
                DD_plan_array.append([ previous_DM, D_DM, DM_step, nDM_step, timeres, downsample, nsub ])

            previous_DM = D_DM

        #Double time res to account for incoherent dedispersion
        timeres *= 2.
        downsample *= 2

    #Check no lines have more Nsteps than max_dms_per_job
    new_DD_plan_array = []
    for dd_line in DD_plan_array:
        new_dd_lines = []
        while dd_line[3] > max_dms_per_job:
            # previous_DM, D_DM, DM_step, nDM_step, timeres, downsample, nsub
            new_dd_lines.append([dd_line[0], dd_line[0] + dd_line[2] * max_dms_per_job,
                                dd_line[2], max_dms_per_job, dd_line[4],
                                dd_line[5], dd_line[6]])
            dd_line = [dd_line[0] + dd_line[2] * max_dms_per_job, dd_line[1],
                       dd_line[2], dd_line[3] - max_dms_per_job, dd_line[4],
                       dd_line[5], dd_line[6]]
        new_dd_lines.append(dd_line)
        for n_line in new_dd_lines:
            new_DD_plan_array.append(n_line)

    return new_DD_plan_array