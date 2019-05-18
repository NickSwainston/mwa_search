#! /usr/bin/env python

import argparse
import math


def dd_plan(centrefreq, bandwidth, nfreqchan, timeres, lowDM, highDM):
    
    DD_plan_array = []
    freqres = bandwidth / float(nfreqchan)
    previous_DM = lowDM
    #Loop until you've made a hit your range max
    D_DM = 0.
    while D_DM < highDM:
        #calculate the DM where the current time resolution equals the 
        #dispersion in a frequency channel (a bit of an overkill)
        D_DM = timeres * centrefreq**3 /\
               (8.3 * 10.**6 * freqres)
        
        #difference in DM that will double the effective width (eq 6.4 of pulsar handbook)
        #TODO make this more robust
        DM_step = math.sqrt( (2.*timeres)**2 - timeres**2 )/\
                  (8.3 * 10**6 * bandwidth / centrefreq**3)

        #round to nearest 0.01
        DM_step = round(DM_step, 2)
        if DM_step == 0.0:
            #set DM to 0.01 as a zero DM doesn't make sense
            DM_step = 0.01

        if D_DM > highDM:
            #last one so range from to max
            nDM_step = int((highDM - previous_DM) / DM_step)
            DD_plan_array.append([ previous_DM, highDM, DM_step, nDM_step, timeres ])
            previous_DM = D_DM
        elif D_DM > previous_DM:
            #range from last to new
            nDM_step = int((D_DM - previous_DM) / DM_step)
            DD_plan_array.append([ previous_DM, D_DM, DM_step, nDM_step, timeres ])
            previous_DM = D_DM
        #else don't do anything because it's under your low DM

        #Double time res to account for incoherent dedispersion
        timeres *= 2.

    return DD_plan_array



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    Used to calculate the Dedispersion plan for low-frequency telescopes. Inspired by PRESTO's DDplan.py.
    """)
    parser.add_argument('-f', '--centrefreq', type=float, help='Centre frequency of the observation in MHz.')
    parser.add_argument('-b', '--bandwidth', type=float, help='Bandwidth of the observation in MHz.')
    parser.add_argument('-nf', '--nfreqchan', type=int, help='Number of frequency channels.')
    parser.add_argument('-t', '--timeres', type=float, help='Time resolution in ms.')
    parser.add_argument('-ld', '--lowDM', type=float, help='Lowest DM of the required range.')
    parser.add_argument('-hd', '--highDM', type=float, help='Highest DM of the required range.')
    #parser.add_argument()
    args=parser.parse_args()
    
    DD_plan_array = dd_plan( args.centrefreq, args.bandwidth, args.nfreqchan, args.timeres, args.lowDM, args.highDM)
    print(" low DM | high DM | DeltaDM | Nsteps | Effective time resolution (ms)")
    total_steps = 0
    for d in DD_plan_array:
        print("{0:7.1f} | {1:7.1f} | {2:7.2f} | {3:6d} | {4:7.3f}".\
               format(d[0], d[1], d[2], d[3], d[4]))
        total_steps += d[3]
    print("Total DM steps required: {}".format(total_steps))
