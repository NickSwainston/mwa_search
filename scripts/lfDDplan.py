#! /usr/bin/env python3

import argparse
from mwa_search.dispersion_tools import plot_sensitivity, dd_plan
from vcstools.metadb_utils import get_common_obs_metadata


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
                        help='The minimun DM step size, default 0.02')
    parser.add_argument('--max_DM_step', type=float, default=500.0,
                        help='The maximum DM step size, default 500.0')
    parser.add_argument('-p', '--plot', action='store_true',
                        help='Plot the sensitivty of the DM plan')
    parser.add_argument('--time', type=int, default=4800,
                        help='Time in seconds to calculate the sensitivity')
    parser.add_argument('--max_dms_per_job', type=int, default=5000,
                        help='If Nsteps is greater than this value split it into multiple lines. '
                             'This will cause the search pipeline to submit fewer DMs per job')
    #parser.add_argument()
    args=parser.parse_args()

    if args.obsid:
        #get the centrefreq from the obsid metadata
        beam_meta_data = get_common_obs_metadata(args.obsid)
        obs, ra, dec, dura, [xdelays, ydelays], centrefreq, channels = beam_meta_data

        args.centrefreq = centrefreq


    DD_plan_array = dd_plan(
        args.centrefreq,
        args.bandwidth,
        args.nfreqchan,
        args.timeres,
        args.lowDM,
        args.highDM,
        min_DM_step=args.min_DM_step,
        max_DM_step=args.max_DM_step,
        max_dms_per_job=args.max_dms_per_job
    )
    print(" low DM | high DM | DeltaDM | Nsteps | Downsamp | nsub | workfactor | Effective time resolution (ms) ")
    total_steps = 0
    for dd_line in DD_plan_array:
        dm_min, dm_max, dm_step, ndm, timeres, downsamp, nsub, total_work_factor = dd_line
        print(f"{dm_min:7.1f} | {dm_max:7.1f} | {dm_step:7.2f} | {ndm:6d} | {downsamp:8d} | {nsub:4d} | {total_work_factor:10.2f} | {timeres:7.3f}")
        total_steps += ndm
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
