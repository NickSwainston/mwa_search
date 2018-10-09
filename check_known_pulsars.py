#! /usr/bin/env python

import argparse
import find_pulsar_in_obs as fpio
import os
import glob
import subprocess
import numpy as np
import blindsearch_pipeline as blind_pipe

def beamform_and_fold(obsid, DI_dir, all_check, args):
    
    
    cmd_line='file_maxmin.py {0}'.format(obsid)
    cmd=subprocess.Popen(cmd_line,shell=True,stdout=subprocess.PIPE)
    for line in cmd.stdout:
        if line.split()[0]=='First':
            obsbeg=int(line.split()[2])
        if line.split()[0]=='Last':
            obsend=int(line.split()[2])
    obsdur=obsend-obsbeg

    #wrapping for find_pulsar_in_obs.py
    names_ra_dec = np.array(fpio.grab_source_alog())
    fpio.find_sources_in_obs([obsid], names_ra_dec)
    known_pulsar_file = "{0}_analytic_beam.txt".format(obsid)

    if all_check:
        #looks through the comined files to use the max and min
        #TODO have some sort of check to look for gaps
        combined_files = glob.glob("/group/mwaops/vcs/{0}/combined/{0}*_ics.dat".format(obsid))
        comb_times = []
        for comb in combined_files:
            comb_times.append(int(comb.split("_")[1]))
        psrbeg = min(comb_times)
        psrend = max(comb_times)
    
    #looping through each puslar
    
    with open(known_pulsar_file) as kpf:
        for pulsar_line in kpf:
            if pulsar_line.startswith("J"):
                PSRJ = pulsar_line.split()[0]
                if len(PSRJ) < 11 or PSRJ[-1] == 'A':
                    inpc = float(pulsar_line.split()[1])
                    otpc = float(pulsar_line.split()[2])
                    if not all_check:
                        psrbeg = int(obsbeg+obsdur*inpc)
                        psrend = int(obsbeg+obsdur*otpc)
                    temp = fpio.get_psrcat_ra_dec(pulsar_list=[PSRJ])
                    temp = fpio.format_ra_dec(temp, ra_col = 1, dec_col = 2)
                    jname, raj, decj = temp[0]
                    print PSRJ, raj, decj

                    fits_dir = '/group/mwaops/vcs/{0}/pointings/{1}_{2}/'.format(obsid, raj, decj)
                    blind_pipe.process_vcs_wrapper(obsid, psrbeg, psrend, [raj,decj], args, DI_dir,\
                                        fits_dir, 0, False, None, None, pulsar_check=jname)

    os.remove(known_pulsar_file)
    return



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    Beamforms, folds on all known pulsars for an observation. If a pulsar is detected it's uploaded to the MWA Pulsar Database.

    Based on a script written by Mengyao Xue.
    """)
    parser.add_argument('-o','--obsid',type=str,help='The observation ID of the fits file to be searched')
    group_beamform = parser.add_argument_group('group_beamform','Beamforming Options')
    group_beamform.add_argument("--DI_dir", default=None, help="Directory containing either Direction Independent Jones Matrices (as created by the RTS) or calibration_solution.bin as created by Andre Offringa's tools.[no default]")
    group_beamform.add_argument('--cal_obs', '-O', type=int, help="Observation ID of calibrator you want to process.", default=None)
    group_beamform.add_argument("--pulsar_file", default=None, help="Location of a file containting the pointings to be processed. Made using grid.py.")
    group_beamform.add_argument("-a", "--all", action="store_true",  help="Uses all of the combined data available. If the options isn't used it'll only use the start and stops times that are recommened by find_pulsar_in_obs.py.")
    args=parser.parse_args()
    
    #option parsing
    if not args.obsid:
        print "Please input observation id by setting -o or --obsid. Exiting"
        quit()
    
    if not args.cal_obs:
        print "Please input calibration observation id by setting -O or --cal_obs. Exiting"
        quit()

    if not args.DI_dir:
        args.DI_dir = "/group/mwaops/vcs/{0}/cal/{1}/rts/".format(args.obsid, args.cal_obs)
        print "No DI_dir given so assuming {0} is the directory".format(args.DI_dir)
    
    beamform_and_fold(args.obsid, args.DI_dir, args.all, args)

