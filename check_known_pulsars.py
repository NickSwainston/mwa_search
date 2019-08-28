#! /usr/bin/env python3

import argparse
import find_pulsar_in_obs as fpio
import os
import glob
import subprocess
import numpy as np

from astropy.coordinates import SkyCoord
import astropy.units as u


import mwa_search_pipeline as search_pipe
from mwa_metadb_utils import get_common_obs_metadata as get_meta
from mwa_metadb_utils import obs_max_min
import config
from grid import get_grid

def find_beg_end(obsid, base_path="/group/mwaops/vcs/"):

    #looks through the comined files of the input obsid and returns the max and min in gps time
    #TODO have some sort of check to look for gaps
    if glob.glob("{0}/{1}/combined/{1}*_ics.dat".format(base_path, obsid)):
        combined_files = glob.glob("{0}/{1}/combined/{1}*_ics.dat".format(base_path, obsid))
    else:
        meta_data = get_meta(obsid)
        channels = meta_data[-1]
        combined_files = glob.glob("{0}/{1}/combined/{1}*_ch{2}.dat".\
                                   format(base_path, obsid, channels[-1]))
    comb_times = []
    for comb in combined_files:
        comb_times.append(int(comb.split("_")[1]))
    beg = min(comb_times)
    end = max(comb_times)

    return beg, end

def beamform_and_fold(obsid, DI_dir, cal_obs, args, psrbeg, psrend,
                      vdif_check=False, product_dir='/group/mwaops/vcs',
                      mwa_search_version='master'):


    #obsbeg, obsend, obsdur = file_maxmin.print_minmax(obsid)

    #wrapping for find_pulsar_in_obs.py
    names_ra_dec = np.array(fpio.grab_source_alog(max_dm=250))
    obs_data, meta_data = fpio.find_sources_in_obs([obsid], names_ra_dec, dt_input=100)
    channels = meta_data[-1]

    pointing_list = []
    jname_list = []
    for pulsar_line in obs_data[obsid]:
        PSRJ = pulsar_line[0]
        if len(PSRJ) < 11 or PSRJ[-1] == 'A' or PSRJ[-2:] == 'aa':
            temp = fpio.get_psrcat_ra_dec(pulsar_list=[PSRJ])
            temp = fpio.format_ra_dec(temp, ra_col = 1, dec_col = 2)
            jname, raj, decj = temp[0]
            #get pulsar period
            cmd = ['psrcat', '-c', 'p0', jname]
            output = subprocess.Popen(cmd,stdout=subprocess.PIPE).communicate()[0].decode()
            period = output.split('\n')[4].split()[1] #in s

            if '*' in period:
                print("WARNING: Period not found in ephermeris for {0}".format(jname))
                period=0
            else:
                period = float(period)*1000.
            print("{0:12} RA: {1} Dec: {2} Period: {3:8.2f} (ms) Begin {4} End {5}".format(
                   PSRJ, raj, decj, period, psrbeg, psrend))

            jname_temp_list = []
            if PSRJ[-1] == 'A' or PSRJ[-2:] == 'aa':
                #Got to find all the pulsar J names with other letters
                vdif_check = True
                for pulsar_l in pulsar_lines:
                    if pulsar_l.startswith(PSRJ[:-2]):
                        jname_temp_list.append(pulsar_l.split()[0])
            else:
                jname_temp_list.append(jname)
                #if b'*' in period:
                #    continue
                #if float(period) < .05 :
                #    vdif_check = True

            #convert to radians
            coord = SkyCoord(raj, decj, unit=(u.hourangle,u.deg))
            rar = coord.ra.radian #in radians
            decr = coord.dec.radian

            #make a grid around each pulsar
            grid_sep = np.radians(0.3 * 0.6) #TODO work this out for each obs
            rads, decds = get_grid(rar, decr, grid_sep, 2)

            #convert back to sexidecimals
            coord = SkyCoord(rads,decds,unit=(u.deg,u.deg))
            rajs = coord.ra.to_string(unit=u.hour, sep=':')
            decjs = coord.dec.to_string(unit=u.degree, sep=':')
            temp = []
            for raj, decj in zip(rajs, decjs):
                temp.append([raj, decj])
            pointing_list_list = fpio.format_ra_dec(temp, ra_col = 0, dec_col = 1)
            for prd in pointing_list_list:
                jname_list.append(jname_temp_list)
                pointing_list.append("{0} {1}".format(prd[0], prd[1]))
    #Setting vdif to false since multi-pixel doesn't have vdif working yet
    vdif_check = False
    relaunch_script = "mwa_search_pipeline.py -o {0} -O {1} --DI_dir {2} -b {3} -e {4} --channels".format(obsid, cal_obs, DI_dir, psrbeg, psrend)
    for ch in channels:
        relaunch_script = "{0} {1}".format(relaunch_script, ch)
    search_opts = search_pipe.search_options_class(obsid, cal_id=cal_obs,
                              begin=psrbeg, end=psrend, channels=channels,
                              args=args, DI_dir=DI_dir, relaunch_script=relaunch_script,
                              search_ver=mwa_search_version)
    search_pipe.beamform(search_opts, pointing_list, pulsar_list_list=jname_list)



    return



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    Beamforms, folds on all known pulsars for an observation. If a pulsar is detected it's uploaded to the MWA Pulsar Database.

    Based on a script written by Mengyao Xue.
    """)
    parser.add_argument('-o','--obsid',type=str,help='The observation ID of the fits file to be searched')
    parser.add_argument("--DI_dir", default=None, help="Directory containing either Direction Independent Jones Matrices (as created by the RTS) or calibration_solution.bin as created by Andre Offringa's tools.[no default]")
    parser.add_argument('-O','--cal_obs', type=int, help="Observation ID of calibrator you want to process.", default=None)
    parser.add_argument("-v", "--vdif", action="store_true",  help="Create vdif files for all pulsars. Default is to only create vdif files for pulsar with a period shorter than 50 ms.")
    parser.add_argument("-b", "--begin", type=int, help="First GPS time to process [no default]")
    parser.add_argument("-e", "--end", type=int, help="Last GPS time to process [no default]")
    parser.add_argument("-a", "--all", action="store_true", default=False, help="Perform on entire observation span. Use instead of -b & -e")
    parser.add_argument('--mwa_search_version', type=str, default='master',
                    help="The module version of mwa_search to use. Default: master")
    args=parser.parse_args()

    #option parsing
    if not args.obsid:
        print("Please input observation id by setting -o or --obsid. Exiting")
        quit()

    if not args.cal_obs:
        print("Please input calibration observation id by setting -O or --cal_obs. Exiting")
        quit()

    comp_config = config.load_config_file()
    if not args.DI_dir:
        args.DI_dir = "{0}/{1}/cal/{2}/rts/".format(comp_config['base_product_dir'],
                                                    args.obsid, args.cal_obs)
        print("No DI_dir given so assuming {0} is the directory".format(args.DI_dir))

    if args.begin and args.end:
        beg = args.begin
        end = args.end
    elif args.all:
        beg, end = obs_max_min(args.obsid)
    else:
        find_beg_end(args.obsid, base_path=comp_config['base_product_dir'])


    beamform_and_fold(args.obsid, args.DI_dir, args.cal_obs, args, beg, end,
                      vdif_check=args.vdif, product_dir=comp_config['base_product_dir'],
                      mwa_search_version=args.mwa_search_version)



