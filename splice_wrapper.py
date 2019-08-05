#! /usr/bin/env python3

import subprocess
import os
import argparse
import time
import mwa_metadb_utils as meta
import glob
from shutil import copy
import socket

parser = argparse.ArgumentParser(description="""
Wraps the splice_psrfits.sh script to automate it. Should be run from the foulder containing the files.
""")
parser.add_argument('-o','--observation',type=str,help='The observation ID to be used.')
parser.add_argument('-i','--incoh',action='store_true',help='Use this check if there are and incoh files from the beamformer.')
parser.add_argument('-d','--delete',action='store_true',help='This will cause the script to remove the unspliced files if splice_psrfits.sh succeeds (error code of 0).')
parser.add_argument('-w','--work_dir',type=str,help='Working directory of the vcs files.', default="./")
parser.add_argument('-c', '--channels', type=int, nargs=24, help='A list of the observations channel IDs for example "-c 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132". If this option is not used a metadata call will be used to find the channel IDs.')
args=parser.parse_args()


obsid = args.observation

# Check if already spliced
if glob.glob('{0}*fits'.format(args.observation)) and \
   not glob.glob('*_{0}*fits'.format(args.observation)):
    print('All files are already spliced so exiting')
    exit()

# Get frequency channels
if args.channels:
    channels = args.channels
else:
    print("Obtaining metadata from http://mwa-metadata01.pawsey.org.au/metadata/ for OBS ID: " + str(obsid))
    beam_meta_data = meta.getmeta(service='obs', params={'obs_id':obsid})
    channels = beam_meta_data[u'rfstreams'][u"0"][u'frequencies']
print("Chan order: {}".format(channels))

hostname = socket.gethostname()
if hostname.startswith('john') or hostname.startswith('bryan'):
    #If on ozstar use their SSD to improve I/O
    SSD_file_dir = '{}/'.format(os.environ['JOBFS'])
    print(SSD_file_dir)
else:
    SSD_file_dir = ''


# Move into working dir
old_dir = os.getcwd()
if hostname.startswith('john') or hostname.startswith('bryan'):
    os.chdir(SSD_file_dir)
else:
    os.chdir(args.work_dir)


# Getting number of files list
if args.incoh:
    n_fits_file = glob.glob('{}/*{}_incoh_ch{}_*fits'.format(args.work_dir, obsid, channels[-1]))
else:
    n_fits_file = glob.glob('{}/*{}*_ch{}*fits'.format(args.work_dir, obsid, channels[-1]))
n_fits = []
for file_name in n_fits_file:
    n_fits.append(int(file_name[-9:-5]))
n_fits.sort()
print('Fits number order: {}'.format(n_fits))

for n in n_fits:
    # List unspliced files
    unspliced_files = []
    for ch in channels:
        if args.incoh:
            unspliced_files.append(glob.glob('{}/*{}_incoh_ch{:03d}_{:04d}.fits'.format(args.work_dir,
                                                   obsid, ch, n))[0])
        else:
            unspliced_files.append(glob.glob('{}/*{}*_ch{:03d}_{:04d}.fits'.format(args.work_dir,
                                                  obsid, ch, n))[0])

    
    if hostname.startswith('john') or hostname.startswith('bryan'):
        print("Moving each channel onto $JOBFS")
        for us_file in unspliced_files:
            copy(us_file, SSD_file_dir)

    # Create splice command and submit 
    submit_line = 'splice_psrfits '
    for us_file in unspliced_files:
        if hostname.startswith('john') or hostname.startswith('bryan'):
            file_on_SSD = "{}{}".format(SSD_file_dir, us_file.replace(args.work_dir, ''))
            submit_line += '{} '.format(file_on_SSD)
        else:
            submit_line += '{} '.format(us_file)
    submit_line += '{}temp_{}'.format(SSD_file_dir, n)

    print(submit_line)
    submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
    out_lines = submit_cmd.stdout
    for l in out_lines:
        print(l.decode()[:-1])
    time.sleep(1)


    print('Finished {}_{:04d}.fits'.format(obsid, n))
    if hostname.startswith('john') or hostname.startswith('bryan'):
        #Use copy because it's faster than mv
        copy('{}temp_{}_0001.fits'.format(SSD_file_dir, n),
             '{}/{}_{:04d}.fits'.format(args.work_dir, obsid, n))
        ssd_files = glob.glob('{}/*'.format(SSD_file_dir))
        for sf in ssd_files:
            os.remove(sf)
    else:
        os.rename('temp_{}_0001.fits'.format(n), '{}_{:04d}.fits'.format(obsid, n))

    #wait to get error code
    (output, err) = submit_cmd.communicate()
    p_status = submit_cmd.wait()

    print("exit code: " + str(submit_cmd.returncode))
    if args.delete and int(submit_cmd.returncode) == 0:
        for us_file in unspliced_files:
            print("Deleting: " + str(us_file))
            if os.path.isfile(us_file):
                 os.remove(us_file)
