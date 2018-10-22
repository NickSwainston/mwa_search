#! /usr/bin/env python

import subprocess
import os
import argparse
import urllib
import urllib2
import json
import time
import mwa_metadb_utils as meta
import glob
import csv

parser = argparse.ArgumentParser(description="""
Wraps the splice_psrfits.sh script to automate it. Should be run from the foulder containing the files.
""")
parser.add_argument('-o','--observation',type=str,help='The observation ID to be used.')
parser.add_argument('-i','--incoh',action='store_true',help='Use this check if there are and incoh files from the beamformer.')
parser.add_argument('-d','--delete',action='store_true',help='This will cause the script to remove the unspliced files if splice_psrfits.sh succeeds (error code of 0).')
parser.add_argument('-w','--work_dir',type=str,help='Working directory of the vcs files.', default="./")
args=parser.parse_args()

    
obsid = args.observation
chan_st = range(24) 

with open('/group/mwaops/nswainston/code/blindsearch_scripts/obs_meta.csv', 'rb') as csvfile:
    spamreader = csv.reader(csvfile)
    next(spamreader, None)
    obsid_meta_file = []
    for row in spamreader:
        obsid_meta_file.append(row)

obs_foun_check = False
for omi in range(len(obsid_meta_file)):
  if int(obsid) == int(obsid_meta_file[omi][0]):
      print "getting obs metadata from obs_meta.csv"
      ob, ra, dec, duration, delays,centrefreq, channels =\
                  obsid_meta_file[omi]
      channels = map(int,channels[1:-1].split(","))
      obs_foun_check = True
if not obs_foun_check:
  ob, ra, dec, duration, delays,centrefreq, channels =\
          meta.get_common_obs_metadata(obsid)  
  
  with open('obs_meta.csv', 'a') as csvfile:
      spamwriter = csv.writer(csvfile)
      spamwriter.writerow([ob, ra, dec, duration, delays,centrefreq, channels])

print "Obtaining metadata from http://mwa-metadata01.pawsey.org.au/metadata/ for OBS ID: " + str(obsid)
beam_meta_data = meta.getmeta(service='obs', params={'obs_id':obsid})
channels = beam_meta_data[u'rfstreams'][u"0"][u'frequencies']

chan_order=[]
for i in range(24):
    chan_order.append(channels[i])

    #if channels[i] < 129:
    #    chan_order.append(channels[i])
    #else:
    #    chan_order.append( channels[chan_st[-(channels[i] -128)]])
print "Chan order: {}".format(chan_order)

#move into working fdir
old_dir = os.getcwd()
os.chdir(args.work_dir)

if args.incoh:
    submit_line = 'ls *'+str(args.observation)+'_ch'+str(channels[-1])+'_incoh* | wc -l'
else:
    submit_line = 'ls *'+str(args.observation)+'_ch'+str(channels[-1])+'* | wc -l'
print submit_line
submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
out_lines = submit_cmd.stdout
for l in out_lines:
    n_fits = l
    
if int(n_fits) == 0:
    new_bf = False
    chan_order=[]
    for i in range(24):
        if channels[i] < 129:
            chan_order.append(i+1)
        else:
            chan_order.append( chan_st[-(channels[i] -128)] + 1)
    print chan_order

    submit_line = 'ls *'+str(args.observation)+'_01* | wc -l'
    submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
    out_lines = submit_cmd.stdout
    for l in out_lines:
        n_fits = l
else:
    new_bf = True

for n in range(int(n_fits)):
    submit_line = submit_line_incoh = 'splice_psrfits '
    for ch in chan_order:
        submit_line += '*{}_ch{:03d}_{:04d}.fits '.format(obsid, ch, int(n)+1)
        submit_line_incoh +=  '*{}_ch{:03d}_incoh_{:04d}.fits '.format(obsid, ch, int(n)+1)
    submit_line += 'temp_'+str(int(n)+1)
    submit_line_incoh += 'temp_incoh_'+str(int(n)+1)
    
    if args.incoh:
        print submit_line_incoh
        submit_cmd = subprocess.Popen(submit_line_incoh,shell=True,stdout=subprocess.PIPE)
        out_lines = submit_cmd.stdout        
        for l in out_lines:
            print l,
        time.sleep(5)
        print 'temp_incoh_'+str(int(n)+1)+'_0001.fits', '{}_incoh_{:04d}.fits'.format(obsid, int(n)+1)
        os.rename('temp_incoh_'+str(int(n)+1)+'_0001.fits',
                  '{}_incoh_{:04d}.fits'.format(obsid, int(n)+1))
        
        #wait to get error code
        (output, err) = submit_cmd.communicate()  
        p_status = submit_cmd.wait()
        print "exist code: " + str(submit_cmd.returncode)
        if args.delete and int(submit_cmd.returncode) == 0:
            for fd in submit_line_incoh[15:].split(" ")[:-1]:
                fd = glob.glob(fd)[0]
                print "Deleting: " + str(fd)
                if os.path.isfile(fd):
                     os.remove(fd)
    print submit_line
    submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
    out_lines = submit_cmd.stdout
    for l in out_lines:
        print l,
    time.sleep(5)
    print 'temp_'+str(int(n)+1)+'_0001.fits', '{}_{:04d}.fits'.format(obsid, int(n)+1)
    os.rename('temp_'+str(int(n)+1)+'_0001.fits',
              '{}_{:04d}.fits'.format(obsid, int(n)+1))
    
    #wait to get error code
    (output, err) = submit_cmd.communicate()  
    p_status = submit_cmd.wait()
    
    print "exist code: " + str(submit_cmd.returncode)
    if args.delete and int(submit_cmd.returncode) == 0:
        for fd in submit_line[15:].split(" ")[:-1]:
            fd = glob.glob(fd)[0]
            print "Deleting: " + str(fd)
            if os.path.isfile(fd):
                 os.remove(fd)
