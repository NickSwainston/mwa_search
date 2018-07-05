#!/usr/bin/env python

from scipy.interpolate import griddata
import matplotlib
import matplotlib.pyplot as plt
import sys
import subprocess

# numerical and maths modules
import numpy as np
from astropy.coordinates import SkyCoord,EarthLocation,AltAz,angles
from astropy.time import Time
from astropy.io import fits
import astropy.units as u
from astropy.constants import c,k_B
from math import cos,sin,acos,asin

#utility and processing modules
import os,sys
#from mpi4py import MPI
import argparse
from mwapy import ephem_utils,metadata
from mwapy.pb import primary_beam as pb
from mwapy.pb.mwa_tile import h2e
import urllib
import urllib2
import json


def getmeta(service='obs', params=None):
    """
    Given a JSON web service ('obs', find, or 'con') and a set of parameters as
    a Python dictionary, return the RA and Dec in degrees from the Python dictionary.
    
    getmeta(service='obs', params=None)
    """
    BASEURL = 'http://mwa-metadata01.pawsey.org.au/metadata/'
    if params:
        data = urllib.urlencode(params)  # Turn the dictionary into a string with encoded 'name=value' pairs
    else:
        data = ''
    #Validate the service name
    if service.strip().lower() in ['obs', 'find', 'con']:
        service = service.strip().lower()
    else:
        print "invalid service name: %s" % service
        return
    #Get the data
    try:
        result = json.load(urllib2.urlopen(BASEURL + service + '?' + data))
    except urllib2.HTTPError as error:
        print "HTTP error from server: code=%d, response:\n %s" % (error.code, error.read())
        return
    except urllib2.URLError as error:
        print "URL or network error: %s" % error.reason
        return
    #Return the result dictionary
    return result


def get_obstime_duration(obsid,fdir="."):
    """
    Funciton to grab the recorded start-time and duration of the observation
    
    Input:
      obsid - the GPS observation ID

    Return:
      a list containing the folloing two items (in order):
        list[0] = observation starting time in UTC
        list[1] = observation duration in seconds
    """
    # metafits file will already have been downloaded
    f = fits.open('{0}_metafits_ppds.fits'.format(obsid))
    
    return [f[0].header['DATE-OBS'],f[0].header['EXPOSURE']]
    

def FWHM(X,Y):
    half_max = max(Y) / 2.
    #find when function crosses line half_max (when sign of diff flips)
    #take the 'derivative' of signum(half_max - Y[])
    d = np.sign(half_max - np.array(Y[0:-1])) - np.sign(half_max - np.array(Y[1:]))
    #plt.plot(X,d) #if you are interested
    #find the left and right most indexes
    for i in range(len(d)):
        if d[i] > 0.0:
            left_idx = i
        if d[i] < 0.0:
            right_idx = i 
    return X[right_idx] - X[left_idx] #return the difference (full width)
    
 
def getTargetAZZA(ra,dec,time,lat=-26.7033,lon=116.671,height=377.827):
    """
    Function to get the target position in alt/az at a given EarthLocation and Time.
    
    Default lat,lon,height is the centre of  MWA.

    Input:
      ra - target right ascension in astropy-readable format
      dec - target declination in astropy-readable format
      time - time of observation in UTC (i.e. a string on form: yyyy-mm-dd hh:mm:ss.ssss)
      lat - observatory latitude in degrees
      lon - observatory longitude in degrees

    Returns:
      a list containing four elements in the following order:
        list[0] = target azimuth in radians
        list[1] = target zenith angle in radians
        list[2] = target azimuth in degrees
        list[3] = target zenith angle in degrees
    """
    #print "Creating EarthLocation for: lat = {0} deg, lon = {1} deg".format(lat,lon)
    location = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=height*u.m)
    
    #print "Creating SkyCoord for target at (RA,DEC) = ({0},{1})".format(ra,dec)
    coord = SkyCoord(ra,dec,unit=(u.hourangle,u.deg))
    #print "Target at: ({0},{1}) deg".format(coord.ra.deg,coord.dec.deg)
    
    obstime = Time(time)
    #print "Observation time: {0}".format(obstime.iso)
    
    #print "Converting to Alt,Az for given time and observatory location..."
    altaz = coord.transform_to(AltAz(obstime=obstime,location=location))
    #print "Target (Alt,Az) = ({0},{1}) deg".format(altaz.alt.deg,altaz.az.deg)
    
    #print "Converting to (Az,ZA)"
    az = altaz.az.rad 
    azdeg = altaz.az.deg
     
    za = np.pi/2 - altaz.alt.rad
    zadeg = 90 - altaz.alt.deg
    
    #print "Target (Az,ZA) = ({0},{1}) deg".format(azdeg,zadeg)

    return [az,za,azdeg,zadeg]
    
    
def getTargetradec(az,za,time,lst,lat=-26.7033,lon=116.671,height=377.827):
    """
    Function to get the target position in ra dec at a given EarthLocation and Time.
    
    Default lat,lon,height is the centre of  MWA.

    Input:
      az - target aximuth in radians
      za - target zenith in radians
      time - time of observation in UTC (i.e. a string on form: yyyy-mm-dd hh:mm:ss.ssss)
      lat - observatory latitude in degrees
      lon - observatory longitude in degrees

    Returns:
      a list containing four elements in the following order:
        list[0] = target ra in degrees
        list[1] = target dec in degrees
    """
    
    ha,dec = h2e(az,za,lat) #hour angle and dec in degrees
    ra = lst-ha
    

    return [ra,dec]
    
    
def get_freqs(obsid):
    print "Obtaining metadata from http://mwa-metadata01.pawsey.org.au/metadata/ for OBS ID: " + str(args.obsid)
    beam_meta_data = getmeta(service='obs', params={'obs_id':obsid})
    channels = beam_meta_data[u'rfstreams'][u"0"][u'frequencies']
    minfreq = float(min(channels))
    maxfreq = float(max(channels))
    centrefreq = 1.28 * (minfreq + (maxfreq-minfreq)/2) #in MHz
    return [channels, minfreq, maxfreq, centrefreq]
    

def two_floats(value):
    values = value.split()
    if len(values) != 2:
        raise argparse.ArgumentError
    return values  
    
    
#gird movements all in rad
def left(ra_in, dec_in, fwhm):
    dec_out = dec_in 
    ra_out = ra_in - acos( (cos(fwhm) - cos(dec_in)**2) / (sin(dec_in)**2 ) )
    return [ra_out,dec_out]
    
def right(ra_in, dec_in, fwhm):
    dec_out = dec_in 
    ra_out = ra_in + acos( (cos(fwhm) - cos(dec_in)**2) / (sin(dec_in)**2 ) )
    return [ra_out,dec_out]
    
    
def up_left(ra_in, dec_in, fwhm):
    #fwhm will be slighty off due to it being measured from the starting point
    half_fwhm_approx = acos( (cos(fwhm) - cos(dec_in)**2) / (sin(dec_in)**2 ) )
    dec_out = dec_in + acos( cos(half_fwhm_approx)/cos(fwhm/2.) )
    ra_out = ra_in - acos( (cos(fwhm/2.) - cos(dec_out)**2) / (sin(dec_out)**2 ) )
    return [ra_out,dec_out]
    
def up_right(ra_in, dec_in, fwhm):
    half_fwhm_approx = acos( (cos(fwhm) - cos(dec_in)**2) / (sin(dec_in)**2 ) )
    dec_out = dec_in + acos( cos(half_fwhm_approx)/cos(fwhm/2.) ) 
    ra_out = ra_in + acos( (cos(fwhm/2.) - cos(dec_out)**2) / (sin(dec_out)**2 ) )
    return [ra_out,dec_out]
    
def down_left(ra_in, dec_in, fwhm):
    half_fwhm_approx = acos( (cos(fwhm) - cos(dec_in)**2) / (sin(dec_in)**2 ) )
    dec_out = dec_in - acos( cos(half_fwhm_approx)/cos(fwhm/2.) )
    ra_out = ra_in - acos( (cos(fwhm/2.) - cos(dec_out)**2) / (sin(dec_out)**2 ) )
    return [ra_out,dec_out]
    
def down_right(ra_in, dec_in, fwhm):
    half_fwhm_approx = acos( (cos(fwhm) - cos(dec_in)**2) / (sin(dec_in)**2 ) )
    dec_out = dec_in - acos( cos(half_fwhm_approx)/cos(fwhm/2.) )  
    ra_out = ra_in + acos( (cos(fwhm/2.) - cos(dec_out)**2) / (sin(dec_out)**2 ) )
    return [ra_out,dec_out]
    
    
def calc_fwhm(obsid, pointing, opts_string):
    channels, minfreq, maxfreq, centrefreq = get_freqs(args.obsid)
    
    #create a batch script to one pabeam along a strip to calculate fwhm
    with open('Beam_calc_' + str(args.obsid) + '.batch','w') as batch_file:
        if args.gpuq:
            batch_line = "#!/bin/bash -l\n" +\
                        "#SBATCH --partition=gpuq\n"
            batch_file.write(batch_line)
        else:
            batch_line = "#!/bin/bash -l\n" +\
                        "#SBATCH --partition=workq\n"
            batch_file.write(batch_line)
        batch_line ="#SBATCH --time=1:00:00\n" +\
                    "#SBATCH --output=Beam_calc_" + str(obsid) + ".out\n" +\
                    "#SBATCH --nodes=1\n" +\
                    "#SBATCH --account=mwaops\n" +\
                    "nnodes=1 # number of requested nodes - make sure this matches the --nodes SBATCH variable\n" +\
                    'nprocesses=`bc -l <<< "${nnodes}*16"` # assuming 16 parallel processes available per node\n' +\
                    "obsid="+str(obsid) +" # observation ID\n" +\
                    "freq="+str(centrefreq)+" # centre observing frequency in Hz\n" +\
                    "eff=0.98037 # antenna efficiency\n" +\
                    'ra='+str(pointing[0])+' # target RA\n' +\
                    "dec="+str(pointing[1])+" # target DEC\n" +\
                    'flags="23 71" # flagged tiles (as in RTS solutions)\n' +\
                    "tres=0.01 # resolution for ZA grid\n" +\
                    "pres=0.01 # resolution for Azimuth grid\n" +\
                    "obstime="+str(obsid) +" # GPS time to evaluate the beam\n" +\
                    'odir="/group/mwaops/nswainston/pabeam/output/" # where to dump the outputs\n' +\
                    'radius="1"\n' +\
                    'echo "aprun -n $nprocesses python /group/mwaops/nswainston/pabeam/pabeam.py -o ${obsid} -f ${freq} -t ${obstime} -e ${eff} -p ${ra} ${dec} --flagged_tiles ${flags} --grid_res ${tres} ${pres} --out_dir ${odir} -r ${radius} --write --strip" \n' +\
                    'time aprun -n $nprocesses python /group/mwaops/nswainston/pabeam/pabeam.py --obsid ${obsid} -f ${freq} -t ${obstime} -e ${eff} -p "${ra} ${dec}" --flagged_tiles ${flags} --grid_res ${tres} ${pres} --out_dir ${odir} -r ${radius} --out_name ' + str(obsid) + "_fwhm --write --strip\n" +\
                    "python /group/mwaops/nswainston/pabeam/grid.py -m g " + opts_string
        
        batch_file.write(batch_line)   
    submit_line = 'sbatch Beam_calc_' + str(args.obsid) + '.batch'
    submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
    print submit_cmd.communicate()[0],
    print "Code will continue when job completes. Should take around 15 mins assuming no queue time."


def hex_grid(ra0,dec0,centre_fwhm, dec_lim_arg):
    #start location list [loop number][shape corner (6 for hexagon 4 for square)][number from corner]
    #each item has [ra,dec,fwhm] in radians
    pointing_list = [[[[ra0,dec0,centre_fwhm]]]]
    orig_centre_fwhm = centre_fwhm
    print "Calculating the tile positions"

    if dec_lim_arg:
        dec_lim_matrix = np.empty([len(del_lim_arg)/3,3])
        for i in range(len(dec_lim_arg)):
            dec_lim_matrix[i/3,i%3] = dec_lim_arg[i]


    for l in range(1,args.loop+1):
        #different step for each corner
        loop_temp = []
        if l%50 == 0:
            print "Calculating loop number: " +str(l)
        for c in range(6):
            corner_temp = []
            if l == 1:
                if c == 0:
                    ra,dec =left(ra0,dec0,centre_fwhm)
                if c == 1:
                    ra,dec =up_left(ra0,dec0,centre_fwhm)
                if c == 2:
                    ra,dec =up_right(ra0,dec0,centre_fwhm)
                if c == 3:
                    ra,dec =right(ra0,dec0,centre_fwhm)
                if c == 4:
                    ra,dec =down_right(ra0,dec0,centre_fwhm)
                if c == 5:
                    ra,dec =down_left(ra0,dec0,centre_fwhm)
                corner_temp.append([ra,dec])
            else:
                for n in range(l):
                    if l == (n + 1):
                        if dec_lim_arg:
                            centre_fwhm = centre_fwhm_orig
                            for d in dec_lim_matrix:
                                if pointing_list[l-1][c+1][0][1] > dec_lim_matrix[d][1] and\
                                        pointing_list[l-1][c+1][0][1] < dec_lim_matrix[d][2]:
                                    centre_fwhm = dec_lim_matrix[d][0]
                        #change the 2 for each loop
                        #uses next corner
                        #TODO use updated fwhm
                        if c == 0:
                            ra,dec =left(pointing_list[l-1][c+1][0][0],
                                         pointing_list[l-1][c+1][0][1],centre_fwhm)
                        if c == 1:
                            ra,dec =up_left(pointing_list[l-1][c+1][0][0],
                                         pointing_list[l-1][c+1][0][1],centre_fwhm)
                        if c == 2:
                            ra,dec =up_right(pointing_list[l-1][c+1][0][0],
                                         pointing_list[l-1][c+1][0][1],centre_fwhm)
                        if c == 3:
                            ra,dec =right(pointing_list[l-1][c+1][0][0],
                                         pointing_list[l-1][c+1][0][1],centre_fwhm)
                        if c == 4:
                            ra,dec =down_right(pointing_list[l-1][c+1][0][0],
                                         pointing_list[l-1][c+1][0][1],centre_fwhm)
                        if c == 5:
                            ra,dec =down_left(pointing_list[l-1][0][0][0],
                                         pointing_list[l-1][0][0][1],centre_fwhm)
                        
                    else:
                        if dec_lim_arg:
                            centre_fwhm = centre_fwhm_orig
                            for d in dec_lim_matrix:
                                if pointing_list[l-1][c][n][1] > dec_lim_matrix[d][1] and\
                                        pointing_list[l-1][c][n][1] < dec_lim_matrix[d][2]:
                                    centre_fwhm = dec_lim_matrix[d][0]

                        if c == 0:
                            ra,dec =left(pointing_list[l-1][c][n][0],
                                         pointing_list[l-1][c][n][1],centre_fwhm)
                        if c == 1:
                            ra,dec =up_left(pointing_list[l-1][c][n][0],
                                            pointing_list[l-1][c][n][1],centre_fwhm)
                        if c == 2:
                            ra,dec =up_right(pointing_list[l-1][c][n][0],
                                             pointing_list[l-1][c][n][1],centre_fwhm)
                        if c == 3:
                            ra,dec =right(pointing_list[l-1][c][n][0],
                                          pointing_list[l-1][c][n][1],centre_fwhm)
                        if c == 4:
                            ra,dec =down_right(pointing_list[l-1][c][n][0],
                                               pointing_list[l-1][c][n][1],centre_fwhm)
                        if c == 5:  
                            ra,dec =down_left(pointing_list[l-1][c][n][0],
                                              pointing_list[l-1][c][n][1],centre_fwhm)
                    corner_temp.append([ra,dec])
            loop_temp.append(corner_temp)
        pointing_list.append(loop_temp)
    return pointing_list
    
parser = argparse.ArgumentParser(description="""
Makes a hexogonal grid pattern around a pointing for a MWA VCS observation.
python grid.py -p 00:24:30.00 \"-72:04:30.00\" -o 1163853320 \n
python /group/mwaops/nswainston/pabeam/grid.py -o 1166459712 -p "06:25:31.20 -36:40:48.0" -f 0.8
""")
parser.add_argument('-o', '--obsid',type=str,help='Observation ID')
parser.add_argument('-p', '--pointing',type=two_floats,help='Centre pointing in hh:mm:ss.sss dd\"mm\'ss.ssss')
parser.add_argument('--plot',action="store_true",help='Plots the output')
parser.add_argument('--gpuq',action="store_true",help='Uses the gpuq (used for when the workq is full).')
parser.add_argument('-f', '--fraction',type=float,help='Fraction of the full width half maximum to use as the distance between beam centres',default=0.85)
parser.add_argument('-d', '--deg_fwhm',type=float,help='Sets the FWHM in degrees. The script will not calculate the FWHM',default=0.3098)
parser.add_argument('--dec_range_fwhm',type=float,nargs='+',help='A list of FWHM and ranges in the order of: "FWHM1 decmin1 decmax1 FWHM2 decmin2 decmax2"')
parser.add_argument('-t', '--type',type=str,help='Can be put in either "hex" or "square" tiling mode. Default is hex.',default='hex')
parser.add_argument('-m', '--mode',type=str,help='Program mode used internally by code (needed to start program again after finishing on the slurm queue). "f" is to calculate the fwhm and the default, "g" Calc grid positions and "p" will plot.',default='f')
parser.add_argument('-l', '--loop',type=int,help='Number  of "loops" around the centre pointing the code will calculate. Default is 1',default=1)
parser.add_argument('--dec_range',type=float,nargs='+',help='Dec limits: "decmin decmax". Default -90 90', default=[-90,90])
parser.add_argument('--ra_range',type=float,nargs='+',help='RA limits: "ramin ramax". Default 0 390', default=[0,360])
parser.add_argument('-v','--verbose_file',action="store_true",help='Creates a more verbose output file with more information than make_beam.c can handle.')

args=parser.parse_args()

opts_string = ""
for k in args.__dict__:
    if args.__dict__[k] is not None:
        if k == "pointing":
            opts_string = opts_string + ' --' + str(k) + ' "' + str(args.__dict__[k][0]) +\
                     ' ' + str(args.__dict__[k][1]) + '"'
        else:
            opts_string = opts_string + ' --' + str(k) + ' ' + str(args.__dict__[k])
        
if not args.plot:
    channels, minfreq, maxfreq, centrefreq = get_freqs(args.obsid)
    
#get fwhm in radians
if args.mode == 'f':
    if args.deg_fwhm:
        centre_fwhm = np.radians(args.deg_fwhm)
    else:
        calc_fwhm(args.obsid, args.pointing,opts_string)

#calc grid positions
if args.mode == 'c' or args.deg_fwhm:
    if not args.deg_fwhm:
        oname = "/group/mwaops/nswainston/pabeam/output/"+str(args.obsid)+"_fwhm.0.dat"
        #oname = "1163853320_1163853320.0_0.00018496MHz_00:24:30.00_-72:04:30.00.0.dat"
        print "loading data"
        theta,phi,beam = np.genfromtxt(oname,comments=("#","*"),skip_header=14,usecols=(0,1,8),unpack=True)
        print "done"

        centre_fwhm = np.radians(FWHM(phi,beam))#*args.fraction #in radians
        print 'FWHM: ',centre_fwhm
    
    coord = SkyCoord(args.pointing[0],args.pointing[1],unit=(u.hourangle,u.deg))
    ra = coord.ra.radian #in radians
    dec = coord.dec.radian
    
    
    if args.type == 'hex':
        pointing_list = hex_grid(ra, dec, centre_fwhm*args.fraction, args.dec_range_fwhm)
    #TODO add square

    time = Time(float(args.obsid),format='gps')
    ra_decs = []      
    ras = []; decs = []; theta = []; phi = []; rads = []; decds = []
    
    print "Converting to ra dec"                
    for l in range(len(pointing_list)):
        if l%50 == 0:
            print "Calculating loop number: " +str(l)
        for c in range(len(pointing_list[l])):
            for n in range(len(pointing_list[l][c])):
                #format grid pointings
                rad = np.degrees(pointing_list[l][c][n][0])
                decd = np.degrees(pointing_list[l][c][n][1])
                
                if decd > 90.:
                    decd = decd - 180.
                coord = SkyCoord(rad,decd,unit=(u.deg,u.deg))
                rag = coord.ra.to_string(unit=u.hour, sep=':')
                decg = coord.dec.to_string(unit=u.degree, sep=':')
                #format the ra dec strings 
                if len(rag) > 11:
                    rag = rag[:11]
                if len(decg) > 12:
                    decg = decg[:12]
                    
                if len(rag) == 8:
                    rag = rag + '.00'
                if len(decg) == 9:
                    decg = decg + '.00'


                if args.verbose_file:
                    az,za,azd,zad = getTargetAZZA(rag,decg,time)
                else:
                    az,za,azd,zad = [0,0,0,0]
                
                ras.append(rag)
                decs.append(decg)
                theta.append(az)
                phi.append(za)
                rads.append(rad)
                decds.append(decd)
                #ra_decs.append([rag,decg,az,za,rad,decd])
    print "Recording the poisitons in grid_positions.txt"            
    with open('grid_positions.txt','w') as out_file:
        if args.verbose_file:
            out_line = "#ra   dec    az     za\n" 
            out_file.write(out_line)
        for i in range(len(rads)):
            if args.verbose_file:
                out_line = str(ras[i])+" "+str(decs[i])+" "+str(theta[i])+" "\
                            +str(phi[i])+" "+str(rads[i])+" "\
                            +str(decds[i])+"\n" 
            else:
                out_line = str(ras[i])+" "+str(decs[i])+"\n" 

            out_file.write(out_line) 
           
    #some ra and dec string for radec limited degreees
    radls = []
    decdls = []
    print "Recording the dec limited poisitons in grid_positions_dec_limited.txt"            
    with open('grid_positions_ra_dec_limited.txt','w') as out_file:
        if args.verbose_file:
            out_line = "#ra   dec    az     za\n" 
            out_file.write(out_line)
        for i in range(len(rads)):
            if  (args.dec_range[0] < float(decds[i]) < args.dec_range[1] ) and \
                (args.ra_range[0]  < float(rads[i]) < args.ra_range[1]):
                    if args.verbose_file:
                        out_line = str(ras[i])+" "+str(decs[i])+" "+str(theta[i])+" "\
                                    +str(phi[i])+" "+str(rads[i])+" "\
                                    +str(decds[i])+"\n" 
                    else:
                        out_line = str(ras[i])+" "+str(decs[i])+"\n" 
                    radls.append(rads[i])
                    decdls.append(decds[i])
                    out_file.write(out_line)
            
    #ras, decs, theta, phi, rads, decds
    #ras, decs, theta, phi, rads, decds = np.genfromtxt('/home/nswainst/blindsearch_scripts/output/grid_position.txt',unpack=True)    
    #ras, decs, theta, phi, rads, decds = np.genfromtxt('/group/mwaops/nswainston/code/blindsearch_scripts/grid_positions.txt',unpack=True)

    matplotlib.use('Agg')
    
    fig = plt.figure(figsize=(7, 7))
    
    plt.xlabel("ra (degrees)")
    plt.ylabel("dec (degrees)")
    ax = plt.gca()
    for i in range(len(ras)):
        if  -0.00001 < sin(np.radians(decds[i]))**2 < 0.00001:
            fwhm_circle = acos( (cos(centre_fwhm) - cos(np.radians(decds[i]))**2) / (0.00001+sin(np.radians(decds[i]))**2 ) ) /2.
        else:
            fwhm_circle = acos( (cos(centre_fwhm) - cos(np.radians(decds[i]))**2) / (sin(np.radians(decds[i]))**2 ) ) /2.
        #fwhm_circle = acos( cos(np.radians(decds[i]))**2 + sin(np.radians(decds[i]))**2 *cos(centre_fwhm) ) /2.
        circle = plt.Circle((rads[i],decds[i]),np.degrees(fwhm_circle),color='r',fill=False)
        ax.add_artist(circle)
    plt.axes().set_aspect('equal')
    plt.plot(rads,decds,'ko',ms=1)
    #plt.plot([114.25,97.5,115.5],[-30.65,-28.56,-28.36],'go',ms=5)
    #plt.plot([114.25,115.5],[-30.65,-28.36],'go',ms=5)
    #plt.axis([113.0, 116.7, -31.5, -27.5])
    plt.savefig('grid_positions_'+str(args.obsid)+'_ra_dec_f'+str(args.fraction)+'_l'+str(args.loop)+'.png',bbox_inches='tight', dpi =1000)
    
    """
    fwhm_circle = acos( (cos(centre_fwhm) - cos(np.radians(decds[i]))**2) / (sin(np.radians(decds[i]))**2 ) ) /2.
    plt.clf()
    
    plt.ylabel("phi (radians)")
    plt.xlabel("theta (radians)")
    ax = plt.gca(projection='polar')
    for i in range(len(theta)):
        circle = plt.Circle((theta[i],phi[i]),fwhm_circle, transform=ax.transData._b,color='r',alpha=0.4)#fill=False)
        ax.add_artist(circle)
    if args.loop < 6:
        plt.plot(theta,phi,'ko')
    else:
        plt.plot()
    plt.savefig('grid_positions_'+str(args.obsid)+'_alt_az_f'+str(args.fraction)+'_l'+str(args.loop)+'.png',bbox_inches='tight')
    
    
    #ras, decs, theta, phi, rads, decds = np.genfromtxt('/group/mwaops/nswainston/pabeam/output/grid_positions_dec_limited.txt',unpack=True) 
    """
    
    #ras, decs, theta, phi, rads, decds = np.genfromtxt('./grid_positions_ra_dec_limited.txt',unpack=True) 
    rads = radls
    decds = decdls

    #dec limited plot
    fig = plt.figure(figsize=(7, 7))
    plt.xlabel("ra (degrees)")
    plt.ylabel("dec (degrees)")
    ax = plt.gca()
    for i in range(len(rads)):
        fwhm_circle = acos( (cos(centre_fwhm) - cos(np.radians(decds[i]))**2) / (sin(np.radians(decds[i]))**2 ) ) /2.
        #fwhm_circle = acos( cos(np.radians(decds[i]))**2 + sin(np.radians(decds[i]))**2 *cos(centre_fwhm) ) /2.
        circle = plt.Circle((rads[i],decds[i]),np.degrees(fwhm_circle),color='r',fill=False)
        ax.add_artist(circle)
    plt.axes().set_aspect('equal')
    plt.plot(rads,decds,'ko',ms=1)
    #plt.plot([114.25,97.5,115.5],[-30.65,-28.56,-28.36],'go',ms=1)
    #plt.plot([114.25,115.5],[-30.65,-28.36],'go',ms=2)
    
    plt.savefig('grid_positions_'+str(args.obsid)+'_ra_dec_limited_f'+str(args.fraction)+'_l'+str(args.loop)+'.png',bbox_inches='tight', dpi =1000)
    
    
    
    
    print "Number of pointings: " + str(len(theta))
    
   
    
    
    
    if args.plot:
        job_id_list = []
        for i in range(len(ras)):
            with open('Beam_calc_' + str(args.obsid) + '_grid_' + str(i) + '.batch','w') as batch_file:
                batch_line = "#!/bin/bash -l\n" +\
                    "#SBATCH --partition=workq\n" +\
                    "#SBATCH --time=12:00:00\n" +\
                    "#SBATCH --output=Beam_calc_" + str(args.obsid) + '_grid_' + str(i) + ".out\n" +\
                    "#SBATCH --nodes=1\n" +\
                    "#SBATCH --account=mwaops\n" +\
                    "nnodes=1 # number of requested nodes - make sure this matches the --nodes SBATCH variable\n" +\
                    'nprocesses=`bc -l <<< "${nnodes}*16"` # assuming 16 parallel processes available per node\n' +\
                    "obsid="+str(args.obsid) +" # observation ID\n" +\
                    "freq="+str(centrefreq)+" # centre observing frequency in Hz\n" +\
                    "eff=0.98037 # antenna efficiency\n" +\
                    'ra='+str(ra)+' # target RA\n' +\
                    "dec="+str(dec)+" # target DEC\n" +\
                    'flags="23 71" # flagged tiles (as in RTS solutions)\n' +\
                    "tres=0.01 # resolution for ZA grid\n" +\
                    "pres=0.01 # resolution for Azimuth grid\n" +\
                    "obstime="+str(args.obsid) +" # GPS time to evaluate the beam\n" +\
                    'odir="/group/mwaops/nswainston/pabeam/output/" # where to dump the outputs\n' +\
                    'radius="1"\n' +\
                    "time aprun -n $nprocesses python /group/mwaops/nswainston/pabeam/pabeam.py -o ${obsid} -f ${freq} -t ${obstime} -e ${eff} -p \"${ra} ${dec}\" --flagged_tiles ${flags} --grid_res ${tres} ${pres} --out_dir ${odir} -r ${radius} --write --plot_centre \"" + str(ra_decs[i][0]) + " " + str(ra_decs[i][1]) + "\" --out_name " + str(args.obsid) + '_grid_' + str(i) + '\n' +\
                    'bash /group/mwaops/nswainston/pabeam/concat.sh ' + str(args.obsid) + '_grid_' + str(i) + ' ' + str(args.obsid) + '_grid_' + str(i) + '.dat 16 '


                batch_file.write(batch_line)   
            submit_line = 'sbatch Beam_calc_' + str(args.obsid) + '_grid_' + str(i) + '.batch'
            submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
            for line in submit_cmd.stdout:
                print line,
                if "Submitted" in line:
                    (word1,word2,word3,jobid) = line.split()
                    job_id_list.append(jobid)
                    
        job_id_str = ""
        for i in job_id_list:
            job_id_str += ":" + str(i)
        with open('dependancy.batch','w') as batch_file:
            batch_line = "#!/bin/bash -l\n" +\
                         "#SBATCH --job-name=dependancy\n" +\
                         "#SBATCH --output=dependancy.out\n" +\
                         "#SBATCH --export=NONE\n" +\
                         "#SBATCH --partition=gpuq\n" +\
                         "#SBATCH --time=0:05:00\n" +\
                         "#SBATCH --gid=mwaops\n" +\
                         "#SBATCH --account=mwaops\n" +\
                         "#SBATCH --nodes=1\n" +\
                         "#SBATCH --dependency=afterok" + job_id_str + "\n" +\
                         "aprun python /group/mwaops/nswainston/pabeam/grid.py " + opts_string + ' -m p'
            batch_file.write(batch_line)
            
                
        submit_line = 'sbatch dependancy.batch'
        submit_cmd = subprocess.Popen(submit_line,shell=True,stdout=subprocess.PIPE)
        for line in submit_cmd.stdout:
                print line,
        print "Sent off dependancy"
        
        print "code will continue when job completes"
    
        
#optional plot all the powers    
if args.mode == 'p' and args.plot:
    import glob
    dat_file_list = glob.glob(str(args.obsid) + '_grid_?.dat')
    dat_file_list.sort()
    print dat_file_list
    for dat_file in dat_file_list:
        if dat_file.endswith('0.dat'):
            theta,phi,beam = np.genfromtxt(dat_file,
                                    comments=("#","*"),skip_header=14,usecols=(0,1,8),unpack=True)
            
            import matplotlib
            matplotlib.use('Agg')
            fig = plt.figure(figsize=(7, 7))
            plt.ylabel("phi (radians)")
            plt.xlabel("theta (radians)")
            ax = fig.add_subplot(111)#,polar=True)
            v = np.linspace(min(beam), max(beam), 15, endpoint=True)
            #ax.tricontourf(np.radians(phi),np.radians(theta),beam,v,cmap=plt.cm.jet)
            mesh = ax.tricontourf(np.radians(phi),np.radians(theta),beam,v,cmap=plt.cm.jet)
            plt.colorbar(mesh, ax=ax)
            plt.savefig(dat_file[:-4]+'.png',bbox_inches='tight',dpi=900)
        else:                            
            theta,phi,beam_temp = np.genfromtxt(dat_file,
                                    comments=("#","*"),skip_header=14,usecols=(0,1,8),unpack=True)
            fig = plt.figure(figsize=(7, 7))
            plt.ylabel("phi (radians)")
            plt.xlabel("theta (radians)")
            ax = fig.add_subplot(111)#,polar=True)
            v = np.linspace(min(beam), max(beam), 15, endpoint=True)
            #ax.tricontourf(np.radians(phi),np.radians(theta),beam,v,cmap=plt.cm.jet)
            mesh = ax.tricontourf(np.radians(phi),np.radians(theta),beam_temp,v,cmap=plt.cm.jet)
            plt.colorbar(mesh, ax=ax)
            plt.savefig(dat_file[:-4]+'.png',bbox_inches='tight',dpi=900)
            
            print max(beam_temp)
            
            for i in range(len(theta)):
                if beam[i] < beam_temp[i]:
                    beam[i] = beam_temp[i]
                    
    print "plotting in polar coordinates"
    print "beam max: " + str(max(beam))

    fig = plt.figure(figsize=(7, 7))
    plt.ylabel("phi (radians)")
    plt.xlabel("theta (radians)")
    ax = fig.add_subplot(111)#,polar=True)
    v = np.linspace(min(beam), max(beam), 15, endpoint=True)
    #ax.tricontourf(np.radians(phi),np.radians(theta),beam,v,cmap=plt.cm.jet)
    mesh = ax.tricontourf(np.radians(phi),np.radians(theta),beam,v,cmap=plt.cm.jet)
    plt.colorbar(mesh, ax=ax)
    plt.savefig(str(args.obsid)+'_grid.png',bbox_inches='tight',dpi=900)

    """
    #BASIC TEST TO MAKE SURE DATA ARE LOADED ETC
    fig = plt.figure()
    ax = fig.add_subplot(111,polar=True)
    ax.scatter(theta,phi,c=beam,edgecolors='none')
    plt.show()
    """
    
    
    
