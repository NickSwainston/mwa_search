#! /usr/bin/env python

import subprocess
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import math

from astropy.io import fits

import mwa_metadb_utils as meta

#from pabeam import getTileLocations
def getTileLocations(obsid, flags=[], fdir="."):
    """
    Function grab the MWA tile locations for any given observation ID. Downloads the relevant metafits file from the database, saves it as <obdID>_metafits_ppds.fits.
    
    Input:
      obsid - the GPS observation ID
      flags - RTS tile flags (i.e. the first entry in the metafits correspond to "tile 0", irrespective of what the antenna name is) 
    Return:
      a list of lists containing the following:
        list[0] = a list of tile positions East of the array centre
        list[1] = a list of tile positions North of the array centre
        list[2] = a list of tile heights about sea-level 
    """

    f = fits.open('{0}/{1}_metafits_ppds.fits'.format(fdir,obsid))      
    
    east = f[1].data['East'][::2]
    north = f[1].data['North'][::2]
    height = f[1].data['Height'][::2] # height above sea-level
    
    # MWA array centre height above sea-level
    mwacentre_h = 377.827
    height = height - mwacentre_h
    
    tiles = range(128)#change if we ever get to phase 3

    # flag the tiles from the x,y,z positions
    east = np.delete(east,flags)
    north = np.delete(north,flags)
    height = np.delete(height,flags)
    tiles = np.delete(tiles,flags)

    return east,north,height,tiles

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""
            Can flag an input number of outermost tiles for an MWA observation and simulate the effect it'll have on signal to noise and FWHM
            """)
    parser.add_argument('-o','--obsid',type=str,help='The observation ID of the fits file to be searched')
    parser.add_argument("--out_dir",type=str,action='store',help="Location (full path) to write the output data files",default=".")
    parser.add_argument("--flagged_tiles",type=str,nargs='+',action='store',metavar="tile",\
                            help="The tiles flagged as in when running the RTS. Must be a list of space separated tile numbers, e.g. 0 1 2 5",default=None)
    parser.add_argument("-n","--numflag",type=int,help="Number of tiles to flag",default=10)
    args=parser.parse_args()

    if args.flagged_tiles is None:
        flags = []
    else:
        flags = args.flagged_tiles

    print "gathering required data"
    if not os.path.exists('{0}/{1}_metafits_ppds.fits'.format(args.out_dir,args.obsid)):
        os.system('wget -O {0}/{1}_metafits_ppds.fits mwa-metadata01.pawsey.org.au/metadata/fits?obs_id={1}'.format(args.out_dir, args.obsid))
    xpos, ypos, zpos, tile_n = getTileLocations(args.obsid, flags, fdir=args.out_dir)
    xpos = xpos-np.mean(xpos)
    ypos = ypos-np.mean(ypos)
    orig_tn = len(xpos)

    dis = []
    for i in range(len(xpos)):
        dis.append(np.sqrt(xpos[i]**2 + ypos[i]**2))
        if math.isnan(dis[i]):
            print xpos[i],ypos[i]

    #dec_list = [x for _,x in sorted(zip(ra_list,dec_list))]
    xpos = [x for _,x in sorted(zip(dis,xpos))]
    ypos = [x for _,x in sorted(zip(dis,ypos))]
    tile_n = [x for _,x in sorted(zip(dis,tile_n))]
    dis = sorted(dis)
    #print dis

    #get freq from metadata
    beam_meta_data = meta.getmeta(service='obs', params={'obs_id':args.obsid})
    channels = beam_meta_data[u'rfstreams'][u"0"][u'frequencies']
    minfreq = float(min(channels))
    maxfreq = float(max(channels))
    centrefreq = 1.28 * 10**6 * (minfreq + (maxfreq-minfreq)/2) #in Hz
    wavelength = 2.997*10**8 / centrefreq

    #simulate change in FWHM
    FWHMs = []
    print range(1,len(tile_n))
    print range(len(xpos)-1)
    print range(1,len(ypos)-1)
    #for fn in range(1,len(tile_n)-1):
    for fn in range(1,60):
        temp_xpos = xpos[:-fn]
        temp_ypos = ypos[:-fn]

        #loop through for baselines calcs
        baselines = []
        for i in range(len(temp_xpos)):
            for j in range(i+1,len(temp_ypos)):#does this need a plus 1?
                temp_base = math.sqrt( (temp_xpos[i] - temp_xpos[j])**2 +\
                                           (temp_ypos[i] - temp_ypos[j])**2 )
                #print temp_base
                baselines.append(temp_base)
        #print baselines
        FWHMs.append(wavelength/max(baselines)/3.14*180)
    
    plt.clf()
    plt.plot(FWHMs)
    plt.xlabel("number of tiles flagged (furthest first)")
    plt.ylabel("FHWM in deg")
    plt.show()

    #remove N furthes tiles
    flagged_xpos = xpos[-(args.numflag):]
    flagged_ypos = ypos[-(args.numflag):]
    flagged_tiles = tile_n[-(args.numflag):]
    print "Tiles to flag: ",flagged_tiles

    xpos = xpos[:-(args.numflag)]
    ypos = ypos[:-(args.numflag)]
    dis = dis[:-(args.numflag)]

    plt.scatter(xpos,ypos, s=0.2)
    plt.scatter(flagged_xpos,flagged_ypos, s=0.2)
    plt.show()

    #calculate relative signal to noise
    tile_n_range = range(1,orig_tn)
    sn_range = []
    for t in tile_n_range:
        sn_range.append(np.sqrt(t*(t-1))/orig_tn)

    plt.clf()
    plt.plot(tile_n_range, sn_range)
    plt.show()

    print "Will have a reltaive sensitivity of "+str(sn_range[-(args.numflag)])

   
