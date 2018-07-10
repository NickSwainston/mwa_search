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
    parser.add_argument("--flag_file",type=str,action='store',help="Location of a flagged_tiles.txt file to alter to include outer tiles.")
    parser.add_argument("--flagged_tiles",type=str,nargs='+',action='store',metavar="tile",\
                            help="The tiles flagged as in when running the RTS. Must be a list of space separated tile numbers, e.g. 0 1 2 5",default=None)
    parser.add_argument("-n","--numflag",type=int,help="Number of tiles to flag",default=10)
    parser.add_argument("-p","--plot",action='store_true',help="Creates several plots")
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
    print "Tiles from furtherst to closest: ",
    for i in reversed(tile_n):
        print i,
    print "\n",
    #get freq from metadata
    beam_meta_data = meta.getmeta(service='obs', params={'obs_id':args.obsid})
    channels = beam_meta_data[u'rfstreams'][u"0"][u'frequencies']
    minfreq = float(min(channels))
    maxfreq = float(max(channels))
    centrefreq = 1.28 * 10**6 * (minfreq + (maxfreq-minfreq)/2) #in Hz
    print centrefreq
    wavelength = 2.997*10**8 / centrefreq

    #simulate change in FWHM
    FWHMs = []
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
        FWHMs.append(np.degrees(1.22*wavelength/max(baselines))*60.)
    print FWHMs[0]
    print FWHMs[15]
    
       #remove N furthes tiles
    flagged_xpos = xpos[-(args.numflag):]
    flagged_ypos = ypos[-(args.numflag):]
    flagged_tiles = tile_n[-(args.numflag):]
    print "Tiles to flag: ",flagged_tiles
    
    if args.flag_file:
        with open(args.flag_file, 'a+') as f:
            lines = f.readlines()
            for ft in flagged_tiles:
                if not (str(ft)+"\n") in lines:
                    f.write(str(ft)+"\n")

    xpos = xpos[:-(args.numflag)]
    ypos = ypos[:-(args.numflag)]
    dis = dis[:-(args.numflag)]

    if args.plot:
        from mpl_toolkits.axes_grid1 import host_subplot
        import mpl_toolkits.axisartist as AA
        import matplotlib.pyplot as plt
        
       
        #calculate relative signal to noise
        tile_n_range = range(1,orig_tn)
        sn_range = []
        for t in tile_n_range:
            sn_range.append(1.-np.sqrt(t*(t-1))/orig_tn)
                
        """
        host = host_subplot(111, axes_class=AA.Axes)
        plt.subplots_adjust(right=0.75)

        par1 = host.twinx()

        offset = 60


        host.set_xlim(0, 2)
        host.set_ylim(0, 2)

        host.set_xlabel("Distance")
        host.set_ylabel("Density")
        par1.set_ylabel("Temperature")

        p1, = host.plot([0, 1, 2], [0, 1, 2], label="Density")
        p2, = par1.plot([0, 1, 2], [0, 3, 2], label="Temperature")

        par1.set_ylim(0, 4)

        host.legend()

        host.axis["left"].label.set_color(p1.get_color())
        par1.axis["right"].label.set_color(p2.get_color())

        plt.draw()
        plt.show()
        
        """
        plt.clf()
        host = host_subplot(111, axes_class=AA.Axes)
        par1 = host.twinx()

        offset = 0
        new_fixed_axis = par1.get_grid_helper().new_fixed_axis
        par1.axis["right"] = new_fixed_axis(loc="right",
                                            axes=par1,
                                            offset=(offset, 0))

        par1.axis["right"].toggle(all=True)

        
        host.set_ylim(min(FWHMs), max(FWHMs))
        par1.set_ylim(min(sn_range[:len(FWHMs)]), max(sn_range[:len(FWHMs)]))
        
        if int(args.obsid) == 1117643248 and centrefreq == 118400000.0:
            #use some simulated values from pabeam.py
            
            host.set_ylim(0.014667480523, max(FWHMs))
            p3, = host.plot([0,1,3,5,9,12,16],np.array([0.014667480523,0.0151320483657,0.0157325732253,\
                            0.0162094738779,0.01782310605,0.0194199249546,0.0198223907006])*60.,\
                            label="Simulated FWHM")


        #plt.plot(FWHMs)
        host.set_xlabel("Number of tiles flagged (furthest first)")
        host.set_ylabel("FHWM in arcmin")
        par1.set_ylabel("Relative signal to noise")
        
        print len(tile_n_range),len(FWHMs), len(sn_range)
        p1, = host.plot(tile_n_range[:len(FWHMs)], FWHMs, label="Calculated FWHM")
        p2, = par1.plot(tile_n_range[:len(FWHMs)], sn_range[:len(FWHMs)], label="Relative S/N")
        

        host.legend()

        host.axis["left"].label.set_color(p1.get_color())
        par1.axis["right"].label.set_color(p2.get_color())
        
        plt.draw()
        plt.savefig("outer_tiles_flagged_plot.png")
        plt.show()
        

        #plot tile positions
        dot_size = 4
        plt.scatter(xpos,ypos, s=dot_size)
        plt.scatter(flagged_xpos,flagged_ypos, color='r',s=dot_size)
        plt.show()
 


    print "Will have a reltaive sensitivity of "+str(sn_range[args.numflag])

   
