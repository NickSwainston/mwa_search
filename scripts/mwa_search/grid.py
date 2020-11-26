#!/usr/bin/env python3
import argparse
from astropy.coordinates import SkyCoord,EarthLocation,AltAz
from astropy.time import Time
import astropy.units as u
from math import cos
import numpy as np

from matplotlib import use
use('Agg')
import matplotlib.pyplot as plt
from matplotlib import patches

from mwa_pb.mwa_tile import h2e

# vcstools imports
import mwa_metadb_utils as meta
import find_pulsar_in_obs as fpio
from find_pulsar_in_obs import get_psrcat_ra_dec
from vcstools.pointing_utils import sex2deg, format_ra_dec

# mwa_search imports
from mwa_search.obs_tools import getTargetAZZA
from mwa_search.grid_tools import get_grid

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    Makes a hexogonal grid pattern around a pointing for a MWA VCS observation.
    grid.py -o 1166459712 -p "06:25:31.20 -36:40:48.0" -d 0.6 -l 1
    """)
    parser.add_argument('-o', '--obsid',type=str,help='Observation ID')
    parser.add_argument('-p', '--pointing',type=str,help='Centre pointing in hh:mm:ss.ss_dd\"mm\'ss.ss')
    parser.add_argument('--aitoff',action="store_true",help='Plots the output in aitoff (may make it hard to analyise).')
    parser.add_argument('-f', '--fraction',type=float,help='Fraction of the full width half maximum to use as the distance between beam centres',default=0.85)
    parser.add_argument('-d', '--deg_fwhm',type=float,help='Sets the FWHM at zenith in degrees (best to test along dec). The script will not calculate the FWHM',default=0.3098)
    parser.add_argument('--dec_range_fwhm',type=float,nargs='+',help='A list of FWHM and ranges in the order of: "FWHM1 decmin1 decmax1 FWHM2 decmin2 decmax2"')
    parser.add_argument('-t', '--type',type=str,help='Can be put in either "hex" or "square" tiling mode. Default is hex.',default='hex')
    parser.add_argument('-l', '--loop',type=int,help='Number  of "loops" around the centre pointing the code will calculate. Default is 1',default=1)
    parser.add_argument('--fill',type=float,help='Calculate the number of loops required to fill a circle of the input radius in degrees.')
    parser.add_argument('-a','--all_pointings',action="store_true",help='Will calculate all the pointings within the FWHM of the observations tile beam.')
    parser.add_argument('-b', '--begin',type=int,help='Begin time of the obs for the --all_pointings options')
    parser.add_argument('-e', '--end',type=int,help='End time of the obs for the --all_pointings options')
    parser.add_argument('--dec_range',type=float,nargs='+',help='Dec limits: "decmin decmax". Default -90 90', default=[-90,90])
    parser.add_argument('--ra_range',type=float,nargs='+',help='RA limits: "ramin ramax". Default 0 360', default=[0,360])
    parser.add_argument('-v','--verbose_file',action="store_true",help='Creates a more verbose output file with more information than make_beam.c can handle.')
    parser.add_argument('--pulsar',type=str,nargs='+',help='A list of pulsar to mark on the plot')
    parser.add_argument('-n', '--n_pointings', type=int, default=None, help='Number of pointings per output file.')
    parser.add_argument('--out_file_name', type=str, help='The output file name.')

    args=parser.parse_args()

    opts_string = ""
    for k in args.__dict__:
        if args.__dict__[k] is not None:
            if k == "pointing":
                opts_string = opts_string + ' --' + str(k) + ' "' + str(args.__dict__[k][0]) +\
                         ' ' + str(args.__dict__[k][1]) + '"'
            else:
                opts_string = opts_string + ' --' + str(k) + ' ' + str(args.__dict__[k])

    if args.obsid:
        obs, ra, dec, duration, xdelays, centrefreq, channels = \
                meta.get_common_obs_metadata(args.obsid)

    #get fwhm in radians
    centre_fwhm = np.radians(args.deg_fwhm)

    #all_pointing parsing
    if (args.loop != 1) and args.all_pointings:
        print("Can't use --loop and --all_poinitings as all_pointings calculates the "
              "loops required. Exiting.")
        quit()
    if (args.loop != 1) and args.fill:
        print("Can't use --loop and --fill as --fill calculates the "
              "loops required. Exiting.")
        quit()
    if args.pointing and args.all_pointings:
        print("Can't use --pointing and --all_poinntings as --all_pointings calculates "
              "the pointing. Exiting.")
        quit()
    if args.pointing and args.pulsar:
        print("Can't use --pointing and --pulsar as --pulsar calculates the pointing. Exiting.")
        quit()
    if args.pulsar and args.all_pointings:
        print("Can't use --pulsar and --all_poinntings as --all_pointings calculates "
              "the pointing. Exiting.")
        quit()

    if args.fill:
        args.loop = int( (args.fill - args.deg_fwhm/2.) / (args.deg_fwhm*args.fraction) )
        print("Using {} loops to fill {} degrees".format(args.loop, args.fill ))
    #calculate pointing
    if args.all_pointings:
        #calculating loop number
        fudge_factor = 2.
        tile_fwhm = np.degrees(fudge_factor * (3*10**8/(centrefreq*10**6))/6.56 )
        #account for the "increase" in tile beam size due to drifting
        tile_fwhm += duration/3600.*15.
        args.loop = int(tile_fwhm/2./(args.deg_fwhm*args.fraction))

        #calculating pointing from metadata
        if int(obs) < 1252177700 :
            #the ra used to mean the start of the obs so it had to be corrected for
            ra = np.radians(ra + duration/3600.*15./2)
        else:
            ra = np.radians(ra)
        dec = np.radians(dec)
    elif args.pulsar:
        temp = fpio.get_psrcat_ra_dec(pulsar_list=args.pulsar)
        _, raj, decj = format_ra_dec(temp, ra_col = 1, dec_col = 2)[0]
        coord = SkyCoord(raj, decj, unit=(u.hourangle,u.deg))
        ra = coord.ra.radian #in radians
        dec = coord.dec.radian
    elif args.pointing:
        coord = SkyCoord(args.pointing.split("_")[0],args.pointing.split("_")[1],
                         unit=(u.hourangle,u.deg))
        ra = coord.ra.radian #in radians
        dec = coord.dec.radian
    else:
        print("Please use either --pointing, --pulsar or --all_pointings. Exiting.")
        quit()

    #calculate grid
    rads, decds = get_grid(ra, dec, centre_fwhm*args.fraction, args.loop, grid_type=args.type)

    #remove pointings outside of ra or dec range
    if args.dec_range != [-90,90] or args.ra_range != [0, 360]:
        print("Removing pointings outside of ra dec ranges")
        radls = []
        decdls = []
        for i in range(len(rads)):
            if  (args.dec_range[0] < float(decds[i]) < args.dec_range[1] ) and \
                (args.ra_range[0]  < float(rads[i]) < args.ra_range[1]):
                    radls.append(rads[i])
                    decdls.append(decds[i])
        rads = radls
        decds = decdls

    if args.all_pointings:
        #calculate powers
        obeg, oend = meta.obs_max_min(obs)
        if args.begin:
            start_time = obeg - args.begin
        else:
            start_time = 0
        if args.end and args.begin:
            duration = args.end - args.begin
        elif args.end:
            duration = args.end - obeg
        obs_metadata = [obs, ra, dec, duration, xdelays, centrefreq, channels]
        names_ra_dec = []
        for ni in range(len(rads)):
            if float(decds[ni]) < -90.:
                continue
            names_ra_dec.append(["name", rads[ni], decds[ni]])
        names_ra_dec = np.array(names_ra_dec)
        power = fpio.get_beam_power_over_time(obs_metadata,
                                              names_ra_dec, degrees=True)

        #check each pointing is within the tile beam
        radls = []
        decdls = []
        tFWHM = np.amax(power)/2. #assumed half power point of the tile beam
        for ni in range(len(names_ra_dec)):
            if max(power[ni]) > tFWHM:
                radls.append(rads[ni])
                decdls.append(decds[ni])
        rads = radls
        decds = decdls

    print("Using skycord to convert ra dec")
    #Use skycoord to get asci
    coord = SkyCoord(rads,decds,unit=(u.deg,u.deg))
    #unformated
    rags_uf = coord.ra.to_string(unit=u.hour, sep=':')
    decgs_uf = coord.dec.to_string(unit=u.degree, sep=':')

    ras = []; decs = []; theta = []; phi = []
    time = Time(float(args.obsid),format='gps')
    print("Formating the outputs")
    #format the ra dec strings
    for i in range(len(rags_uf)):
        rag = rags_uf[i]
        decg = decgs_uf[i]

        temp = format_ra_dec([[rag,decg]])
        rag = temp[0][0]
        decg = temp[0][1]

        if args.verbose_file:
            az,za,azd,zad = getTargetAZZA(rag,decg,time)
        else:
            az,za,azd,zad = [0,0,0,0]

        ras.append(rag)
        decs.append(decg)
        theta.append(az)
        phi.append(za)

    if args.out_file_name:
        out_file_name = args.out_file_name
    else:
        if args.obsid:
            out_file_name = str(args.obsid)
        else:
            out_file_name = ''
        if args.pulsar:
            out_file_name = '{0}_{1}'.format(out_file_name, args.pulsar[0])
        out_file_name += '_grid_positions'
        if args.dec_range != [-90,90] or args.ra_range != [0, 360]:
            out_file_name += '_ra_dec_limited'
        out_file_name = '{0}_f{1}_d{2}_l{3}'.format(out_file_name, args.fraction,
                                                    args.deg_fwhm, args.loop)

    #Writing file
    if args.n_pointings is None:
        print("Recording the dec limited positons in {0}.txt".format(out_file_name))
        with open('{0}.txt'.format(out_file_name),'w') as out_file:
            if args.verbose_file:
                out_line = "#ra   dec    az     za\n"
                out_file.write(out_line)
            for i in range(len(rads)):
                if args.verbose_file:
                    out_line = str(ras[i])+" "+str(decs[i])+" "+str(theta[i])+" "\
                                +str(phi[i])+" "+str(rads[i])+" "\
                                +str(decds[i])+"\n"
                else:
                    out_line = str(ras[i])+"_"+str(decs[i])+"\n"
                out_file.write(out_line)
    else:
        ra_chunks = [ras[x:x+args.n_pointings] for x in range(0, len(ras), args.n_pointings)]
        dec_chunks = [decs[x:x+args.n_pointings] for x in range(0, len(decs), args.n_pointings)]
        for ci in range(len(ra_chunks)):
            first_id = ci * args.n_pointings + 1
            last_id  = ci * args.n_pointings + len(ra_chunks[ci])
            print("Recording the dec limited positons in {0}_{1}_{2}.txt".format(out_file_name, first_id, last_id))
            with open('{0}_{1}_{2}.txt'.format(out_file_name, first_id, last_id),'w') as out_file:
                for i in range(len(ra_chunks[ci])):
                    out_file.write("{0}_{1}\n".format(ra_chunks[ci][i], dec_chunks[ci][i]))

    #matplotlib.use('Agg')
    print("Plotting")
    fig = plt.figure(figsize=(7, 7))
    if args.aitoff:
        fig.add_subplot(111)
        print("changing axis")
        ax = plt.axes(projection='mollweide')
        rads = -(np.radians(np.array(rads)))+ np.pi
        decds = np.radians(np.array(decds))
    else:
        plt.axes().set_aspect('equal')
        ax = plt.gca()
        #ax.axis([325., 345., -9., 0.])

    plt.xlabel("ra (degrees)")
    plt.ylabel("dec (degrees)")

    for i in range(len(ras)):
        if args.aitoff:
            fwhm_circle = centre_fwhm/cos(decds[i]) / 2.
            circle = plt.Circle((rads[i],decds[i]),fwhm_circle,
                                color='r', lw=0.1,fill=False)
            ax.add_artist(circle)
        else:
            fwhm_vert = np.degrees(centre_fwhm/cos(np.radians(decds[i] + 26.7))**2)
            fwhm_horiz = np.degrees(centre_fwhm/cos(np.radians(decds[i])) )

            ellipse = patches.Ellipse((rads[i],decds[i]), fwhm_horiz, fwhm_vert,
                                          linewidth=0.3, fill=False, edgecolor='green')
            ax.add_patch(ellipse)
            #fwhm_circle = centre_fwhm/cos(np.radians(decds[i])) / 2.
            #circle = plt.Circle((rads[i],decds[i]),np.degrees(fwhm_circle),
            #                     color='r', lw=0.1,fill=False)
    plt.scatter(rads,decds,s=0.1,c='black')

    #add some pulsars
    if args.pulsar:
        ra_PCAT = []
        dec_PCAT = []
        pulsar_list = get_psrcat_ra_dec(pulsar_list = args.pulsar)
        for pulsar in pulsar_list:
            ra_temp, dec_temp = sex2deg(pulsar[1], pulsar[2])
            ra_PCAT.append(ra_temp)
            dec_PCAT.append(dec_temp)
        ax.scatter(ra_PCAT, dec_PCAT, s=15, color ='r', zorder=100)

    plt.savefig('{0}.png'.format(out_file_name), bbox_inches='tight', dpi =1000)




    print("Number of pointings: " + str(len(rads)))
    #times out and segfaults after this so I'm going to exit here
    exit()
