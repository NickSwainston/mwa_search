#! /usr/bin/env python3

import os
import math
import argparse
import numpy as np
import csv
from scipy.interpolate import UnivariateSpline
from math import radians, degrees

#vcstools
from vcstools.beam_calc import get_beam_power_over_time
from vcstools.catalogue_utils import get_psrcat_ra_dec
from vcstools.pointing_utils import sex2deg, deg2sex
from vcstools.metadb_utils import find_obsids_meta_pages, get_common_obs_metadata

#matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import rcParams
rcParams['font.family'] = 'monospace'

def SMART_obs_calc(degree_overlap, manual_overlap):
    """
    Work out how many observations are required to cover the southern sky
    """

    #setting up the dec ranges
    dec_range = [-72., -55., -40.5, -26.7, -13., +1.6, +18.3] #Gleam pointings
    delays_range = [[0,0,0,0,6,6,6,6,12,12,12,12,18,18,18,18],\
                    [0,0,0,0,4,4,4,4,8,8,8,8,12,12,12,12],\
                    [0,0,0,0,2,2,2,2,4,4,4,4,6,6,6,6],\
                    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],\
                    [6,6,6,6,4,4,4,4,2,2,2,2,0,0,0,0],\
                    [12,12,12,12,8,8,8,8,4,4,4,4,0,0,0,0],\
                    [18,18,18,18,12,12,12,12,6,6,6,6,0,0,0,0]]

    print("Using GLEAM dec range: {}".format(dec_range))
    """
    sweet_dec_range = [-82.8,-71.4,-63.1,-55.,-47.5,-40.4,-33.5,-26.7,-19.9,-13.,-5.9,1.6,9.7,18.6,29.4,44.8]
    sweet_delays_range= [[0,0,0,0,7,7,7,7,14,14,14,14,21,21,21,21],\
                         [0,0,0,0,6,6,6,6,12,12,12,12,18,18,18,18],\
                         [0,0,0,0,5,5,5,5,10,10,10,10,15,15,15,15],\
                         [0,0,0,0,4,4,4,4,8,8,8,8,12,12,12,12],\
                         [0,0,0,0,3,3,3,3,6,6,6,6,9,9,9,9],\
                         [0,0,0,0,2,2,2,2,4,4,4,4,6,6,6,6],\
                         [0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3],\
                         [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],\
                         [3,3,3,3,2,2,2,2,1,1,1,1,0,0,0,0],\
                         [6,6,6,6,4,4,4,4,2,2,2,2,0,0,0,0],\
                         [9,9,9,9,6,6,6,6,3,3,3,3,0,0,0,0],\
                         [12,12,12,12,8,8,8,8,4,4,4,4,0,0,0,0],\
                         [15,15,15,15,10,10,10,10,5,5,5,5,0,0,0,0],\
                         [18,18,18,18,12,12,12,12,6,6,6,6,0,0,0,0],\
                         [21,21,21,21,14,14,14,14,7,7,7,7,0,0,0,0],\
                         [24,24,24,24,16,16,16,16,8,8,8,8,0,0,0,0]]

    dec_range = []
    delays_range =[]
    sweet_spots_range = [0,2,4,7,10,12,14]
    for i in sweet_spots_range:
      dec_range.append(sweet_dec_range[i])
      delays_range.append(sweet_delays_range[i])
    print dec_range
    """

    #Going to work out how many pointings are needed
    #setting up some metadata requirements
    time = 4800 #one hour 20 min
    channels = range(107,131)
    minfreq = float(min(channels))
    maxfreq = float(max(channels))
    centrefreq = 1.28 * (minfreq + (maxfreq-minfreq)/2) #in MHz

    start_obsid = '1117624530'
    start_ra = 180.
    Dec_FWHM_calc = []
    RA_FWHM_calc = []
    for i in range(-89,89,1):
        for j in range(0,361,1):
            Dec_FWHM_calc.append(i)
            RA_FWHM_calc.append(j)

    observations = []
    ra_list =[]
    dec_list =[]
    delays_list = []
    FWHM = []
    FWHM_Dec = []
    pointing_count = 0
    for i in range(len(dec_range)):
        #calculating the FWHM at this dec
        ra_sex, deg_sex = deg2sex(start_ra, dec_range[i])
        cord = [start_obsid, str(ra_sex), str(deg_sex), 1, delays_range[i],centrefreq, channels]
        #powout=get_beam_power(cord, zip(RA_FWHM_calc,Dec_FWHM_calc), dt=600)
        names_ra_dec = np.column_stack((['source']*len(RA_FWHM_calc), RA_FWHM_calc, Dec_FWHM_calc))
        powout = get_beam_power_over_time(cord, names_ra_dec, dt=600, degrees = True)
        powout_RA_line = []
        powout_Dec_line = []
        RA_line = []
        Dec_line = []
        for p in range(len(powout)):
            #print(int(y[i]/np.pi*180.), int(dec) )
            if int(Dec_FWHM_calc[p]) == int(dec_range[i]):
                powout_RA_line.append(float(powout[p]))
                RA_line.append(float(RA_FWHM_calc[p]))
            if int (RA_FWHM_calc[p]) == int(start_ra):
                powout_Dec_line.append(float(powout[p]))
                Dec_line.append(float(Dec_FWHM_calc[p]))

        print("\nValues for Dec " + str(dec_range[i]))
        #work out RA FWHM (not including the drift scan, 0sec observation)
        if args.fwhm:
            spline = UnivariateSpline(RA_line, powout_RA_line-np.max(powout_RA_line)/2., s=0)
        else:
            spline = UnivariateSpline(RA_line, powout_RA_line-np.full(len(powout_RA_line),0.5), s=0)
        try:
            r1, r2 = spline.roots()
        except ValueError:
            print("No FWHM for " + str(dec_range[i]) + " setting to 1000 to skip")
            FWHM.append(1000.)
            pointing_count -=1
        else:
            FWHM.append(float(r2-r1))
            print("FWHM along RA at dec "+ str(dec_range[i]) + ": " + str(FWHM[i]))

        #work out Dec FWHM
        if args.fwhm:
            spline = UnivariateSpline(Dec_line, powout_Dec_line-np.max(powout_Dec_line)/2., s=0)
            r1, r2 = spline.roots()
            FWHM_Dec.append(float(r2-r1))
            print("FWHM along Dec at dec "+ str(dec_range[i]) + ": " + str(FWHM_Dec[i]))

        deg_move = total_angle = FWHM[i] - degree_overlap*math.cos(math.radians(dec_range[i])) + \
                    float(time)/3600.*15.*math.cos(math.radians(dec_range[i]))
        if manual_overlap is not None:
            point_num_this_deg = manual_overlap[i]
        else:
            point_num_this_deg = int(360./deg_move) + 1
        print("Number for this dec: " +str(point_num_this_deg))
        deg_move = 360. / point_num_this_deg
        overlap_true = FWHM[i] + float(time)/3600.*15.*math.cos(math.radians(dec_range[i])) -\
                       360./point_num_this_deg
        print("True overlap this dec: " + str(overlap_true))

        # offset every second dec range by half a FWHM in RA
        for x in range(point_num_this_deg):
            if i % 2 == 0:
                temp_ra = start_ra + x * deg_move
                observations.append(str(int(start_obsid) + int(x*deg_move*240)))
            else:
                temp_ra = start_ra + x * deg_move +\
                          deg_move / math.cos(math.radians(dec_range[i]))
                observations.append(str(int(start_obsid) + int(x*deg_move*240) +\
                                        int(deg_move*120)))
            if temp_ra > 360.:
               temp_ra = temp_ra -360.
            ra_list.append(temp_ra)
            dec_list.append(dec_range[i])
            delays_list.append(delays_range[i])
            total_angle += deg_move
            pointing_count+=1

    #Sort by ra
    dec_list =     [x for _,x in sorted(zip(ra_list,dec_list))]
    delays_list =  [x for _,x in sorted(zip(ra_list,delays_list))]
    observations = [x for _,x in sorted(zip(ra_list,observations))]
    ra_list = sorted(ra_list)

    return observations, dec_list, ra_list, delays_list



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    A ploting script tha can be used to plot MWA tile beams, pulsars and used to work out the SMART observations to best cover the southern sky.

    #SMART survey update example
    plot_obs_pulsar.py -m 6 9 11 11 11 11 11 --smart -f --contour --pulsar J2241-5236 J2145-0750 J2222-0137 J2248-0101 J2330-2005 J0034-0721 J0133-6957 J2324-6054 J0206-4028 J0051+0423
    #SMART sensitivity example
    plot_obs_pulsar.py -m 6 9 11 11 11 11 11 --sens --smart -f
    #All MWA tile beam example
    plot_obs_pulsar.py -f --contour --all_obsids
    """)
    obs_group = parser.add_argument_group('Observation Options')
    obs_group.add_argument('--obsid_list', type=str, nargs='*',
                           help='Instead of calculating which positions to use the script will use the input obsids. eg: "1088850560 1090249472"')
    obs_group.add_argument('--all_obsids', action='store_true',
                           help='Uses all VCS obsids on the MWA metadatabase.')
    obs_group.add_argument('--incoh', action='store_true',
                           help='Calculates the sensitivity assuming thast all observations are incoherent')
    obs_group.add_argument('--smart', action='store_true',
                           help='Cover the Southern sky with observations for the SMART survey.')
    obs_group.add_argument('-m', '--manual', nargs='+', type=int,
                           help='Used with the --smart option to manually decide the obs at each declination, input them as 1 2 3 4 5 6 7')
    obs_group.add_argument('-d', '--degree_overlap', type=float, default=10.,
                           help='Used with the --smart option to manually set degrees overlap in RA of the observations')

    obs_plot_group = parser.add_argument_group('Observation Plot Types Options')
    obs_plot_group.add_argument('-s', '--sens', action='store_true',
                                help='Plots sensitivity of each observation')
    obs_plot_group.add_argument('-o', '--overlap', action='store_true',
                                help='Plots sensitivity by summing the power of overlapping observations.')
    obs_plot_group.add_argument('-c', '--contour', action='store_true',
                                help='Plots the contour of each observations power')
    obs_plot_group.add_argument('-l', '--lines', action='store_true',
                                help='Includes the min decs of other telescopes in plots')

    add_group = parser.add_argument_group('Extra Plot Layers Options')
    add_group.add_argument('--pulsar', type=str, nargs='+',
                           help='A list of pulsar to mark on the plot.')
    add_group.add_argument('--pulsar_detected', action='store_true',
                           help='Plots all pulsars detected by the MWA and uploaded to the pulsar database.')
    add_group.add_argument('--pulsar_all', action='store_true',
                           help='Plots all known pulsars.')
    add_group.add_argument('--pulsar_discovered', action='store_true',
                           help='Plots all pulsars discovered with the MWA.')
    add_group.add_argument('--fill', action='store_true',
                           help='Shades the area the MWA can view.')
    add_group.add_argument('--shade', type=str, nargs='+',
                           choices=['red','green','purple','darkorange','blue'],
                           help='Shades the chosen colour observations group')
    add_group.add_argument('--shade_temp', nargs='*', default='black',
                           help='First input is the colour and then the observation IDs you want to shade')

    plot_group = parser.add_argument_group('Plotting Options')
    plot_group.add_argument('-f', '--fwhm', action='store_true',
                            help='if this options is used the FWHM of each pointing is used. If it is not chosen the FWHM of a zenith pointing is used.')
    plot_group.add_argument('-r', '--resolution', type=int, default=3,
                            help='The resolution in degrees of the final plot (must be an integer). Default = 1')
    plot_group.add_argument('--square', action='store_true',
                            help='Plots a square grid instead of aitoff.')
    plot_group.add_argument('-p', '--plot_type', type=str,
                            help='Determines the output plot type, Default="png".',default='png')
    plot_group.add_argument('--ra_offset', action='store_true',
                            help='Offsets the RA by 180 so that 0h is in the centre')
    args=parser.parse_args()

    #Setting up some of the plots
    fig = plt.figure(figsize=(6, 4))
    plt.rc("font", size=8)
    if args.square:
        fig.add_subplot(111)
        ax = plt.axes()
    else:
        fig.add_subplot(111)
        ax = plt.axes(projection='mollweide')

    #levels = np.arange(0.25, 1., 0.05)
    colors= ['0.5' for _ in range(50)] ; colors[0]= 'blue'
    linewidths= [0.4 for _ in range(50)] ; linewidths[0]= 1.0
    alpha = 0.5


    #setting up some default metadata requirements
    time = 4800 #one hour 20 min
    channels = range(107,131)
    minfreq = float(min(channels))
    maxfreq = float(max(channels))
    centrefreq = 1.28 * (minfreq + (maxfreq-minfreq)/2) #in MHz


    #setting up RA Dec ranges for power calculations
    res = args.resolution
    map_dec_range = range(-90,91,res)
    map_ra_range = range(0,361,res)
    RA=[]; Dec=[]; x = []; y = []
    for i in map_dec_range:
        for j in map_ra_range:
            Dec.append(i)
            RA.append(j)
    for c in range(len(RA)):
        if args.ra_offset:
            if RA[c] > 180:
                x.append(-RA[c]/180.*np.pi+2*np.pi)
            else:
                x.append(-RA[c]/180.*np.pi)
        else:
            x.append(-RA[c]/180.*np.pi +np.pi)
        y.append(Dec[c]/180.*np.pi)
    ny = np.array(y)
    nx = np.array(x)

    #Working out the observations required -----------------------------------------------
    if args.all_obsids:
        observations = find_obsids_meta_pages(params={'mode':'VOLTAGE_START','cenchan':145})
        #print(observations)
        #observations = [o for o in observations if o <=  1104109216]
        #print(observations)
    elif args.obsid_list:
        observations = args.obsid_list
        pointing_count = len(observations)
    elif args.smart:
        observations, dec_list, ra_list, delays_list = SMART_obs_calc(args.degree_overlap, args.manual)
    else:
        print("No observation options selected. No observations will be plotted")
        observations = []
    pointing_count = len(observations)


    nz_sens_overlap = np.zeros(len(RA))
    nz_shade_colour = {'red'        : np.zeros(len(RA)),
                       'green'      : np.zeros(len(RA)),
                       'purple'     : np.zeros(len(RA)),
                       'darkorange' : np.zeros(len(RA)),
                       'blue'       : np.zeros(len(RA))}
    nz_shade_colour_temp = np.zeros(len(RA))
    #nz_sens = np.full(len(RA), 50.)
    nz_sens = np.zeros(len(RA))
    nz_sens[:] = np.nan
    max_ra_list = []
    RA_FWHM_atdec =[]

    #Print colour group files
    if args.smart:
        colour_groups = ['red','green','purple','darkorange','blue']
        for c in range(len(colour_groups)):
            f = open(str(colour_groups[c]) + '_group_file.txt','w')
            f.write('RA\tDec\n')
            f.close()


    # a little hack to save metadata to speed up repeated calls
    if args.obsid_list or args.all_obsids:
        dec_list = []
        ra_list = []
        delays_list = []
        if not os.path.exists('obs_meta.csv'):
            os.mknod('obs_meta.csv')
        with open('obs_meta.csv', 'r') as csvfile:
            spamreader = csv.reader(csvfile)
            next(spamreader, None)
            obsid_meta_file = []
            for row in spamreader:
                obsid_meta_file.append(row)
        for i, ob in enumerate(observations):
            obs_foun_check = False
            for omf in obsid_meta_file:
                if int(ob) == int(omf[0]):
                    print("getting obs metadata from obs_meta.csv")
                    ob, ra, dec, time, delays,centrefreq, channels = omf
                    ob = int(ob)
                    time = int(time)
                    delaysx = list(map(int,delays[2:-2].split("], [")[0].split(",")))
                    delaysy = list(map(int,delays[2:-2].split("], [")[1].split(",")))
                    delays = [delaysx, delaysx]
                    centrefreq = float(centrefreq)
                    channels = list(map(int,channels[1:-1].split(",")))
                    obs_foun_check = True
            if not obs_foun_check:
                print("Getting metadata for {}".format(ob))
                ob, ra, dec, time, delays,centrefreq, channels =\
                    get_common_obs_metadata(ob)

                with open('obs_meta.csv', 'a') as csvfile:
                    spamwriter = csv.writer(csvfile)
                    spamwriter.writerow([ob, ra, dec, time, delays,centrefreq, channels])
            cord = [ob, ra, dec, time, delays,centrefreq, channels]
            ra_list.append(ra)
            dec_list.append(dec)
            delays_list.append(delays)

    #Loop over observations and calc beam power
    for i, ob in enumerate(observations):
        print("Calculating obs {0}/{1}".format(i + 1, len(observations)))
        ra = ra_list[i]
        dec = dec_list[i]
        delays = delays_list[i]

        cord = [ob, ra, dec, time, delays, centrefreq, channels]
        z=[] ; z_sens =[]

        #print(max(Dec), min(RA), Dec.dtype)
        time_intervals = 600 # seconds
        names_ra_dec = np.column_stack((['source']*len(RA), RA, Dec))
        powout = get_beam_power_over_time(cord, names_ra_dec, dt=time_intervals, degrees = True)
        #grab a line of beam power for the pointing declination
        #if i == 0:
        #    print("len powers list: " + str(powout.shape))
        for c in range(len(RA)):
            temppower = 0.
            temppower_sense = 0.
            for t in range(powout.shape[1]):
                power_ra = powout[c,t,0]
                temppower_sense += power_ra #average power kinds
                nz_sens_overlap[c] += power_ra * math.cos(ny[c])
                if power_ra > temppower:
                    temppower = power_ra
            z_sens.append(temppower_sense)
            z.append(temppower)

        nz=np.array(z)

        #calculates sensitiviy and removes zeros -------------------------
        nz_sense_obs = []
        for zsi in range(len(z_sens)):
            if nz[zsi] < 0.001:
                nz_sense_obs.append(np.nan)
            else:
                nz_sense_obs.append(4.96/np.sqrt(z_sens[zsi]))

        for zi, zs in enumerate(nz_sense_obs):
            if math.isnan(nz_sens[zi]):
                nz_sens[zi] = zs
            elif nz_sens[zi] > zs:
                #append if larger
                nz_sens[zi] = zs

        if args.fwhm:
            levels = np.arange(0.5*max(nz), max(nz), 0.5/6.)
        else:
            levels = np.arange(0.5, 1., 0.05)

        # Fill group files ------------------------------------------
        if args.smart:
            #find middle ra for each pointing
            powout_RA_line = []
            RA_line = []
            for p in range(len(nz)):
                if args.ra_offset:
                    if abs(ny[p]*180/np.pi + 0.001 - dec) < 0.5*float(res):
                        powout_RA_line.append(float(nz[p]))
                        temp_ra_line = - float(nx[p])*180/np.pi
                        if temp_ra_line <= 0:
                            temp_ra_line += 360.
                        RA_line.append(temp_ra_line)

                else:
                    if abs(ny[p]*180/np.pi + 0.001 - dec) < 0.5*float(res):
                        powout_RA_line.append(float(nz[p]))
                        RA_line.append(180. - float(nx[p])*180/np.pi)

            #if ra_offset it needs to be restarted because it'll start at 180 not 0
            if args.ra_offset:
                powout_RA_line = [x for _,x in sorted(zip(RA_line,powout_RA_line))]
                RA_line = sorted(RA_line)
                #janky fix because there's two 360 values at the end
                RA_line = [0.] + RA_line[:-1]
                powout_RA_line = [powout_RA_line[-1]] + powout_RA_line[:-1]

            spline = UnivariateSpline(RA_line, powout_RA_line-np.max(powout_RA_line)/2., s=0)
            if len(spline.roots()) != 2:
                #print(spline.roots())
                #print(ra,dec)
                r1 = spline.roots()[-1]
                r2 = spline.roots()[0]
            else:
                r1 = spline.roots()[0]
                r2 = spline.roots()[1]

            diff = r2 - r1
            if diff > 180. and dec != -72.0:
                diff = r1 - (r2 -360)
                max_ra = r1 - (diff)/2.
            else:
                max_ra = r1 + (diff)/2.

            #max_ra = 180.-max_ra*180/np.pi
            if max_ra < 0.:
                max_ra += 360.
            if max_ra > 360.:
                max_ra -= 360.

            #if abs(max_ra - ra) > 180.:
            #    max_ra += 180.
            # I can't hunt down this error but this is what it needs
            if i == 69:
                max_ra -= 180.
            #if i in [0,69]:
            #    print(max_ra, dec)
            for c in range(len(colour_groups)):
                # Split the colour ranges into 5 ra ranges
                fudge_factor = 35.
                min_lim = 72.*c + fudge_factor
                max_lim = 72.*(c+1) + fudge_factor
                if max_lim >= 360.:
                    max_lim -= 360.
                    max_check = True
                else:
                    max_check = False

                # This was a temp bit of code to write the obs index on the plot
                #if args.ra_offset:
                #    if max_ra > 180:
                #        ra_text = -max_ra/180.*np.pi+2*np.pi
                #    else:
                #        ra_text = -max_ra/180.*np.pi
                #else:
                #    ra_text = -max_ra/180.*np.pi+np.pi
                #dec_text = dec/180.*np.pi
                #ax.text(ra_text, dec_text, str(i), fontsize=12, ha='center', va='center')

                # Check if this obs max power RA is in this colours group range
                if  (min_lim <= max_ra < max_lim) or \
                    ( max_check and ((min_lim <= max_ra < 360.   ) or \
                                     (0.      <= max_ra < max_lim)) ):
                        #print("Hit colour {}:".format(colour_groups[c]))
                        colors = ['0.5' for _ in range(50)]
                        colors[0] = colour_groups[c]

                        f = open(str(colour_groups[c]) + '_group_file.txt','a+')
                        f.write(str(max_ra) + '\t' + str(dec) + '\n')
                        f.close()
                        #plt.scatter(-max_ra/180*np.pi + np.pi, dec/180*np.pi, 1.5,\
                        #            lw=0, marker='o', color=colour_groups[c])


                        if args.shade:
                            if colour_groups[c] in args.shade:
                                #or ("blue" in args.shade and i in [0, 69]):
                                #sum powers for this colour to be shaded when plotting
                                for zi in range(len(nz)):
                                    if nz[zi] >= levels[0]:
                                        nz_shade_colour[colour_groups[c]][zi] = nz[zi]

                        # This is a temp feature that I'll delete to shade certain obs
                        if args.shade_temp and str(i) in args.shade_temp[1:]:
                            #print(i, args.shade_temp[1:])
                            for zi in range(len(nz)):
                                if nz[zi] >= levels[0]:
                                    #print(nz_shade_colour_temp[zi], nz[zi])
                                    nz_shade_colour_temp[zi] = nz[zi]

        # plot contours ---------------------------------------
        if args.contour:
            #print("plotting colour {}".format(colors[0]))
            plt.tricontour(nx, ny, nz, levels=[levels[0]], alpha = 0.6,
                           colors=colors,
                           linewidths=linewidths)


    # plot sens -------------------------------------------------------
    if args.sens:
        if args.overlap:
            for zi in range(len(nz)):
                if nz_sens_overlap[zi] < 0.5:
                    nz_sens_overlap[zi] = np.nan
            nz = 1.5*4.96/np.sqrt(nz_sens_overlap)
            #nz = nz_sens_overlap
        else:
            if args.incoh:
                nz = nz_sens * 11.3 #(sqrt128)
            else:
                nz = nz_sens

        with open('obs_plot_data.csv', 'w') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter=',')
            spamwriter.writerow(['RA','Dec','Sens mJy'])
            for ni in range(len(nx)):
                if args.ra_offset:
                    if RA[c] > 180:
                        x.append(-RA[c]/180.*np.pi+2*np.pi)
                    else:
                        x.append(-RA[c]/180.*np.pi)
                else:
                    x.append(-RA[c]/180.*np.pi +np.pi)
                y.append(Dec[c]/180.*np.pi)
                ra_temp = -math.degrees(nx[ni])
                if ra_temp < 0.:
                    ra_temp = ra_temp + 360.
                spamwriter.writerow([ra_temp, math.degrees(ny[ni]), nz[ni]])

        nx.shape = (len(map_dec_range),len(map_ra_range))
        ny.shape = (len(map_dec_range),len(map_ra_range))
        nz.shape = (len(map_dec_range),len(map_ra_range))
        if args.ra_offset:
            roll_by = len(map_ra_range)//2
            nx = np.roll(nx, roll_by)
            ny = np.roll(ny, roll_by)
            nz = np.roll(nz, roll_by)
        dec_limit_mask = ny > np.radians(63.3)
        nz[dec_limit_mask] = np.nan
        import matplotlib.colors as colors
        colour_map = 'plasma_r'
        if args.incoh:
            plt.pcolor(nx, ny, nz, cmap=colour_map, vmin=20, vmax=90)
        else:
            plt.pcolor(nx, ny, nz, cmap=colour_map, vmin=2., vmax=10.)
        plt.colorbar(spacing='uniform', shrink = 0.65, #ticks=[2., 10., 20., 30., 40., 50.],
                     label=r"Detection Sensitivity, 10$\sigma$ (mJy)")

    #Add extra plot layers ---------------------------------------

    #shades only the selected colout
    if (args.shade or args.shade_temp) and args.smart:
        for c in colour_groups:
            # Use the correct shade
            if args.shade:
                if c in args.shade:
                    nz = nz_shade_colour[c]
            elif args.shade_temp:
                if c == args.shade_temp[0]:
                    nz = nz_shade_colour_temp

            if args.shade:
                if c in args.shade or c == args.shade_temp[0]:
                    #choose lighter equivalent colour
                    if   c == 'red':
                        ecolour = 'lightcoral'
                    elif c == 'blue':
                        ecolour = 'skyblue'
                    else:
                        ecolour = c

                    if 'green' == c and args.ra_offset:
                        # add extra ra to fix shading issues
                        map_dec_range = np.arange(radians(-90),radians(91),radians(res))
                        map_ra_range = np.arange(radians(-220),radians(-180-res),radians(res))
                        #print(degrees(map_ra_range))
                        RA=[]; Dec=[]
                        for i in map_dec_range:
                            for j in map_ra_range:
                                Dec.append(i)
                                RA.append(j)
                        # Fill the nz
                        nz_temp = []
                        for dec in Dec:
                            if dec < radians(8):
                                nz_temp.append(1.0)
                            else:
                                nz_temp.append(0.)
                        # Append the new values
                        nx = np.append(np.array(RA),      nx)
                        ny = np.append(np.array(Dec),     ny)
                        nz = np.append(np.array(nz_temp), nz)
                        
                        print(c)
                        print("shapes nx: {} ny: {} nz: {}".format(nx.shape, ny.shape, nz.shape))
                        # add extra ra to fix shading issues
                        map_dec_range = np.arange(radians(-90),radians(91),radians(res))
                        map_ra_range = np.arange(radians(181),radians(221),radians(res))
                        #print(degrees(map_ra_range))
                        RA=[]; Dec=[]
                        for i in map_dec_range:
                            for j in map_ra_range:
                                Dec.append(i)
                                RA.append(j)
                        # Fill the nz
                        nz_temp = []
                        for dec in Dec:
                            if dec < radians(8):
                                nz_temp.append(1.0)
                            else:
                                nz_temp.append(0.)
                        # Append the new values
                        nx = np.append(nx, np.array(RA))
                        ny = np.append(ny, np.array(Dec))
                        nz = np.append(nz, np.array(nz_temp))
                        np.arange(radians(181+res),radians(221),radians(res))
                        print("shapes nx: {} ny: {} nz: {}".format(nx.shape, ny.shape, nz.shape))
                        print("len    nx: {} ny: {} nz: {}".format(len(nx), len(ny), len(nz)))
                    cs = plt.tricontour(nx.flatten(), ny.flatten(), nz.flatten(), levels=[levels[0]], alpha=0.0)
                    cs0 = cs.collections[0]
                    cspaths = cs0.get_paths()
                    for cspath in cspaths:
                        spch_0 = patches.PathPatch(cspath, facecolor=ecolour,
                                                edgecolor='gray',lw=0.5, alpha=0.45)
                        ax.add_patch(spch_0)


    #add lines of other surveys
    if args.lines:
        plt.plot(np.radians(np.array(map_ra_range)) - np.pi,
                 np.full(len(map_ra_range),np.radians(30.)),
                 'r',  label=r'MWA   ( 80 -    300 MHz)', zorder=130)
        plt.plot(np.array(map_ra_range)/180.*np.pi + -np.pi,
                 np.full(len(map_ra_range),0./180.*np.pi),
                 '--m',label=r'LOFAR ( 10 -    240 MHz)', zorder=130)
        plt.plot(np.array(map_ra_range)/180.*np.pi + -np.pi,
                 np.full(len(map_ra_range),-40./180.*np.pi),
                 '--g',label=r'GBT   (390 - 49,800 MHz)', zorder=130)
        """
        plt.plot(np.array(map_ra_range)/180.*np.pi + -np.pi,
                 np.full(len(map_ra_range),-55./180.*np.pi),
                 linestyle='--', color='orange',
                       label=r'GMRT  ( 50 -  1,500 MHz)', zorder=130)
        """

        plt.legend(bbox_to_anchor=(0.84, 0.85,0.21,0.2), loc='upper left', fontsize=6)



    if args.fill:
        import matplotlib.transforms as mtransforms
        trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
        map_ra_range = range(-80,481,res)
        ff = 30.
        ffa = 28.5
        ax.fill_between(np.array(map_ra_range)/180.*np.pi + -np.pi,
                        np.full(len(map_ra_range),np.radians((-100)/90.*ff+ffa)),
                        np.full(len(map_ra_range),np.radians((34.5)/90.*ff+ffa)),
                        facecolor='0.5', alpha=0.5, transform=trans)


    # Add pulsars to plot
    if args.pulsar_all:
        #add all pulsars on the antf catalogue
        ra_PCAT = []
        dec_PCAT = []
        pulsar_list = get_psrcat_ra_dec()
        for pulsar in pulsar_list:
            ra_temp, dec_temp = sex2deg(pulsar[1], pulsar[2])
            if args.ra_offset:
                if ra_temp > 180.:
                    ra_temp -= 180.
                else:
                    ra_temp += 180.
            if args.square:
                ra_PCAT.append(ra_temp)
                dec_PCAT.append(dec_temp)
            else:
                ra_PCAT.append(-(ra_temp-180.)/180.*np.pi)
                dec_PCAT.append(dec_temp/180.*np.pi)
        #print(min(ra_PCAT), max(ra_PCAT))
        ax.scatter(ra_PCAT, dec_PCAT, s=0.2, color ='b', zorder=90)

    if args.pulsar_detected:
        #add some pulsars
        ra_PCAT = []
        dec_PCAT = []
        from mwa_pulsar_client import client
        auth = (os.environ['MWA_PULSAR_DB_USER'],os.environ['MWA_PULSAR_DB_PASS'])
        pulsar_list_dict = client.pulsar_list('https://pulsar-cat.icrar.uwa.edu.au/', auth)
        pulsar_list = []
        for pulsar in pulsar_list_dict:
            pulsar_list.append(pulsar[u'name'])
        pulsar_pos_list = get_psrcat_ra_dec(pulsar_list=pulsar_list)
        for pulsar in pulsar_pos_list:
            ra_temp, dec_temp = sex2deg(pulsar[1], pulsar[2])
            if args.ra_offset:
                if ra_temp > 180:
                    ra_PCAT.append(-ra_temp/180.*np.pi+2*np.pi)
                else:
                    ra_PCAT.append(-ra_temp/180.*np.pi)
            else:
                ra_PCAT.append(-ra_temp/180.*np.pi+np.pi)
            dec_PCAT.append(dec_temp/180.*np.pi)
        ax.scatter(ra_PCAT, dec_PCAT, s=5, color ='r', zorder=100)

    if args.pulsar:
        #add some pulsars
        ra_PCAT = []
        dec_PCAT = []
        print("{} input pulsars".format(len(args.pulsar)))
        raw_pulsar_list = list(dict.fromkeys(args.pulsar))
        print("{} distinct pulsars".format(len(raw_pulsar_list)))
        pulsar_list = get_psrcat_ra_dec(pulsar_list=raw_pulsar_list)
        for pulsar in pulsar_list:
            ra_temp, dec_temp = sex2deg(pulsar[1], pulsar[2])
            if args.ra_offset:
                if ra_temp > 180:
                    ra_PCAT.append(-ra_temp/180.*np.pi+2*np.pi)
                else:
                    ra_PCAT.append(-ra_temp/180.*np.pi)
            else:
                ra_PCAT.append(-ra_temp/180.*np.pi+np.pi)
            dec_PCAT.append(dec_temp/180.*np.pi)
        ax.scatter(ra_PCAT, dec_PCAT, s=5, color ='purple', zorder=100)

    if args.pulsar_discovered:
        #add some pulsars
        ra_PCAT = []
        dec_PCAT = []
        print("{} input pulsars".format(len(args.pulsar)))
        raw_pulsar_list = list(dict.fromkeys(args.pulsar))
        print("{} distinct pulsars".format(len(raw_pulsar_list)))
        pulsar_list = [["J0036-1033", "00:36:14.58", "-10:33:16.40"]]
        for pulsar in pulsar_list:
            ra_temp, dec_temp = sex2deg(pulsar[1], pulsar[2])
            if args.ra_offset:
                if ra_temp > 180:
                    ra_PCAT.append(-ra_temp/180.*np.pi+2*np.pi)
                else:
                    ra_PCAT.append(-ra_temp/180.*np.pi)
            else:
                ra_PCAT.append(-ra_temp/180.*np.pi+np.pi)
            dec_PCAT.append(dec_temp/180.*np.pi)
        ax.scatter(ra_PCAT, dec_PCAT, s=10, color ='r', zorder=100)

    plt.xlabel("Right Ascension")
    plt.ylabel("Declination")

    #xtick_labels = ['0h','2h','4h','6h','8h','10h','12h','14h','16h','18h','20h','22h']
    if args.ra_offset:
        xtick_labels = ['10h', '8h', '6h', '4h', '2h', '0h', '22h', '20h', '18h', '16h', '14h']
        xticks = [150., 120., 90., 60., 30., 0., 330., 300., 270., 240., 210. ]
    else:
        xtick_labels = [ '22h', '20h', '18h', '16h', '14h','12h','10h', '8h', '6h', '4h', '2h']
        xticks = [330., 300., 270., 240., 210., 180., 150., 120., 90., 60., 30.]

    if args.square:
        plt.xticks(xticks, tuple(xtick_labels))
    else:
        ax.set_xticklabels(xtick_labels, zorder=150)
    print("plotting grid")
    plt.grid(True, color='gray', lw=0.5, linestyle='dotted')


    # Creates a plot name --------------------------
    plot_name = 'mwa_obs_n{}_res{}'.format(pointing_count, res)

    if args.sens:
        plot_name += '_sens'
    if args.contour:
        plot_name += '_contour'
    if args.lines:
        plot_name += '_minlines'
    if args.obsid_list:
        plot_name += '_obslist'
    if args.incoh:
        plot_name += '_incoh'
    if args.fill:
        plot_name += '_fill'
    if args.pulsar:
        plot_name += '_pulsar_n{}'.format(len(raw_pulsar_list))
    if args.pulsar_detected:
        plot_name += '_pulsar_detected'
    if args.pulsar_discovered:
        plot_name += '_pulsar_discovered'
    plot_type = args.plot_type
    #plt.title(plot_name)
    print("saving {}.{}".format(plot_name, plot_type))
    fig.savefig(plot_name + '.' + plot_type, format=plot_type, dpi=1000, bbox_inches='tight')
    #plt.show()


    if args.smart:
        #sort the output into the right order
        import glob
        from operator import itemgetter
        for g in glob.glob("./*group*"):
            with open(g) as f:
                lines = [line.split("\t") for line in f]
                lines = lines[1:]
                lines = sorted(lines, key=itemgetter(0))
            with open(g, 'w') as csvfile:
                spamwriter = csv.writer(csvfile, delimiter=',')
                spamwriter.writerow(['RA','Dec'])
                for l in lines:
                    spamwriter.writerow(["("+str(round(float(l[0]),1)),l[1][:-1]+")"])
