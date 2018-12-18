#! /usr/bin/env python

import os
import sys
import math
import argparse
import urllib
import urllib2
import json
import subprocess
import numpy as np
import csv
from scipy.interpolate import UnivariateSpline
import ephem

#astropy
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table
from astropy.time import Time

#matplotlib
import matplotlib.pyplot as plt
import matplotlib.path as Path
import matplotlib.patches as patches
import matplotlib.tri as tri
import matplotlib.cm as cm
#plt.rc("text",usetex=True)

#vcstools/mwapy
from mwapy.pb import primary_beam
from mwapy import ephem_utils,metadata
import find_pulsar_in_obs as fpio
import mwa_metadb_utils as meta


if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="""
  A ploting script originaly written by Mengyao but editted by Nick to estimate the number of pointings required to cover the southern sky.
  southern_survey_sim.py -d 10 -f -r 3 -m 6 9 11 11 11 11 11 --semester_ra -a
  southern_survey_sim.py -d 10 -f  -r 1 -a -m 6 9 11 11 11 11 11 -c -o -l
  southern_survey_sim.py -d 10 -f -r 1 -a -c -o  --obsid_list 1088850560 1090249472
  """)
  parser.add_argument('-f','--fwhm',action='store_true',help='if this options is used the FWHM of each pointing is used. If it is not chosen the FWHM of a zenith pointing is used.')
  parser.add_argument('-d','--degree',type=float,help='The degrees overlap in RA of the observations')
  parser.add_argument('-r','--resolution',type=int,help='The resolution in degrees of the final plot (must be an integer). Default = 1', default=1)
  parser.add_argument('-a','--aitoff',action='store_true',help='Plots it in aitoff.')
  parser.add_argument('-s','--sens',action='store_true',help='Plots sensitivity')
  parser.add_argument('-o','--sens_overlap',action='store_true',help='Plots sensitivity that overlaps between observations.')
  parser.add_argument('-c','--colour',action='store_true',help='Plots sensitivity plots in colour instead of contour')
  parser.add_argument('-l','--lines',action='store_true',help='Includes the min decs of other telescopes in plots')
  parser.add_argument('--fill',action='store_true',help='Shades and area') 
  parser.add_argument('-m','--manual', nargs='+', type=int, help='Makes the pointing numbers manual, input them as 1 2 3 4 5 6 7')
  parser.add_argument('--semester', action='store_true', help='Changed the colours of the FWHM for each semester')
  parser.add_argument('--semester_ra', action='store_true', help='Similar to semester but uses a RA cut off (changes number per group)')
  parser.add_argument('-p','--plot_type',type=str,help='Determines the output plot type, Default="png".',default='png')
  parser.add_argument('--obsid_list',type=str,nargs='+',help='Instead of calculating which positions to use the script will use the input obsids. eg: "1088850560 1090249472"')
  parser.add_argument('--pulsar',type=str,nargs='+',help='A list of pulsar to mark on the plot')
  parser.add_argument('--all_obsids',action='store_true', help='Uses all VCS obsids on the MWA metadatabase.')
  parser.add_argument('--ra_offset',action='store_true', help='Offsets the RA by 180 so that 0h is in the centre')
  args=parser.parse_args()

  #Setting up some of the plots
  fig = plt.figure(figsize=(6, 4))
  plt.rc("font", size=8)
  if args.aitoff:
    fig.add_subplot(111)
    print "changing axis"
    ax = plt.axes(projection='mollweide')
  else:
    fig.add_subplot(111)
    ax = plt.axes()



  #levels = np.arange(0.25, 1., 0.05)
  colors= ['0.5' for _ in xrange(50)] ; colors[0]= 'blue'
  linewidths= [0.4 for _ in xrange(50)] ; linewidths[0]= 1.0
  
  #txtfile=open('/group/mwaops/xuemy/incoh_census/fold_code/get_obs_oblist_test.txt').readlines()
  
  #setting up the dec ranges
  dec_range = [-72., -55., -40.5, -26.7, -13., +1.6, +18.3] #Gleam pointings
  delays_range = [[0,0,0,0,6,6,6,6,12,12,12,12,18,18,18,18],\
                  [0,0,0,0,4,4,4,4,8,8,8,8,12,12,12,12],\
                  [0,0,0,0,2,2,2,2,4,4,4,4,6,6,6,6],\
                  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],\
                  [6,6,6,6,4,4,4,4,2,2,2,2,0,0,0,0],\
                  [12,12,12,12,8,8,8,8,4,4,4,4,0,0,0,0],\
                  [18,18,18,18,12,12,12,12,6,6,6,6,0,0,0,0]]
  
  print "Using GLEAM dec range: {}".format(dec_range)
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

  #setting up RA Dec ranges for power calculations
  res = args.resolution
  map_dec_range = range(-90,91,res)
  map_ra_range = range(0,361,res)
  RA=[] ; Dec=[]
  for i in map_dec_range:
      for j in map_ra_range:
          Dec.append(i)
          RA.append(j)

  #Working out the observations required -----------------------------------------------
  if args.all_obsids:
      #observations= find_pulsar_in_obs.find_obsids_meta_pages()
      observations = fpio.find_obsids_meta_pages(params={'mode':'VOLTAGE_START','cenchan':145})
      pointing_count = len(observations)
      print observations
  elif args.obsid_list:
      observations = args.obsid_list
      pointing_count = len(observations)
  else:
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
              
      manual_point_num = args.manual   
              
      observations = []
      ra_list =[]
      dec_list =[]
      FWHM = []
      FWHM_Dec = []
      pointing_count = 0
      for i in range(len(dec_range)):
        #calculating the FWHM at this dec
        ra_sex, deg_sex = fpio.deg2sex(start_ra, dec_range[i]).split()
        cord = [start_obsid, str(ra_sex), str(deg_sex), 1, delays_range[i],centrefreq, channels]
        #powout=get_beam_power(cord, zip(RA_FWHM_calc,Dec_FWHM_calc), dt=600)
        names_ra_dec = np.column_stack((['source']*len(RA_FWHM_calc), RA_FWHM_calc, Dec_FWHM_calc))
        powout = fpio.get_beam_power_over_time(cord, names_ra_dec, dt=600, degrees = True)
        powout_RA_line = []
        powout_Dec_line = []
        RA_line = []
        Dec_line = []
        for p in range(len(powout)):
          #print int(y[i]/np.pi*180.), int(dec) 
          if int(Dec_FWHM_calc[p]) == int(dec_range[i]):
              powout_RA_line.append(float(powout[p]))
              RA_line.append(float(RA_FWHM_calc[p]))
          if int (RA_FWHM_calc[p]) == int(start_ra):
              powout_Dec_line.append(float(powout[p]))
              Dec_line.append(float(Dec_FWHM_calc[p]))
        
        print "\nValues for Dec " + str(dec_range[i])
        #work out RA FWHM (not including the drift scan, 0sec observation)
        if args.fwhm:
            spline = UnivariateSpline(RA_line, powout_RA_line-np.max(powout_RA_line)/2., s=0)
        else:
            spline = UnivariateSpline(RA_line, powout_RA_line-np.full(len(powout_RA_line),0.5), s=0)
        try:
            r1, r2 = spline.roots()
        except ValueError:
            print "No FWHM for " + str(dec_range[i]) + " setting to 1000 to skip"
            FWHM.append(1000.)
            pointing_count -=1
        else:
            FWHM.append(float(r2-r1))
            print "FWHM along RA at dec "+ str(dec_range[i]) + ": " + str(FWHM[i])
        
        #work out Dec FWHM
        if args.fwhm:
            spline = UnivariateSpline(Dec_line, powout_Dec_line-np.max(powout_Dec_line)/2., s=0)
            print spline.roots(), max(powout_Dec_line)
            r1, r2 = spline.roots()
            FWHM_Dec.append(float(r2-r1))
            print "FWHM along Dec at dec "+ str(dec_range[i]) + ": " + str(FWHM_Dec[i])
        
        deg_move = total_angle = FWHM[i] - args.degree*math.cos(math.radians(dec_range[i])) + \
                    float(time)/3600.*15.*math.cos(math.radians(dec_range[i]))
        if args.manual:
            point_num_this_deg = manual_point_num[i]
        else:
            point_num_this_deg = int(360./deg_move) + 1
        print "Number for this dec: " +str(point_num_this_deg)
        deg_move = 360. / point_num_this_deg
        overlap_true = FWHM[i] + float(time)/3600.*15.*math.cos(math.radians(dec_range[i])) -\
                       360./point_num_this_deg
        print "True overlap this dec: " + str(overlap_true)
        
        # offset every second dec range by half a FWHM in RA
        if i % 2 == 0:
            ra_list.append(start_ra)
            observations.append(start_obsid)
        else:
            ra_list.append(start_ra+deg_move/math.cos(math.radians(dec_range[i])))
            observations.append(str(int(start_obsid)+int(deg_move*120)))
        dec_list.append(dec_range[i])
        pointing_count += 1
        for x in range(point_num_this_deg-1):
          temp_ra = ra_list[-1]+deg_move
          if temp_ra > 360.:
             temp_ra = temp_ra -360.
          ra_list.append(temp_ra)
          dec_list.append(dec_range[i])
          observations.append(str(int(observations[-1])+int(deg_move*240)))
          total_angle += deg_move
          pointing_count+=1
      
      #Sort by ra
      dec_list = [x for _,x in sorted(zip(ra_list,dec_list))]
      observations = [x for _,x in sorted(zip(ra_list,observations))]
      ra_list = sorted(ra_list)
      
  s_overlap_z = np.zeros(len(RA))
  sens_colour_z =[]
  max_ra_list = []
  RA_FWHM_atdec =[]
  with open('obs_meta.csv', 'rb') as csvfile: 
      spamreader = csv.reader(csvfile)
      next(spamreader, None)
      obsid_meta_file = []
      for row in spamreader:
          obsid_meta_file.append(row)
  for i, ob in enumerate(observations):
      if args.obsid_list or args.all_obsids:
          obs_foun_check = False
          for omi in range(len(obsid_meta_file)):
              if int(ob) == int(obsid_meta_file[omi][0]):
                  print "getting obs metadata from obs_meta.csv"
                  ob, ra, dec, time, delays,centrefreq, channels =\
                              obsid_meta_file[omi]
                  ob = int(ob)
                  time = int(time)
                  delays = map(int,delays[1:-1].split(","))
                  centrefreq = float(centrefreq)
                  channels = map(int,channels[1:-1].split(","))
                  obs_foun_check = True
          if not obs_foun_check:
              ob, ra, dec, time, delays,centrefreq, channels =\
                  meta.get_common_obs_metadata(ob)  
              
              with open('obs_meta.csv', 'a') as csvfile:
                  spamwriter = csv.writer(csvfile)
                  spamwriter.writerow([ob, ra, dec, time, delays,centrefreq, channels])
          cord = [ob, ra, dec, time, delays,centrefreq, channels]
      else:
        ra = ra_list[i]
        dec = dec_list[i]
        temp_dec_index = dec_range.index(dec)
        delays = delays_range[dec_range.index(dec)]
      
        cord = [ob, ra, dec, time, delays,centrefreq, channels]
  
      z=[] ; z_sens =[] ; x=[] ; y=[]
      
      
      

      #print max(Dec), min(RA), Dec.dtype
      time_intervals = 600 # seconds
      names_ra_dec = np.column_stack((['source']*len(RA), RA, Dec))
      powout = fpio.get_beam_power_over_time(cord, names_ra_dec, dt=time_intervals, degrees = True)
      #grab a line of beam power for the pointing declination
      if i == 0:
          print "len powers list: " + str(powout.shape)
      for c in range(len(RA)):
          s_overlap_z[c] += powout[c,0,0]*math.cos(Dec[c]/180.*np.pi)

          temppower=powout[c,0,0]
          temppower_sense=powout[c,0,0]
          for t in range(1,powout.shape[1]):
              power_ra = powout[c,t,0]
              temppower_sense += power_ra #average power kinds
              s_overlap_z[c] += powout[c,t,0]*math.cos(Dec[c]/180.*np.pi)
              if power_ra > temppower:
                  temppower = power_ra
          z_sens.append(temppower_sense)
          z.append(temppower)
          """
          for t in range(0,(time/time_intervals)):
              if i%(len(map_ra_range)) >= t:
                  power_ra = powout[i-t,0,0] #ra+t*15./3600 3deg
              else : 
                  power_ra = powout[i+len(map_ra_range)-t,0,0]
              if args.sens:
                  temppower += power_ra #average power kinds
              else:
                  if power_ra > temppower:
                      temppower = power_ra
              #print temppower, power_ra
          z.append(temppower)
          """ 
          if args.ra_offset:
              if RA[c] > 180:
                  x.append(-RA[c]/180.*np.pi+2*np.pi)
              else:
                  x.append(-RA[c]/180.*np.pi)
          else:
              x.append(-RA[c]/180.*np.pi +np.pi)
          y.append(Dec[c]/180.*np.pi)
      
      nx=np.array(x) ; ny=np.array(y); nz=np.array(z)
      if args.sens:
          nz_sense = []
          for zsi in range(len(z_sens)):
              if nz[zsi] < 0.01:
                  nz_sense.append(np.nan)
              else:
                  nz_sense.append(4.96/np.sqrt(z_sens[zsi]))
          nz_sense = np.array(nz_sense)
          sens_min = None
      if args.fwhm:
          if args.sens:
              #find the power at the centre (in Ra/time) of the FWHM
              max_index = np.argmax(nz_sense) #index of the centre of the beam
              if res > 2:
                  round_to = 1
              else:
                  round_to = 2
              max_dec = ny[max_index] + (FWHM_Dec[temp_dec_index]/2.)/180.*np.pi
              max_ra = nx[max_index]
              for m in range(len(nz)):
                  if nx[m] == max_ra and abs((ny[m] - max_dec)*180/np.pi) < (res/2.):
                      #print (ny[m] - max_dec)*180/np.pi
                      sens_min = nz_sense[m]
              levels = np.arange(sens_min, max(nz_sense), (max(nz_sense)-sens_min))#/10.)

              #plot in different colours
              colors= ['0.5' for _ in xrange(50)] ; colors[0]= 'g'
              ax.scatter(max_ra,max_dec, 1.5, lw=0, marker='o', color ='blue')
              plt.tricontour(nx, ny, nz_sense, levels=levels, alpha = 0.3,
                             colors=colors,
                             linewidths=linewidths)
              colors= ['0.5' for _ in xrange(50)] ; colors[0]= 'r'
          #else:
          levels = np.arange(0.5*max(nz), max(nz), 0.5/6.)
      else:
          levels = np.arange(0.5, 1., 0.05)
      
      if (args.semester or args.semester_ra) and i == 0:
          colour_groups = ['red','green','purple','darkorange','blue']
          for c in range(len(colour_groups)):
              f = open(str(colour_groups[c]) + '_group_file.txt','w')
              f.write('RA\tDec\n')
              f.close()
      
      if args.semester or args.semester_ra:
          #find middle ra for each pointing
          powout_RA_line = []
          RA_line = []
          for p in range(len(nz)):
            #print int(y[i]/np.pi*180.), int(dec) 
            
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
          r1, r2 = spline.roots()
      
          diff = r2 - r1
          if diff > 180.:
              diff = r1 - (r2 -360)
              #print r1,r2
              #print diff
              max_ra = r1 - (diff)/2.
          else:
              max_ra = r1 + (diff)/2.
          #max_ra = 180.-max_ra*180/np.pi 
          if max_ra < 0.:
              max_ra += 360.
          if max_ra > 360.:
              max_ra -= 360.
          #max_ra = max_ra*180/np.pi
          #print max_ra 
          #print str(ra)
      
      
      if args.semester:
          for c in range(len(colour_groups)):
              if 14*c <= i and i < 14*(c+1):
                  colors = [colour_groups[c] for _ in xrange(50)]
                  
                  f = open(str(colour_groups[c]) + '_group_file.txt','a+')
                  f.write(str(max_ra) + '\t' + str(dec) + '\n')
                  f.close()
                      
          alpha = 0.6
      elif args.semester_ra:
          for c in range(len(colour_groups)):
              fudge_factor = 35.
              min_lim = 72.*c + fudge_factor
              max_lim = 72.*(c+1) + fudge_factor
              if max_lim > 360.:
                  max_lim -= 360.
                  max_check = True
              else:
                  max_check = False
              
              if i == 0:
                  print min_lim,max_lim
              if c == (len(colour_groups)-1) and max_check:
                  if (min_lim < max_ra and max_ra < 360.) or (0. < max_ra and max_ra <= max_lim):
                      colors = ['0.5' for _ in xrange(50)]
                      colors[0] = colour_groups[c]

                      f = open(str(colour_groups[c]) + '_group_file.txt','a+')
                      f.write(str(max_ra) + '\t' + str(dec) + '\n')
                      f.close()
                      #plt.scatter(-max_ra/180*np.pi + np.pi, dec/180*np.pi, 1.5,\
                      #            lw=0, marker='o', color=colour_groups[c])

              if min_lim < max_ra and max_ra <= max_lim:
                  colors = ['0.5' for _ in xrange(50)]
                  colors[0] = colour_groups[c]

                  f = open(str(colour_groups[c]) + '_group_file.txt','a+')
                  f.write(str(max_ra) + '\t' + str(dec) + '\n')
                  f.close()
                  #plt.scatter(-max_ra/180*np.pi + np.pi, dec/180*np.pi, 1.5,\
                  #            lw=0, marker='o', color=colour_groups[c])
           
          
          alpha = 0.5

      else:
          alpha = 0.3
      
      if not args.sens_overlap:
          if args.colour:
              if i == 0:
                  sens_colour_z = nz_sense
              else:
                  for zi, zs in enumerate(nz_sense):
                      if math.isnan(sens_colour_z[zi]):
                          sens_colour_z[zi] = zs
                      elif sens_colour_z[zi] > zs: #TODO change back after sensitivity plot over
                          #append if larger
                          sens_colour_z[zi] = zs
            
          else:
              plt.tricontour(nx, ny, nz, levels=levels, alpha = 0.6, 
                     colors=colors,
                     linewidths=linewidths)
          
      
      #shades only the blue ones
      if colors[0] == 'blue':
          print "check"
          cs = plt.tricontour(nx, ny, nz, levels=levels[0],alpha=0.0)
          cs0 = cs.collections[0]
          cspaths = cs0.get_paths()
          spch_0 = patches.PathPatch(cspaths[0], facecolor='skyblue', 
                                     edgecolor='gray',lw=0.5, alpha=0.55)
          ax.add_patch(spch_0)

  if args.sens_overlap:
      print "making np arrays"
      nz=np.sqrt(np.array(s_overlap_z))
      #nz = nz/max(nz)
      if args.colour:
        colour_map = 'plasma_r'
        nx.shape = (len(map_dec_range),len(map_ra_range))
        ny.shape = (len(map_dec_range),len(map_ra_range))
        nz.shape = (len(map_dec_range),len(map_ra_range))
        print "colour plotting"
        plt.pcolor(nx, ny, nz, cmap=colour_map)
        plt.colorbar(spacing='uniform', shrink = 0.65)
      else:
        levels = np.arange(0.5*max(nz), max(nz), max(nz)/20.)
        print "plotting"
        plt.tricontour(nx, ny, nz, levels=levels, alpha = 0.3,
                                   colors=colors,
                                   linewidths=linewidths)
   
  if args.sens and args.colour:
      nz=sens_colour_z
      colour_map = 'plasma_r'
      nx.shape = (len(map_dec_range),len(map_ra_range))
      ny.shape = (len(map_dec_range),len(map_ra_range))
      nz.shape = (len(map_dec_range),len(map_ra_range))
      dec_limit_mask = ny > np.radians(30)
      nz[dec_limit_mask] = np.nan
      plt.pcolor(nx, ny, nz, cmap=colour_map)
      plt.colorbar(spacing='uniform', shrink = 0.65, label=r"Detection Sensitivity, 10$\sigma$ (mJy)")
        
  
  if args.semester or args.semester_ra:
      #sort the output into the right order
      import glob
      from operator import itemgetter
      import csv
      for g in glob.glob("./*group*"):
          with open(g) as f:
            lines = [line.split("\t") for line in f]
            lines = lines[1:]
            lines = sorted(lines, key=itemgetter(0))
          with open(g, 'wb') as csvfile:
             spamwriter = csv.writer(csvfile, delimiter=',')
             spamwriter.writerow(['RA','Dec'])
             for l in lines:
                 spamwriter.writerow(["("+str(round(float(l[0]),1)),l[1][:-1]+")"])

  
  plt.xlabel("Right Ascension")
  plt.ylabel("Declination")
      
  #xtick_labels = ['0h','2h','4h','6h','8h','10h','12h','14h','16h','18h','20h','22h']
  if args.ra_offset:
      xtick_labels = ['10h', '8h', '6h', '4h', '2h', '0h', '22h', '20h', '18h', '16h', '14h']
  else:
      xtick_labels = [ '22h', '20h', '18h', '16h', '14h','12h','10h', '8h', '6h', '4h', '2h']

  ax.set_xticklabels(xtick_labels) 
  print "plotting grid"
  plt.grid(True, color='gray', lw=0.5, linestyle='dotted')
  

  #add lines of other surveys
  if args.lines:
    """
    plt.plot(np.array(map_ra_range)/180.*np.pi + -np.pi, np.full(len(map_ra_range),0./180.*np.pi),\
          'r',label='LOFAR limit')
    plt.plot(np.array(map_ra_range)/180.*np.pi + -np.pi, np.full(len(map_ra_range),-40./180.*np.pi),\
          '--g',label='GBT limit')
    plt.plot(np.array(map_ra_range)/180.*np.pi + -np.pi, np.full(len(map_ra_range),-55./180.*np.pi),\
          '--b',label='GMRT limit')
    """
    plt.plot(np.radians(np.array(map_ra_range)) - np.pi, 
             np.full(len(map_ra_range),np.radians(-15.)),
             '--b',label=r'FAST $\delta_{min}$')
    plt.plot(np.radians(np.array(map_ra_range)) - np.pi, 
             np.full(len(map_ra_range),np.radians(30.)),
             '--r',label=r'MWA $\delta_{max}$')
  if args.fill:
      import matplotlib.transforms as mtransforms
      trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
      map_ra_range = range(-20,381,res)
      ff = 30.
      ffa = 28.5
      ax.fill_between(np.array(map_ra_range)/180.*np.pi + -np.pi,
                      np.full(len(map_ra_range),np.radians((-16.5)/90.*ff+ffa)),
                      np.full(len(map_ra_range),np.radians((34.5)/90.*ff+ffa)),
                      facecolor='0.5', alpha=0.5, transform=trans)

  handles, labels = ax.get_legend_handles_labels()
  plt.legend(bbox_to_anchor=(0.84, 0.85,0.21,0.2), loc=3,numpoints=1,
             ncol=1, mode="expand", borderaxespad=0., fontsize=6)
  
  if args.pulsar:
      from find_pulsar_in_obs import get_psrcat_ra_dec, sex2deg
      #add some pulsars
      ra_PCAT = []
      dec_PCAT = []
      pulsar_list = get_psrcat_ra_dec(pulsar_list = args.pulsar)
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
      ax.scatter(ra_PCAT, dec_PCAT, s=15, color ='r', zorder=100)
  
  if args.aitoff:
      plot_name = 'tile_beam_t'+str(time)+'s_o'+str(args.degree)+'deg_res' + str(res) + '_n'+str(pointing_count)+'_aitoff'
  else:
      plot_name = 'tile_beam_t'+str(time)+'s_o'+str(args.degree)+ 'deg_res' + str(res) + '_n'+str(pointing_count)
  
  if args.sens:
      plot_name += '_sens'
  
  if args.sens_overlap:
      plot_name += '_overlap'
      if args.colour:
          plot_name += '_colour_' + str(colour_map)
      else:
          plot_name += '_contour'
  if args.lines:
      plot_name += '_minlines'
  
  if args.semester:
      plot_name += '_semester'
  
  if args.obsid_list:
      plot_name += '_obslist'
  if args.manual:
      plot_name += '_manual'
      for m in args.manual:
          plot_name += str(m) + '-'
      plot_name = plot_name[:-1]
  
  if args.fwhm:
    plot_name += '_ownFWHM'
  else:
    plot_name +='_zenithFWHM'
  
  plot_type = args.plot_type
  #plt.title(plot_name)
  print "saving {}.{}".format(plot_name, plot_type)
  fig.savefig(plot_name + '.' + plot_type, format=plot_type, dpi=1000)
  #plt.show()
  plt.close

