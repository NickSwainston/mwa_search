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
from scipy.interpolate import UnivariateSpline
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table
from astropy.time import Time
from mwapy.pb import primary_beam
import ephem
from mwapy import ephem_utils,metadata
import matplotlib.pyplot as plt
import matplotlib.path as Path
import matplotlib.patches as patches
import matplotlib.tri as tri
import matplotlib.cm as cm
#from mpl_toolkits.basemap import Basemap

def get_beam_power(obsid_data,
                   sources,
                   dt=296,
                   centeronly=True,
                   verbose=False,
                   min_power=0.6):
    """
    obsid_data = [obsid,ra, dec, time, delays,centrefreq, channels]
    sources=[names,coord1,coord2] #astropy table coloumn names

    Calulates the power (gain at coordinate/gain at zenith) for each source and if it is above
    the min_power then it outputs it to the text file.

    """
    #print "Calculating beam power"
    obsid,ra, dec, time, delays,centrefreq, channels = obsid_data
    
    starttimes=np.arange(0,time,dt)
    #starttimes=np.arange(0,time,time)
    stoptimes=starttimes+dt
    stoptimes[stoptimes>time]=time
    Ntimes=len(starttimes)
    midtimes=float(obsid)+0.5*(starttimes+stoptimes)

    mwa = ephem_utils.Obs[ephem_utils.obscode['MWA']]
    # determine the LST
    observer = ephem.Observer()
    # make sure no refraction is included
    observer.pressure = 0
    observer.long = mwa.long / ephem_utils.DEG_IN_RADIAN
    observer.lat = mwa.lat / ephem_utils.DEG_IN_RADIAN
    observer.elevation = mwa.elev

    if not centeronly:
        PowersX=np.zeros((len(sources),
                             Ntimes,
                             len(channels)))
        PowersY=np.zeros((len(sources),
                             Ntimes,
                             len(channels)))
        # in Hz
        frequencies=np.array(channels)*1.28e6
    else:
        PowersX=np.zeros((len(sources),
                             Ntimes,1))
        PowersY=np.zeros((len(sources),
                             Ntimes,1))
        frequencies=[centrefreq]

    RAs=np.array([x[0] for x in sources])
    Decs=np.array([x[1] for x in sources])
    if len(RAs)==0:
        sys.stderr.write('Must supply >=1 source positions\n')
        return None
    if not len(RAs)==len(Decs):
        sys.stderr.write('Must supply equal numbers of RAs and Decs\n')
        return None
    for itime in xrange(Ntimes):
        obstime = Time(midtimes[itime],format='gps',scale='utc')
        observer.date = obstime.datetime.strftime('%Y/%m/%d %H:%M:%S')
        LST_hours = observer.sidereal_time() * ephem_utils.HRS_IN_RADIAN

        HAs = -RAs + LST_hours * 15
        Azs, Alts = ephem_utils.eq2horz(HAs, Decs, mwa.lat)
        # go from altitude to zenith angle
        theta=np.radians(90-Alts)
        phi=np.radians(Azs)
        #az, za = np.meshgrid(sorted(set(phi)), sorted(set(theta)))

        for ifreq in xrange(len(frequencies)):
            rX,rY=primary_beam.MWA_Tile_analytic(theta, phi,
                                                 freq=frequencies[ifreq], delays=delays,
                                                 zenithnorm=True,
                                                 power=True)#pixels_per_deg=1)
            PowersX[:,itime,ifreq]=rX
            PowersY[:,itime,ifreq]=rY

    #Power [#sources, #times, #frequencies]
    Powers=0.5*(PowersX+PowersY)
                 
    return Powers

# Append the service name to this base URL, eg 'con', 'obs', etc.
BASEURL = 'http://mwa-metadata01.pawsey.org.au/metadata/'

def getmeta(service='obs', params=None):
  """Given a JSON web service ('obs', find, or 'con') and a set of parameters as
     a Python dictionary, return a Python dictionary containing the result.
  """
  if params:
    data = urllib.urlencode(params)  # Turn the dictionary into a string with encoded 'name=value' pairs
  else:
    data = ''
  print data
  #             Validate the service name
  if service.strip().lower() in ['obs', 'find', 'con']:
    service = service.strip().lower()
  else:
    print "invalid service name: %s" % service
    return
  #             Get the data
  try:
    result = json.load(urllib2.urlopen(BASEURL + service + '?' + data))
  except urllib2.HTTPError as error:
    print "HTTP error from server: code=%d, response:\n %s" % (error.code, error.read())
    return
  except urllib2.URLError as error:
    print "URL or network error: %s" % error.reason
    return
  #            Return the result dictionary
  return result


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
  parser.add_argument('-m','--manual', nargs='+', type=int, help='Makes the pointing numbers manual, input them as 1 2 3 4 5 6 7')
  parser.add_argument('--semester', action='store_true', help='Changed the colours of the FWHM for each semester')
  parser.add_argument('--semester_ra', action='store_true', help='Similar to semester but uses a RA cut off (changes number per group)')
  parser.add_argument('-p','--plot_type',type=str,help='Determines the output plot type, Default="png".',default='png')
  parser.add_argument('--obsid_list',type=str,nargs='+',help='Instead of calculating which positions to use the script will use the input obsids. eg: "1088850560 1090249472"')
  args=parser.parse_args()

  #Setting up some of the plots
  fig = plt.figure(figsize=(6, 4))
  plt.rc("font", size=8)
  if args.aitoff:
    fig.add_subplot(111)
    ax = plt.axes(projection='mollweide')
  else:
    fig.add_subplot(111)
    ax = plt.axes()


  #levels = np.arange(0.25, 1., 0.05)
  colors= ['0.5' for _ in xrange(50)] ; colors[0]= 'r'
  linewidths= [0.4 for _ in xrange(50)] ; linewidths[0]= 0.7
  
  #txtfile=open('/group/mwaops/xuemy/incoh_census/fold_code/get_obs_oblist_test.txt').readlines()
  
  dec_range = [-72., -55., -40.5, -26.7, -13., +1.6, +18.3] #Gleam pointings
  delays_range = [[0,0,0,0,6,6,6,6,12,12,12,12,18,18,18,18],\
                  [0,0,0,0,4,4,4,4,8,8,8,8,12,12,12,12],\
                  [0,0,0,0,2,2,2,2,4,4,4,4,6,6,6,6],\
                  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],\
                  [6,6,6,6,4,4,4,4,2,2,2,2,0,0,0,0],\
                  [12,12,12,12,8,8,8,8,4,4,4,4,0,0,0,0],\
                  [18,18,18,18,12,12,12,12,6,6,6,6,0,0,0,0]]
  
  #setting up the dec ranges
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
  """
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
  print len(map_dec_range),len(map_ra_range)
  RA=[] ; Dec=[]
  for i in map_dec_range:
      for j in map_ra_range:
          Dec.append(i)
          RA.append(j)

  
  if not args.obsid_list:
      #setting up some metadata requirements
      time = 4800 #one hour 
      channels = range(107,131)
      minfreq = float(min(channels))
      maxfreq = float(max(channels))
      centrefreq = 1.28e6 * (minfreq + (maxfreq-minfreq)/2) #in MHz

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
        cord = [start_obsid, start_ra, dec_range[i], 1, delays_range[i],centrefreq, channels]
        powout=get_beam_power(cord, zip(RA_FWHM_calc,Dec_FWHM_calc), dt=600)

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
      
      dec_list = [x for _,x in sorted(zip(ra_list,dec_list))]
      observations = [x for _,x in sorted(zip(ra_list,observations))]
      ra_list = sorted(ra_list)
  else:
      observations = args.obsid_list
      pointing_count = len(observations)
      
  s_overlap_z = np.zeros(len(RA))
  s_overlap_x =[]
  s_overlap_y =[]
  max_ra_list = []
  RA_FWHM_atdec =[]
  for i in range(len(observations)):
      ob = observations[i]
      if args.obsid_list:
        print "Obtaining metadata from http://mwa-metadata01.pawsey.org.au/metadata/ for OBS ID: "
        beam_meta_data = getmeta(service='obs', params={'obs_id':ob})
        ra = beam_meta_data[u'metadata'][u'ra_pointing']
        dec = beam_meta_data[u'metadata'][u'dec_pointing']
        time = beam_meta_data[u'stoptime'] - beam_meta_data[u'starttime'] #gps time
        delays = beam_meta_data[u'rfstreams'][u'0'][u'xdelays']
             
        minfreq = float(min(beam_meta_data[u'rfstreams'][u"0"][u'frequencies']))
        maxfreq = float(max(beam_meta_data[u'rfstreams'][u"0"][u'frequencies']))
        centrefreq = 1.28e6 * (minfreq + (maxfreq-minfreq)/2) #in MHz
        channels = beam_meta_data[u'rfstreams'][u"0"][u'frequencies']
      else:
        ra = ra_list[i]
        dec = dec_list[i]
        temp_dec_index = dec_range.index(dec)
        delays = delays_range[dec_range.index(dec)]
      
      cord = [ob, ra, dec, time, delays,centrefreq, channels]
  
      z=[] ; z_sens =[] ; x=[] ; y=[]
      
      
      

      #print max(Dec), min(RA), Dec.dtype
      time_intervals = 600 # seconds
      powout=get_beam_power(cord, zip(RA,Dec), dt=time_intervals)
      
      #grab a line of beam power for the pointing declination
      if i == 0:
          print "len powers list: " + str(powout.shape)
      for c in range(len(RA)):
          if i == 0:
              s_overlap_x.append(-RA[c]/180.*np.pi+np.pi)
              s_overlap_y.append(Dec[c]/180.*np.pi)
          
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
          
          if RA[i] > 180:
              x.append(-RA[i]/180.*np.pi+2*np.pi)
              #x.append(-RA[i]+360.)
          else:
              x.append(-RA[i]/180.*np.pi)
              #x.append(-RA[i])
          """
          x.append(-RA[c]/180.*np.pi+np.pi)
          y.append(Dec[c]/180.*np.pi)
          #y.append(Dec[i])
      nx=np.array(x) ; ny=np.array(y); nz=np.array(z)
      if args.sens:
          nz_sense = np.sqrt(np.array(z_sens))
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
      
      #find middle ra for each pointing
      powout_RA_line = []
      RA_line = []
      for p in range(len(nz)):
        #print int(y[i]/np.pi*180.), int(dec) 
        
        if abs(ny[p]*180/np.pi + 0.001 - dec) < 0.5*float(res):
          powout_RA_line.append(float(nz[p]))
          RA_line.append(180. - float(nx[p])*180/np.pi)
      
      #print RA_line
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
        plt.tricontour(nx, ny, nz, levels=levels, alpha = alpha, 
                     colors=colors,
                     linewidths=linewidths)
      """ 
      cs = plt.tricontour(nx, ny, nz, levels=levels[0],alpha=0)
      cs0 = cs.collections[0]
      cspaths = cs0.get_paths()
      spch_0 = patches.PathPatch(cspaths[0], facecolor='k', edgecolor='gray',lw=0.5, alpha=0.1)
      ax.add_patch(spch_0)
      """
  
  if args.sens_overlap:
      nx=np.array(s_overlap_x) ; ny=np.array(s_overlap_y)
      nz=np.sqrt(np.array(s_overlap_z))
      nz = nz/max(nz)
      if args.colour:
        colour_map = 'plasma'
        nx.shape = (len(map_dec_range),len(map_ra_range))
        ny.shape = (len(map_dec_range),len(map_ra_range))
        nz.shape = (len(map_dec_range),len(map_ra_range))
        plt.pcolor(nx, ny, nz, cmap=colour_map)
        plt.colorbar(spacing='uniform', shrink = 0.65)
      else:
        levels = np.arange(0.5*max(nz), max(nz), max(nz)/20.)
        plt.tricontour(nx, ny, nz, levels=levels, alpha = 0.3,
                                   colors=colors,
                                   linewidths=linewidths)
   
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
            print lines
          with open(g, 'wb') as csvfile:
             spamwriter = csv.writer(csvfile, delimiter=',')
             spamwriter.writerow(['RA','Dec'])
             for l in lines:
                 spamwriter.writerow(["("+str(round(float(l[0]),1)),l[1][:-1]+")"])


  #xtick_labels = ['0h','2h','4h','6h','8h','10h','12h','14h','16h','18h','20h','22h']
  xtick_labels = [ '22h', '20h', '18h', '16h', '14h','12h','10h', '8h', '6h', '4h', '2h']
  ax.set_xticklabels(xtick_labels) 
  plt.grid(True, color='gray', lw=0.5, linestyle='dotted')
  
  #p1=ax.scatter(ra_PCAT_N, dec_PCAT_N, 1.5, lw=0, marker='o', color ='gray', label="Known pulsars beyond MWA Dec limitation")
  #p1=ax.scatter(ra_PCAT, dec_PCAT, 1.5, lw=0, marker='o', color ='blue', label="Known pulsars MWA could reach")
  

  #add lines of other surveys
  if args.lines:
    plt.plot(np.array(map_ra_range)/180.*np.pi + -np.pi, np.full(len(map_ra_range),0./180.*np.pi),\
          'r',label='LOFAR limit')
    plt.plot(np.array(map_ra_range)/180.*np.pi + -np.pi, np.full(len(map_ra_range),-40./180.*np.pi),\
          '--g',label='GBT limit')
    plt.plot(np.array(map_ra_range)/180.*np.pi + -np.pi, np.full(len(map_ra_range),-55./180.*np.pi),\
          '--b',label='GMRT limit')

  handles, labels = ax.get_legend_handles_labels()
  plt.legend(bbox_to_anchor=(0.84, 0.85,0.21,0.2), loc=3,numpoints=1,
             ncol=1, mode="expand", borderaxespad=0., fontsize=6)
  
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
  plt.savefig(plot_name + '.' + plot_type, format=plot_type, dpi=1000)
  plt.show()

