



fwhm=0.30
loops=int(7.5/0.85/fwhm)
#loops=58
print loops
plot_08_deg =[]
plot_08_num =[]
plot_085_deg =[]
plot_085_num =[]
plot_09_deg =[]
plot_09_num =[]
plot_1_deg=[]
plot_1_num=[]

num = 1
for l in range(1,loops+1):
    num += 6*l
    plot_08_deg.append(l*fwhm*0.8)
    plot_085_deg.append(l*fwhm*0.85)
    plot_09_deg.append(l*fwhm*0.9)
    plot_1_deg.append(l*fwhm*0.95)
    plot_08_num.append(num)
    plot_085_num.append(num)
    plot_09_num.append(num)
    plot_1_num.append(num)

print plot_085_num[-1]

import pylab as plt  
import matplotlib
#matplotlib.use('Agg')
plt.ylabel("Number of pointings")
plt.xlabel("Degrees from beam centre")
plt.plot(plot_08_deg, plot_08_num, 'r', label='0.8')
plt.plot(plot_085_deg, plot_085_num, 'y', label='0.85')
plt.plot(plot_09_deg, plot_09_num, 'g', label='0.9')
plt.plot(plot_1_deg, plot_1_num, 'b', label='0.95')
plt.legend(loc='upper left')
plt.savefig('number_of_pointings_deg_3_pulsar_test_'+str(loops*fwhm)+'.png')

"""

signals=[["06:30:48.62","-28:34:43.14",2380.172,386.7],
         ["06:30:47.62","-28:34:43.14",2390.862,387.4],
         ["06:30:49.62","-28:34:43.14",2370.056,385.7],
         ["06:30:48.62","-28:34:42.14",2382.483,386.7],
         ["06:30:48.62","-28:34:44.14",2382.056,386.7],
         ["06:30:38.62","-28:34:43.14",2308.548,380.7],
         ["06:30:58.62","-28:34:43.14",2117.749,364.5],
         ["06:30:48.62","-28:34:33.14",2383.561,386.8],
         ["06:30:48.62","-28:34:53.14",2379.422,386.5],
         ["06:30:48.62","-28:34:38.14",2382.926,386.8],
         ["06:30:46.62","-28:34:38.14",2396.263,387.8],
         ["06:30:50.62","-28:34:38.14",2355.528,384.5],
         ["06:30:48.62","-28:34:36.14",2383.310,386.8],
         ["06:30:48.62","-28:34:40.14",2382.797,386.7],
         ["06:30:43.62","-28:34:43.14",2389.176,387.3],
         ["06:30:41.62","-28:34:43.14",2367.050,385.5],
         ["06:30:45.62","-28:34:43.14",2397.150,387.9],
         ["06:30:43.62","-28:34:41.14",2389.136,387.3],
         ["06:30:43.62","-28:34:45.14",2389.191,387.3],
         ["06:30:44.62","-28:34:43.14",2217.044,373.0],
         ["06:30:46.62","-28:34:43.14",2396.012,387.8],
         ["06:30:45.62","-28:33:43.14",2370.017,385.7],
         ["06:30:45.62","-28:35:43.14",2362.847,385.1]]
ra=[]
dec=[]
chi=[]
sigma=[]

from astropy.coordinates import SkyCoord
import astropy.units as u

for s in signals:
    coord = SkyCoord(s[0],s[1],unit=(u.hourangle,u.deg))
    ra.append(coord.ra.deg)
    dec.append(coord.dec.deg)
    chi.append(s[2])
    sigma.append(s[3])
    
print max(chi)

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(ra, dec, chi)

ax.set_xlabel('RA')
ax.set_ylabel('Dec')
ax.set_zlabel('Chi')

plt.savefig("dithering_powers.png")
plt.show()
"""
