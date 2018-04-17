import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser(description="""
plot fwhm from grid.py or pabeam.py
""")
parser.add_argument('-o', '--obsid',type=str,help='Observation ID')
args=parser.parse_args()

theta,phi,beam = np.genfromtxt(str(args.obsid)+'_fwhm.0.dat',
                                    comments=("#","*"),skip_header=14,usecols=(0,1,8),unpack=True)
plt.plot(np.radians(phi),beam)
plt.savefig(str(args.obsid)+'_fwhm.0.png',bbox_inches='tight',dpi=900)
