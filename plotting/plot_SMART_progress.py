#! /usr/bin/env python3

import argparse
import numpy as np
import math

#matplotlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches

#astropy
from astropy.coordinates import SkyCoord
from astropy import units as u

def sex2deg(ra, dec):
    """
    Convert sexagesimal coordinates to degrees.
    sex2deg( ra, dec)
    Args:
        ra: the right ascension in HH:MM:SS
        dec: the declination in DD:MM:SS
    """
    c = SkyCoord( ra, dec, frame='icrs', unit=(u.hourangle,u.deg))

    # return RA and DEC in degrees in degrees
    return [c.ra.deg, c.dec.deg]


def get_psrcat_ra_dec(pulsar_list=None, max_dm=1000., include_dm=False, query=None):
    """
    Uses PSRCAT to return a list of pulsar names, ras and decs. Not corrected for proper motion.
    Removes pulsars without any RA or DEC recorded
    If no pulsar_list given then returns all pulsar on the catalogue
    If include_dm is True then also ouput DM
    get_psrcat_ra_dec(pulsar_list = None)
    Args:
        pulsar_list: A space list of pulsar names eg: [J0534+2200, J0538+2817].
               (default: uses all pulsars)
    return [[Jname, RAJ, DecJ]]
    """
    import psrqpy

    #params = ['JNAME', 'RAJ', 'DECJ', 'DM']
    if query is None:
        query = psrqpy.QueryATNF(params = ['PSRJ', 'RAJ', 'DECJ', 'DM'], psrs=pulsar_list).pandas

    pulsar_ra_dec = []
    for i, _ in enumerate(query["PSRJ"]):
        # Only record if under the max_dm
        dm = query["DM"][i]
        if not math.isnan(dm):
            if float(dm) < max_dm:
                if include_dm:
                    pulsar_ra_dec.append([query["PSRJ"][i], query["RAJ"][i], query["DECJ"][i], dm])
                else:
                    pulsar_ra_dec.append([query["PSRJ"][i], query["RAJ"][i], query["DECJ"][i]])


    return pulsar_ra_dec


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    A plotting script that can be used to plot the SMART surveys progress.
    python plot_SMART_progress.py --ra_offset --shade 1221832280 1222697776 1255197408 --shade_light 1224252736
    """)
    parser.add_argument('--shade', type=int, nargs='+',
                           help='The space seperated obsIDs to be shaded.')
    parser.add_argument('--shade_light', type=int, nargs='*',
                           help='The space seperated obsIDs to be shaded lightly. Used for observations in progress.')
    parser.add_argument('--pulsar', type=str, nargs='+',
                           help='A list of pulsar J names to mark on the plot.')
    parser.add_argument('--pulsar_cand', type=str, nargs='+',
                           help='A list of pulsar cands in the format "HH:MM:SS +DD:MM:SS HH:MM:SS +DD:MM:SS".')
    parser.add_argument('-r', '--resolution', type=int, default=1,
                            help='The resolution in degrees of the final plot (must be an integer). Default = 1')
    parser.add_argument('-p', '--plot_type', type=str,
                            help='Determines the output plot type, Default="png".',default='png')
    parser.add_argument('--ra_offset', action='store_true',
                            help='Offsets the RA by 180 so that 0h is in the centre')
    args=parser.parse_args()

    #Setting up some of the plots
    fig = plt.figure(figsize=(6, 4))
    plt.rc("font", size=8)
    fig.add_subplot(111)
    ax = plt.axes(projection='mollweide')

    #levels = np.arange(0.25, 1., 0.05)
    colors= ['0.5' for _ in range(50)] ; colors[0]= 'blue'
    linewidths= [0.4 for _ in range(50)] ; linewidths[0]= 1.0
    alpha = 0.5

    smart_colours = {'B': {'light': '', 'dark': 'blue'},
                     'R': {'light': '', 'dark': 'red'},
                     'G': {'light': '', 'dark': 'green'},
                     'P': {'light': '', 'dark': 'purple'},
                     'O': {'light': '', 'dark': 'darkorange'}}

    #[id, smart_name, obsid, ra, dec]
    SMART_metadata = [[0,  "B01", 1221399680, 330.4, -55.0],
                      [1,  "B11", 1224859816, 26.8, -40.5],
                      [2,  "B13", 1225462936, 26.7, -13.0],
                      [3,  "B12", 1225118240, 26.7, 18.3],
                      [4,  "B08", 1255444104, 10.3, 1.6],
                      [5,  "B07", 1226062160, 10.3, -26.7],
                      [6,  "R04", 1253991112, 59.6, -40.5],
                      [7,  "R06", 1255197408, 59.5, -13.0],
                      [8,  "R05", 1254594264, 59.5, 18.3],
                      [9,  "R01", 1252177744, 43.1, 1.6],
                      [10, "B09", 1224252736, 10.3, -55.0],
                      [11, "R02", 1252780888, 43.2, -26.7],
                      [12, "R07", 1255803168, 70.7, -72.0],
                      [13, "R11", 1258221008, 92.4, -40.5],
                      [14, "R13", 1259427304, 92.4, -13.0],
                      [15, "R12", 1259685792, 92.3, 18.3],
                      [16, "R08", 1256407632, 75.9, 1.6],
                      [17, "R09", 1257010784, 76.0, -26.7],
                      [18, "R03", 1253471952, 50.5, -55.0],
                      [19, "G03", 1265983624, 125.2, -40.5],
                      [20, "G05", 1266155952, 125.2, -13.0],
                      [21, "G04", 1265725128, 125.1, 18.3],
                      [22, "G01", 1260638120, 108.8, 1.6],
                      [23, "G02", 1261241272, 108.8, -26.7],
                      [24, "G07", 1266932744, 130.8, -72.0],
                      [25, "R10", 1257617424, 90.6, -55.0],
                      [26, "G10", 1266680784, 158.0, -40.5],
                      [27, "G12", 1267283936, 158.0, -13.0],
                      [28, "G11", 1267111608, 158.0, 18.3],
                      [29, "G08", 1264867416, 141.6, 1.6],
                      [30, "G09", 1265470568, 141.6, -26.7],
                      [31, "G06", 1266329600, 130.7, -55.0],
                      [32, "P03", 32, 189.4, -72.0],
                      [33, "P04", 33, 189.7, -40.5],
                      [34, "P01", 34, 189.8, -13.0],
                      [35, "P02", 35, 189.8, 18.3],
                      [36, "G14", 1268063336, 174.4, 1.6],
                      [37, "G15", 1268321832, 174.4, -26.7],
                      [38, "G13", 1267459328, 170.8, -55.0],
                      [39, "P10", 39, 222.6, -40.5],
                      [40, "P08", 40, 222.6, -13.0],
                      [41, "P09", 41, 222.6, 18.3],
                      [42, "P06", 42, 206.2, 1.6],
                      [43, "P05", 43, 206.2, -26.7],
                      [44, "P13", 44, 249.8, -72.0],
                      [45, "O03", 45, 255.4, -40.5],
                      [46, "O01", 46, 255.4, -13.0],
                      [47, "O02", 47, 255.5,  18.3],
                      [48, "P12", 48, 239.0, 1.6],
                      [49, "P11", 49, 239.0, -26.7],
                      [50, "P07", 50, 209.8, -55.0],
                      [51, "O08", 51, 288.2, -40.5],
                      [52, "O06", 52, 288.2, -13.0],
                      [53, "O07", 53, 288.3, 18.3],
                      [54, "O05", 54, 271.8, 1.6],
                      [55, "O04", 55, 271.8, -26.7],
                      [56, "P14", 56, 249.9, -55.0],
                      [57, "O12", 57, 310.0, -72.0],
                      [58, "O15", 58, 321.0, -40.5],
                      [59, "O13", 59, 321.1, -13.0],
                      [60, "O14", 60, 321.1, 18.3],
                      [61, "O11", 61, 304.7, 1.6],
                      [62, "O10", 62, 304.6, -26.7],
                      [63, "O09", 63, 290.0, -55.0],
                      [64, "B06", 1225713560, 353.9, -40.5],
                      [65, "B04", 1222697776, 353.9, -13.0],
                      [66, "B05", 1223042480, 353.9, 18.3],
                      [67, "B02", 1221832280, 337.5, 1.6],
                      [68, "B03", 1222435400, 337.5, -26.7],
                      [69, "B10", 1227009976, 10.6, -72.0]]


    map_dec_range = range(-90, 91, args.resolution)
    map_ra_range = range(0, 361, args.resolution)
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
    smart_nz = []

    with open("SMART_obs_data.npy", 'rb') as f:
        for i in range(70):
            smart_nz.append(np.load(f))

    for sobs in SMART_metadata:
        print(sobs)
        sid, sname, sobsid, sra, sdec = sobs
        nz = np.array(smart_nz[sid])
        
        # plot the contour
        plt.tricontour(nx, ny, nz, levels=[max(nz)/2], alpha = 0.6,
                        colors=smart_colours[sname[0]]['dark'],
                        linewidths=linewidths, zorder=0.25)
        
        #Shade selected obs
        if args.shade_light:
            if sobsid in args.shade_light:
                cs = plt.tricontour(nx, ny, nz, levels=max(nz)/2, alpha=0.0)
                cs0 = cs.collections[0]
                cspaths = cs0.get_paths()
                for cspath in cspaths:
                    spch_0 = patches.PathPatch(cspath, facecolor=smart_colours[sname[0]]['dark'],
                                                edgecolor='gray', lw=0.5, alpha=0.3, zorder=0.5)
                    ax.add_patch(spch_0)

        if args.shade:
            if sobsid in args.shade:
                cs = plt.tricontour(nx, ny, nz, levels=max(nz)/2, alpha=0.0)
                cs0 = cs.collections[0]
                cspaths = cs0.get_paths()
                for cspath in cspaths:
                    spch_0 = patches.PathPatch(cspath, facecolor=smart_colours[sname[0]]['dark'],
                                                edgecolor=smart_colours[sname[0]]['dark'], lw=0.5, alpha=0.85, zorder=0.5)
                    ax.add_patch(spch_0)

    # Add pulsars
    if args.pulsar:
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
        ax.scatter(ra_PCAT, dec_PCAT, s=5, color ='r', zorder=2)

    if args.pulsar_cand:
        ra_PCAT = []
        dec_PCAT = []
        for pulsar in args.pulsar_cand:
            pulsar = pulsar.split('_')
            ra_temp, dec_temp = sex2deg(pulsar[0], pulsar[1])
            if args.ra_offset:
                if ra_temp > 180:
                    ra_PCAT.append(-ra_temp/180.*np.pi+2*np.pi)
                else:
                    ra_PCAT.append(-ra_temp/180.*np.pi)
            else:
                ra_PCAT.append(-ra_temp/180.*np.pi+np.pi)
            dec_PCAT.append(dec_temp/180.*np.pi)
        ax.scatter(ra_PCAT, dec_PCAT, s=5, color ='g', zorder=2)

    plt.xlabel("Right Ascension")
    plt.ylabel("Declination")

    #xtick_labels = ['0h','2h','4h','6h','8h','10h','12h','14h','16h','18h','20h','22h']
    if args.ra_offset:
        xtick_labels = ['10h', '8h', '6h', '4h', '2h', '0h', '22h', '20h', '18h', '16h', '14h']
        xticks = [150., 120., 90., 60., 30., 0., 330., 300., 270., 240., 210. ]
    else:
        xtick_labels = [ '22h', '20h', '18h', '16h', '14h','12h','10h', '8h', '6h', '4h', '2h']
        xticks = [330., 300., 270., 240., 210., 180., 150., 120., 90., 60., 30.]

    ax.set_xticklabels(xtick_labels, zorder=150)
    plt.grid(True, color='gray', lw=0.5, linestyle='dotted')

    plot_type = args.plot_type
    plot_name = "SMART_progress"
    print("saving {}.{}".format(plot_name, plot_type))
    fig.savefig(plot_name + '.' + plot_type, format=plot_type, dpi=1000, bbox_inches='tight')