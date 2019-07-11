#!/usr/bin/env python3

import os
import glob
import logging
import argparse
import sys
import logging
from job_submit import submit_slurm
import binfinder
import subprocess
logger = logging.getLogger(__name__)

def ephem_to_file(pulsar, fname=None)
    if fname==None:
        fname="{0}.eph".format(pulsar)
    
    f = open(fname, "w+")
    subprocess.call(["psrcat", "-e", "{0}".format(pulsar)], stdout=f)
    f.close()


def psrcat_RM(pulsar)

    sp = subprocess.Popen(["psrcat", "-e"," {0}".format(pulsar)], stdout=subprocess.PIPE)
    out = sp.stdout.read()
    out = out.decode("utf-8")
    for i in out:
        measurement = i.split()
        if measurement[0]=="RM"
            RM_found=True
            break
    
    if RM_found==True:
        RM = float(measurement[1])
        RM_err = float(measurement[2])    
        logger.info("Rotation measure found from psrcat ephemeris") 
        return RM, RM_err
    else:
        logger.warn("Rotation measure not stored in psrcat ephemeris")
        RM = None
        RM_err = None
    
    return RM, RM_err


def find_RM_from_archive(archive_file)
    
    sp = subprocess.Popen(["rmfit", "{0}","-t".format(archive_file)], stdout=subprocess.PIPE)
    out = sp.stdout.read()
    out = out.decode("utf-8").split()
  
    if out:
        #An RM was found 
        RM = out[3]
        RM_err = out[5]
    else:
        #RM not found 
        RM = None
        RM_err = None    
        logger.warn("RM could not be generated from archive")

    return RM, RM_err

#def submit_RM()
    #This is not implemented in the pulsar databse yet

def make_archive(pointing_dir, pulsar, nbins=None, subint=10.0)
    
    if nbins==None:
        logger.warn("Number of bins not provided. Attempting to find the best profile")
        try:
            prof_name = binfinder.get_best_profile(pointing_dir, 10.0) 
            nbins = bindinder.bestprof_info(filename=prof_name)["nbins"]
        except RuntimeError as err:
            logger.error("Number of bins could not be found. Full error:\n{0}".format(err))
            logger.error("Exiting...")
            sys.exit(1)
            
    
    #print ephemeris to file
    ephem_to_file(pulsar)
    #make dspsr archive
    subprocess.call(["dspsr", "-cont", "-U", "4000", "-A", "-L", "{0}".format(subint), "-E", "{1}.eph".format(pulsar), "-K", "-b", "{0}".format(nbins), "-O", "{0}_subint_{1}".format(pulsar, subint), "*.fits"
    #Attempt to find rotation measure from archive file
    RM, RM_err = find_RM_from_archive("{0}_subint_{1}.ar".format(pulsar, subint)) 
    if RM==None:
        RM, RM_err = psrcat_RM(pulsar)
        if RM==None:
            logger.info("RM for pulsar {0} not found on record. Cannot continue with polarimetry".format(pulsar))
            sys.exit(0)
    #else: 
       # submit_RM
    
    #Correct for the RM:
    subprocess.call(["pam", "-e", "ar2", "-R", "{0}".format(RM), "{0}_subint_{1}".format(pulsar, subint)])
    #Turn the archive into a readable ascii file
    fname="{0}_pol_prof.txt".format(pulsar)
    f = open(fname, "w+")
    subprocess.call(["pdv", "-FTt", "{0}_subint_{1}.ar2".format(pulsar, subint)   
    f.close()

    return fname

if __name__ == '__main__':
    #dictionary for choosing log-levels
    loglevels = dict(DEBUG=logging.DEBUG,
                    INFO=logging.INFO,
                    WARNING=logging.WARNING,
                    ERROR = logging.ERROR)

    #Arguments
    parser = argparse.ArgumentParser(description="""Folds across stokes IQUV and attempts to find the RM""")

    parser.add_argument("-d", "--pointing_dir", type=str, help="Pointing directory that contains the fits files")
    parser.add_argument("-L", "--loglvl", type=str, default="INFO", help="Logger verbosity level. Default: INFO", choices=loglevels.keys())
    parser.add_argument("-p", "--pulsar", type=str, default=None, help="The name of the pulsar. If not provided, the script will try to get the pulsar from the pointing directory"
    parser.add_argment("-b", "--nbins", typ=int, default=None, help="The number of bins to fold the profile with" )
    args = parser.parse_args()

    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    ch.setLevel(loglevels[args.loglvl])
    formatter = logging.Formatter('%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    os.chdir(args.pointing_dir)

    fname = make_archive(args.pointing_dir, args.pulsar, nbins=args.nbins)
    if fname is not None:
        plot_pol(fname) #TODO: implement this
