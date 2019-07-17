#!/usr/bin/env python3

import os
import glob
import logging
import argparse
import sys
import logging
import subprocess

from job_submit import submit_slurm
import binfinder
import data_process_pipeline
import plotting_toolkit
logger = logging.getLogger(__name__)


class run_params_class:
    
    def __init__(self, pointing_dir=None, pulsar=None,\
                nbins=None, obsid=None, cal_id=None,\
                subint=10.0, loglvl="INFO", noplot=False,\
                vcs_tools="multi-pixel_beamform",\
                mwa_search="master"):

        self.pointing_dir       = pointing_dir
        self.pulsar             = pulsar
        self.nbins              = nbins
        self.obsid              = obsid
        self.cal_id             = cal_id
        self.subint             = subint
        self.loglvl             = loglvl
        self.vcs_tools          = vcs_tools
        self.mwa_search         = mwa_search
        self.noplot             = noplot

        if self.obsid is None:
            self.obsid = data_process_pipeline.info_from_dir(self.pointing_dir)["obsid"]
    

def ephem_to_file(pulsar, fname=None):
    #prints ephemeris to a file
    if fname is None:
        fname="{0}.eph".format(pulsar)
    
    f = open(fname, "w+")
    subprocess.call(["psrcat", "-e", "{0}".format(pulsar)], stdout=f)
    f.close()


def psrcat_RM(pulsar):
    #Tries to get RM from psrcat
    sp = subprocess.Popen(["psrcat", "-e"," {0}".format(pulsar)], stdout=subprocess.PIPE)
    out = sp.stdout.read()
    out = out.decode("utf-8")
    for i in out:
        measurement = i.split()
        if measurement[0]=="RM":
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


def find_RM_from_archive(archive_file):
   #Tries to get RM from a dspsr archive file 
    sp = subprocess.Popen(["rmfit", "{0}","-t".format(archive_file)], stdout=subprocess.PIPE)
    out = sp.stdout.read()
    out = out.decode("utf-8").split()
  
    if out:
        #An RM was found
        #TODO: submit_RM(RM, RM_err) 
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


def make_archive(run_params):
   #Tries to make an RM fixed archive file using dspsr 
    #failsafe:
    if run_params.nbins is None:
        logger.warn("Number of bins not provided. Attempting to find the best profile")
        try:
            prof_name = binfinder.get_best_profile(pointing_dir, 10.0) 
            nbins = bindinder.bestprof_info(filename=prof_name)["nbins"]
        except RuntimeError as err:
            logger.error("Number of bins could not be found. Full error:\n{0}".format(err))
            logger.error("Exiting...")
            sys.exit(1)
            
    
    #print ephemeris to file
    ephem_to_file(run_params.pulsar)
    #make dspsr archive
    subprocess.run(["dspsr", "-h"])
    subprocess.run(["dspsr", "-cont", "-U", "4000", "-A", "-L", "{0}".format(run_params.subint),\
    "-E", "{0}.eph".format(run_params.pulsar), "-K", "-b", "{0}".format(run_params.nbins), "-O",\
    "{0}_subint_{1}".format(run_params.pulsar, run_params.subint), "*.fits"])

    #Attempt to find rotation measure from archive file
    RM, RM_err = find_RM_from_archive("{0}_subint_{1}.ar".format(run_params.pulsar, run_params.subint)) 
    if RM is None:
        RM, RM_err = psrcat_RM(run_params.pulsar)
        if RM is None:
            logger.info("RM for pulsar {0} not found on record. Cannot continue with polarimetry".format(run_params.pulsar))
            sys.exit(0)
    #else: 
       # submit_RM
    
    #Correct for the RM:
    subprocess.run(["pam", "-e", "ar2", "-R", "{0}".format(RM), "{0}_subint_{1}".format(run_params.pulsar, run_params.subint)])
    #Turn the archive into a readable ascii file
    fname="{0}_pol_prof.txt".format(run_params.pulsar)
    f = open(fname, "w+")
    subprocess.run(["pdv", "-FTt", "{0}_subint_{1}.ar2".format(run_params.pulsar, run_params.subint)], stdout=f)   
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

    foldop = parser.add_argument_group("Folding Options:")
    foldop.add_argument("-d", "--pointing_dir", type=str, help="Pointing directory that contains the fits files")
    foldop.add_argument("-p", "--pulsar", type=str, default=None, help="The name of the pulsar. If not provided, the script will try to get the pulsar from the pointing directory")
    foldop.add_argument("-b", "--nbins", type=int, default=128, help="The number of bins to fold the profile with")
    foldop.add_argument("-s", "--subint", type=float, default=10.0, help="The length of the integrations in seconds. Default: 10.0")   
 
    plotop = parser.add_argument_group("Plotting Options:")
    plotop.add_argument("--noplot", action="store_true", help="Use this tag to not output a plot")

    otherop = parser.add_argument_group("Other Options:")
    otherop.add_argument("-L", "--loglvl", type=str, default="INFO", help="Logger verbosity level. Default: INFO", choices=loglevels.keys())
    otherop.add_argument("--vcs_tools", type=str, default="multi-pixel_beamform", help="The version of vcstools to use. Default: multi-pixel_beamform")
    otherop.add_argument("--mwa_search", type=str, default="master", help="The version of mwa_search to use. Default: master")

    args = parser.parse_args()

    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    ch.setLevel(loglevels[args.loglvl])
    formatter = logging.Formatter('%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)


    run_params = run_params_class(pointing_dir=args.pointing_dir, pulsar=args.pulsar,\
                nbins=args.nbins, noplot=args.noplot, loglvl=args.loglvl, subint=args.subint)
    

    if run_params.pulsar is None:
        logger.error("Pulsar name not supplied. Please run again and specify puslar name")
        sys.exit(1)

    fname = make_archive(run_params)
    if fname is not None:
        plotting_toolkit.plot_archive(obsid=run_params.obsid, archive=fname, pulsar=run_params.pulsar, out_dir=run_params.pointing_dir)
