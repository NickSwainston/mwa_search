#!/usr/bin/env python3

import logging
import argparse
from job_submit import submit_slurm

logger = logging.getLogger(__name__)


#----------------------------------------------------------------------
class run_params_class:

    def __init__(self, pointing_dir=None, cal_id=None,\
                obsid=None, pulsar=None, threshold=10.0,\
                stop=False, loglvl="INFO", mode=None,\
                mwa_search="master", vcs_tools="multi-pixel_beamform",\
                nbins=128, subint=10.0):

        self.pointing_dir   = pointing_dir
        self.cal_id         = cal_id
        self.obsid          = obsid
        self.pulsar         = pulsar
        self.stop           = stop
        self.loglvl         = loglvl
        self.threshold      = threshold
        self.mode           = mode
        self.mwa_search     = mwa_search
        self.vcs_tools      = vcs_tools
        self.nbins          = nbins
        self.subint         = subint

        if self.obsid==None:
            mydict=info_from_dir(self.pointing_dir)
            if self.obsid==None:
                self.obsid=mydict["obsid"]

    def single_pointing(self):
        self.pointing_dir=self.pointing_dir[0]
#----------------------------------------------------------------------
def info_from_dir(pointing_dir):

    #given a pointing directory, returns a dictionary containing the obsid and pulsar name
    mydict = {}
    mydict["obsid"] = pointing_dir.split("/")[4]
    pulsar_coords = pointing_dir.split("/")[6].split("_")
    pulsar_1 = pulsar_coords[0].split(":")[0] + pulsar_coords[0].split(":")[1]
    pulsar_2 = pulsar_coords[1].split(":")[0] + pulsar_coords[1].split(":")[1]
    mydict["pulsar"] = "J" + pulsar_1 + pulsar_2

    return mydict

#----------------------------------------------------------------------
def stokes_fold(run_params):

    logger.info("Initilizing stokes fold")
    commands=[]
    commands.append("stokes_fold.py -d {0} -p {1} -b {2} -s {3} -L {4} --vcs_tools {5} --mwa_search {6}"\
            .format(run_params.pointing_dir, run_params.pulsar, run_params.nbins, run_params.subint,\
            run_params.loglvl, run_params.vcs_tools, run_params.mwa_search))
    
    name="Stokes_Fold_{0}_{1}".format(run_params.pulsar, run_params.obsid)
    batch_dir = "/group/mwaops/vcs/{0}/batch/".format(run_params.obsid)

    submit_slurm(name, commands,\
                batch_dir=batch_dir,\
                slurm_kwargs={"time": "2:00:00"},\
                module_list=["mwa_search/{0}".format(run_params.mwa_search)],\
                submit=True, vcstools_version="{0}".format(run_params.vcs_tools))


#----------------------------------------------------------------------
def binfind(run_params):

    if run_params.stop==False:
        cont = "--launch_next"
    else:
        cont = ""

    logger.info("pulsar and obsid: {0}, {1}".format(run_params.pulsar, run_params.obsid))

    #Run binfinder.py
    commands = []
    if run_params.mode=='f':
        commands.append("binfinder.py -m f -d {0} -O {1} -t {2} -p {3} -o {4} -L {5} --mwa_search {6}\
        --vcs_tools {7} {8}"\
        .format(run_params.pointing_dir, run_params.cal_id, run_params.threshold, run_params.pulsar,\
        run_params.obsid, run_params.loglvl, run_params.mwa_search, run_params.vcs_tools, cont))

    elif run_params.mode=='p':
        commands.append("binfinder.py -m c -d {0} -O {1} -t {2} -p {3} -o {4} -L {5} --mwa_search {6}\
        --vcs_tools {7} {8}"\
        .format(run_params.pointing_dir, run_params.cal_id, run_params.threshold, run_params.pulsar,\
        run_params.obsid, run_params.loglvl, run_params.mwa_search, run_params.vcs_tools, cont))

    elif run_params.mode=="m":
        pointing_string=""
        for p in run_params.pointing_dir:
            logger.info("folding on: {0}".format(p))
            pointing_string = pointing_string + p + " "
        commands.append("binfinder.py -m m -d {0} -O {1} -t {2} -p {3} -o {4} -L {5} --mwa_search {6}\
        --vcs_tools {7} {8}"\
        .format(pointing_string, run_params.cal_id, run_params.threshold, run_params.pulsar,\
        run_params.obsid, run_params.loglvl, run_params.mwa_search, run_params.vcs_tools, cont))


    name = "binfinder_{0}_startup".format(run_params.pulsar)
    batch_dir = batch_dir = "/group/mwaops/vcs/{0}/batch/".format(run_params.obsid)
    submit_slurm(name, commands,\
                batch_dir=batch_dir,\
                slurm_kwargs={"time": "00:05:00"},\
                module_list=['mwa_search/{0}'.format(run_params.mwa_search),\
                            'presto/no-python'],\
                submit=True, vcstools_version="{0}".format(run_params.vcs_tools))

    logger.info("Job successfully submitted")

#----------------------------------------------------------------------
if __name__ == '__main__':
    #dictionary for choosing log-levels
    loglevels = dict(DEBUG=logging.DEBUG,
                    INFO=logging.INFO,
                    WARNING=logging.WARNING,
                    ERROR = logging.ERROR)


    #Arguments

    parser = argparse.ArgumentParser(description="""A pipeline for processing calibrated VCS data""")


    obsop = parser.add_argument_group("Observation Options")
    obsop.add_argument("-d", "--pointing_dir", nargs='+', help="The location of the pointing directory/s")
    obsop.add_argument("-o", "--obsid", type=str, help="The obs ID of the data")
    obsop.add_argument("-O", "--cal_id", type=str, help="The ID of the calibrator used to calibrate the data")
    obsop.add_argument("-p", "--pulsar", type=str, help="The J name of the pulsar. e.g. J2241-5236")


    binfindop = parser.add_argument_group("Binfinder Options")
    binfindop.add_argument("-t", "--threshold", type=float, default=10.0, help="The presto sigma value\
                             above which is deemed a detection. If this value is not exceeded in any\
                             of the folds, the pipeline will terminate. However, a value of 5.0 will\
                             always be checked and uploaded to the pulsar database if it exists.\
                             Default: 10.0")

    stokesop = parser.add_argument_group("Stokes Fold Options")
    stokesop.add_argument("-b", "--nbins", type=int, default=128, help="The number of bins for to fold over for the stokes folding script. Default: 128")
    stokesop.add_argument("-s", "--subint", type=float, default=10.0, help="The length of the integrations (in seconds) used for dspsr. Default: 10.0")

    otherop = parser.add_argument_group("Other Options")
    otherop.add_argument("-L", "--loglvl", type=str, default="INFO", help="Logger verbosity level. Default: INFO", choices=loglevels.keys())
    otherop.add_argument("-S", "--stop", action="store_true", help="Use this mode to tell the pipeline not to continue processing data after finishing the desired task")
    otherop.add_argument("--mwa_search", type=str, default="master",  help="The version of mwa_search to use. Default: master")
    otherop.add_argument("--vcs_tools", type=str, default="multi-pixel_beamform", help="The version of vcs_tools to use. Default: master")
    otherop.add_argument("-m", "--mode", type=str, help="The mode in which to run this script:\n\
                        'p' - Folds on a small number of bins in order to check if a pulsar is\n\
                         detected in the given pointing directory (runs the b mode after by default)\n\
                        'm' - Folds on  any numer of pointings and finds the one with the best detection\
                        (runs the b mode after by defualt)\n\
                        'b' - attempts to find an adequate number of bins to fold the pulsar with\n\
                         and outputs a bestprof file (runs the s mode after by default)\n\
                        's' - will fold across stokes IQUV and attempt to find the rotation measure")


    args = parser.parse_args()

    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    ch.setLevel(loglevels[args.loglvl])
    formatter = logging.Formatter('%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)


    run_params = run_params_class(pointing_dir=args.pointing_dir, cal_id=args.cal_id,\
                                pulsar=args.pulsar, obsid=args.obsid, stop=args.stop,\
                                mode=args.mode, mwa_search=args.mwa_search,\
                                vcs_tools=args.vcs_tools, loglvl=args.loglvl,\
                                threshold=args.threshold, nbins=args.nbins,\
                                subint=args.subint)

    if run_params.mode is not 'm':
        run_params.single_pointing()

    if run_params.mode=="b" or run_params.mode=="p" or run_params.mode=="m":
        binfind(run_params)
    elif run_params.mode=="s":
        stokes_fold(run_params)
    else:
        logger.error("Mode not recognized. Please rerun with a valid mode identifer")
