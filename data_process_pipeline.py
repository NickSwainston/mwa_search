#!/usr/bin/env python3

import logging
import argparse
from job_submit import submit_slurm

logger = logging.getLogger(__name__)


#----------------------------------------------------------------------
class run_params_class:

    def __init__(self, pointing_dir, cal_id,
                obsid=None, pulsar=None, threshold=10.0,
                stop=False, loglvl="INFO",
                mode=None):

        self.pointing_dir   = pointing_dir
        self.cal_id         = cal_id
        self.obsid          = obsid
        self.pulsar         = pulsar
        self.stop           = stop
        self.stop_stokes    = stop_stokes
        self.threshold      = threshold
        self.mode           = mode

        if self.obsid==None:
            mydict=info_from_dir(self.pointing_dir)
            if self.obsid==None:
                self.obsid=mydict["obsid"]

    def single_pointing(self):
        self.pointing_dir=pointing_dir[0]
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
#def stokes_fold(pointing_dir, cal_id)



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
        commands.append("binfinder.py -m f -d {0} -O {1} -t {2} -p {3} -o {4} {5}".format(run_params.pointing_dir, run_params.cal_id, run_params.threshold, run_params.pulsar, run_params.obsid, cont))
    elif run_params.mode=='p':
        commands.append("binfinder.py -m c -d {0} -O {1} -t {2} -p {3} -o {4} {5}".format(run_params.pointing_dir, run_params.cal_id, run_params.threshold, run_params.pulsar, run_params.obsid, cont))
    elif run_params.mode=="m":
        commands.append("binfinder.py -m m -d {0} -O {1} -t {2} -p {3} -o {4} {5}".format(run_params.pointing_dir, run_params.cal_id, run_params.threshold, run_params.pulsar, run_params.obsid, cont))


    name = "binfinder_{0}_startup".format(run_params.pulsar)
    batch_dir = batch_dir = "/group/mwaops/vcs/{0}/batch/".format(run_params.obsid)
    submit_slurm(name, commands,
                batch_dir=batch_dir,
                slurm_kwargs={"time": "00:05:00"},
                module_list=['mwa_search/k_smith',
                            'presto/no-python'],
                submit=True, vcstools_version="multi-pixel_beamform")

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
    binfindop.add_argument("-t", "--threshold", type=float, default=10, help="The presto sigma value\
                             above which is deemed a detection. If this value is not exceeded in any\
                             of the folds, the pipeline will terminate. However, a value of 5.0 will\
                             always be checked and uploaded to the pulsar database if it exists.\
                             Default=10.0")


    otherop = parser.add_argument_group("Other Options")
    otherop.add_argument("-L", "--loglvl", type=str, default="INFO", help="Logger verbosity level. Default: INFO", choices=loglevels.keys())
    otherop.add_argument("-S", "--stop", action="store_true", help="Use this mode to tell the pipeline not to continue processing data after finishing the desired task")
    otherop.add_argument("-m", "--mode", type=str, help="The mode in which to run this script:\n\
                        'p' - Folds on a small number of bins in order to check if a pulsar is\n\
                         detected in the given pointing directory (runs the b mode after by default)\n\
                        'm' - Folds on many pointings and finds the one with the best detection\
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


    run_params = run_params_class(args.pointing_dir, args.cal_id, pulsar=args.pulsar,obsid=args.obsid,\
                                 stop=args.stop, mode=args.mode)

    if run_params.mode is not 'm':
        run_params.single_pointing()

    if run_params.mode=="b" or run_params.mode=="p" or run_params.mode=="m":
        binfind(run_params)
    elif run_params.mode=="s":
        stokes_fold(args.pointing_dir, args.cal_id)
    else:
        logger.error("Mode not recognized. Please rerun with a valid mode identifer")
