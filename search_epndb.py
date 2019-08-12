import os
import glob
import logging
import argparse
import sys
from job_submit import submit_slurm
import ast
logger = logging.getLogger(__name__)

def dload_dbase(cores = 8, path = "/group/mwaops/k_smith/"):

    epn_url = "http://www.epta.eu.org/epndb/json/"
    wget_cmd = ""
    for _ in range(cores-1):
        wget_cmd = wget_cmd + "wget -r -np -N {0} & ".format(epn_url)
    wget_cmd = wget_cmd + "wget -r -np -N {0}".format(epn_url)

    #create slurm job:
    commands = []
    commands.append("cd {0}".format(path))
    commands.append('echo "Downloading epn database to {0}"'.format(path))
    commands.append("srun --export=all -N 1 -n 1 {0}".format(wget_cmd))


    slurm_kwargs = {"partition": "gpuq",
                    "time": "12:00:00",
                    "nodes": "1",
                    "ntasks-per-node": "{0}".format(cores),
                    "gres": "gpu:1"}

    name = "epn_dload"
    batch_dir = path + "/wget_batch"
    submit_slurm(name, commands,
                 batch_dir=batch_dir,
                 slurm_kwargs=slurm_kwargs,
                 submit=True, vcstools_version="master")



def get_epn_paths(pulsar_name, path_to_database):
    #first check if the epn dictionary has been made. If not, make it
    if "/epn_database_dictionary.dat" not in glob.glob(path_to_database):
        update_epn_dict(path_to_database)
    #open and read the dictionary into a dictionary variable
    f = open(path_to_database + "/epn_database_dictionary.dat", "r")
    line = f.read()
    #epn_dict = {}
    #epn_dict = eval(line)
    pn_dict = ast.literal_eval(line)
    f.close()
    return epn_dict


def update_epn_dict(path_to_database):

    logger.info("Creating a dictionary of pulsars from the epn database...")

    #begin looking through the databse:
    author_dirs = []
    for adir in glob.glob(path_to_database + "/*/"):
        author_dirs.append(adir)

    #Now check all of the directories we just found for pulsar directories
    pulsar_dirs = []
    for auth_dir in author_dirs:
        for puls_dir in glob.glob(auth_dir + "*/"):
            pulsar_dirs.append(puls_dir)

    #Now find all .json files in every Directory
    json_paths = []
    for puls_dir in pulsar_dirs:
        for filepath in glob.glob(puls_dir + "*"):
            if ".json" in filepath:
                json_paths.append(filepath)

    #For each filepath we also want the associated pulsar's name:
    all_pulsars = []
    for filepath in json_paths:
        pname = os.path.dirname(filepath)
        pname = os.path.basename(pname)
        all_pulsars.append(pname)

    #Now we can create a dictionary of all of the pulsar_dirs
    epn_dict = {}
    for i, pulsar in enumerate(all_pulsars):
        if pulsar in epn_dict:
            paths = epn_dict[pulsar]
            paths.append(json_paths[i])
        else:
            epn_dict[pulsar] = []
            epn_dict[pulsar].append(json_paths[i])

    logger.info("Number of unique pulsars found: {0}".format(len(epn_dict.keys())))
    logger.info("Writing the dictionary to: {0}".format(path_to_database))

    f = open(path_to_database + "/epn_database_dictionary.dat", "w+")
    f.write(str(epn_dict))
    f.close()


if __name__ == '__main__':

    #dictionary for choosing log-levels
    loglevels = dict(DEBUG=logging.DEBUG,
                    INFO=logging.INFO,
                    WARNING=logging.WARNING,
                    ERROR = logging.ERROR)

    #Arguments
    parser = argparse.ArgumentParser(description="""Contains functions to use the epn database""")
    parser.add_argument("-L", "--loglvl", type=str, help="Logger verbosity level. Default: INFO", choices=loglevels.keys(), default="INFO")
    parser.add_argument("-p", "--path_to_database", type=str, help="Where the json directory is located")
    parser.add_argument("-m", "--mode", type=str, help="Mode selsction. 'epn' - creates or updates epn dictionary from the databse at the given output directory.\n'update' - forces an update of the epn database")

    args = parser.parse_args()

    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    ch.setLevel(loglevels[args.loglvl])
    formatter = logging.Formatter('%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    if args.mode=="epn":
        update_epn_dict(args.path_to_database)
    elif args.mode=="update":
        dload_dbase()
    else:
        logger.error("No valid mode chosen. Please choose a mode and try again. For more info use the '-h' command.")
        sys.exit(0)
