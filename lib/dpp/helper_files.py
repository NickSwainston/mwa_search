import logging
from os import symlink, rmdir, unlink, getcwd, chdir, remove
from os.path import exists, join, basename
from shutil import copyfile
from glob import glob

from vcstools.config import load_config_file
from vcstools.general_utils import mdir

comp_config = load_config_file()
logger = logging.getLogger(__name__)


def create_dpp_dir(kwargs):
    dpp_dir = join(comp_config["base_data_dir"], str(kwargs["obsid"]), "dpp")
    mdir(dpp_dir, dpp_dir)


def setup_cfg_dirs(cfg):
    """Creates the necessary folders and symlinks for dpp"""
    # Create pulsar directory
    mdir(cfg["files"]["psr_dir"], cfg["files"]["psr_dir"])
    # Create classify dir
    mdir(cfg["files"]["classify_dir"], cfg["files"]["classify_dir"])
    # Create edited .eph if necessary
    if cfg["source"]["binary"]:
        with open(cfg["source"]["edited_eph_name"], "w") as f:
            f.write(cfg["source"]["edited_eph"])
    # Create symlinks to pointing dirs
    remove_pointings = []
    for pointing in cfg['folds'].keys():
        real = join(comp_config["base_data_dir"], str(cfg["obs"]["id"]), "pointings", pointing)
        sym = join(cfg["files"]["psr_dir"], pointing)
        if exists(real):
            if not exists(sym):
                symlink(real, sym)
        else: # Remove the pointing from the dictionary if the real pointing directory isn't found
            logger.warn(f"Expected pointing directory not found: {pointing}. Skipping")
    for pointing in remove_pointings:
        del cfg["folds"][pointing]


def clean_cfg(cfg):
    """Remove any lists in the cfg that are empty and deletes the directories"""
    if not cfg["folds"]:
        logger.warn(f"No pointings available for {cfg['source']['name']}. Removing")
        cfg = None # Delete directory


def remove_old_results(cfg):
    """Removes old results from previous ppp runs"""
    # Remove everything in classify dir
    for f in glob(join(cfg["files"]["classify_dir"], "*")):
        remove(f)
    # Remove pfds
    for pointing in cfg["folds"].keys():
        pfds = glob(join(cfg["files"]["psr_dir"], f"*{cfg['files']['file_precursor']}*.pfd*"))
        for pfd in pfds:
            remove(pfd)
    # Remove pngs
    pngs = glob(join(cfg["files"]["psr_dir"], f"*{cfg['files']['file_precursor']}*.png"))
    for png in pngs:
        remove(png)
    # Remove Postscripts
    pscripts = glob(join(cfg["files"]["psr_dir"], f"*{cfg['files']['file_precursor']}*.ps"))
    for pscript in pscripts:
        remove(pscript)
    # Remove other various files
    files_to_remove = [
        cfg["files"]["archive"],
        cfg["files"]["archive_ascii"],
        cfg["files"]["converted_fits"],
        cfg["files"]["debased_fits"],
        cfg["files"]["paswing"],
        cfg["files"]["RVM_fit_initial"],
        cfg["files"]["RVM_fit_final"]
    ]
    for f in files_to_remove:
        try:
            remove(f)
        except FileNotFoundError as e:
            pass


def file_precursor(kwargs, psr):
    """
    Creates a common precursor string for files and directories
    using kwargs from observation_processing_pipeline.py and a pulsar name
    """
    label = kwargs["label"]
    if label:
        label = f"_{kwargs['label']}"
    return f"{kwargs['obsid']}{label}_{psr}"


def setup_classify(cfg):
    """Creates the required directories and copies files for the lotaas classifier"""
    owd = getcwd()
    chdir(cfg["files"]["psr_dir"])
    mdir(cfg["files"]["classify_dir"], cfg["files"]["classify_dir"]) # This should already exist but keep it anyway
    for pointing in cfg["folds"].keys():
        init_bins = list(cfg["folds"][pointing]["init"].keys())[0]
        if int(init_bins) not in (50, 100):
            raise ValueError(f"Initial bins for {cfg['source']['name']} is invalid: {init_bins}")
        pfd_name = glob_pfds(cfg, pointing, init_bins, pfd_type=".pfd")[0]
        # Copy pdf file to classify directory
        newfilename=join(cfg["files"]["classify_dir"], basename(pfd_name))
        copyfile(pfd_name, newfilename)
    chdir(owd)


def find_config_files(obsid, label=""):
    """Searches the obsid/dpp directories to find any config (.yaml) files"""
    dpp_dir = join(comp_config["base_data_dir"], str(obsid), "dpp")
    yaml_files = join(dpp_dir, "*", f"{obsid}*{label}.yaml")
    config_pathnames = glob(yaml_files)
    if not config_pathnames:
        raise ValueError(f"No config files found: {yaml_files}")
    return config_pathnames


def glob_pfds(cfg, pointing, bins, pfd_type=".pfd"):
    """Globs the appropriate directory for the given pointing and bins for .pfds and returns the list"""
    # See helper_prepfold.generate_prep_name() for glob dir reference
    glob_dir = join(f"*{cfg['files']['file_precursor']}*{pointing}*b{bins}*{pfd_type}")
    return glob(glob_dir)