from os.path import join as ospj
import logging

from binfinder import NoSuitableProfileError, bestprof_info
from dpp import plotting_toolkit
import data_processing_pipeline as dpp
from config_vcs import load_config_file


comp_config = load_config_file()
logger = logging.getLogger(__name__)


def submit_to_db(pipe, dep_id=None, dep_type="afterok"):
    """Submits the best fold profile to the pulsar database. Will also submit .ppps, .png and .pfd"""
    bin_list = pipe["folds"]["init"].keys().append(
        pipe["folds"]["post"].keys())
    commands = []
    commands.append(f"cd {pipe['run_ops']['my_dir']}")
    for bin_count in bin_list:
        ppps = f"*_b{bin_count}*_{pipe['source']['name'][1:]}*.pfd.ps"
        ppps = glob.glob(ospj(pipe["run_ops"]["my_dir"], ppps))[0]
        bestprof = f"*_b{bin_count}*_{pipe['source']['name'][1:]}*".format(
            best_bins)
        bestprof = glob.glob(ospj(
            pipe["run_ops"]["my_dir"], bestprof))[0]
        png = f"*_b{bin_count}*_{pipe['source']['name'][1:]}*.png".format(
            best_bins)
        png = glob.glob(ospj(pipe["run_ops"]["my_dir"], png))[0]
        pfd = f"*_b{bin_count}*_{pipe['source']['name'][1:]}*.pfd".format(
            best_bins)
        pfd = glob.glob(ospj(pipe["run_ops"]["my_dir"], pfd))[0]
        commands.append(
            f"echo 'Submitting profile to database with {bin_count} bins'")
        commands.append(
            f"submit_to_database.py -o {pipe['obs']['id']} --cal_id {pipe['obs']['cal_id']} -p {pipe['source']['name']} --bestprof {bestprof} --ppps {ppps}")
        # Also Make a nice plot
        plotting_toolkit.plot_bestprof(ospj(
            pipe["run_ops"]["my_dir"], bestprof, out_dir=pipe["run_ops"]["my_dir"]))

    name = f"Submit_db_{pipe['source']['name']}_{pipe['obs']['id']}"
    batch_dir = ospj(
        comp_config['base_data_dir'], pipe['obs']['id'], "batch")
    this_id = submit_slurm(name, commands,
                           batch_dir=batch_dir,
                           slurm_kwargs={"time": "02:00:00"},
                           depend=dep_id,
                           module_list=[
                               f"mwa_search/{pipe['run_ops']['mwa_search']}", f"mwa_search/{pipe['run_ops']['vcstools']}"],
                           submit=True, depend_type=dep_type)

    logger.info(f"Submission script on queue for profile: {bestprof}")
    logger.info(f"Job name: {name}")
    logger.info(f"Dependenices: {dep_id}")
    logger.info(f"Depend type: {dep_type}")
    logger.info(f"Job ID: {this_id}")

    return this_id


def move_product_dir(pipe, dep_id=None, dep_type="afterok"):
    """Creates a new folder for the data products"""
    current_dir = pipe["run_ops"]["dir"]
    move_loc = ospj(load_from_config()[
                    "base_data_dir"], pipe["obs"]["id"], "data_products")
    pointing = [i for i in current_dir.split("/") if i != ""]
    new_dir = os.path.join(move_loc, pointing[-1])

    commands = []
    commands.append(f"echo 'Moving {current_dir} to new location: {new_dir}'")
    if pipe["source"]["name"][-1].isalpha():
        commands.append(f"cp -ru {current_dir} {new_dir}")
    else:
        commands.append(f"mv {current_dir} {new_dir}")

    name = f"Move_{pipe['source']['name']}_{pipe['obs']['id']}"
    batch_dir = os.path.join(
        comp_config['base_data_dir'], pipe['obs']['id'], "batch")
    this_id = submit_slurm(name, commands,
                           batch_dir=batch_dir,
                           slurm_kwargs={"time": "02:00:00"},
                           depend=dep_id,
                           submit=True, depend_type=dep_type)

    logger.info(f"Move script submitted for pointing directory: {current_dir}")
    logger.info(f"Job name: {name}")
    logger.info(f"Dependenices: {dep_id}")
    logger.info(f"Depend type: {dep_type}")
    logger.info(f"Job ID: {this_id}")


def fill_pipe_folds(pipe):
    """Fills the pipe with the bestprof information from all available folds"""
    for bins in pipe["folds"]["init"]:
        bestprof = glob.glob(ospj(pipe["run_ops"]["my_dir"],
            f"*b{bins}**{pipe['source']['name'][1:]}*.bestprof"))[0]
        info = bestprof_info(bestprof)
        pipe["folds"]["init"][bins] = info
    for bins in pipe["folds"]["post"]:
        bestprof = glob.glob(ospj(pipe["run_ops"]["my_dir"],
            f"*b{bins}**{pipe['source']['name'][1:]}*.bestprof"))[0]
        info = bestprof_info(bestprof)
        pipe["folds"]["post"][bins] = info


def get_best_fold(pipe):
    """Finds the profile with the highest bins with acceptable SN and Chi values in a directory"""
    bestprof_names = glob.glob(ospj(pipe["run_ops"]["my_dir"], f"*b*{pulsar[1:]}*.bestprof"))
    # throw all of the information from each bestprof into an array
    bin_order = []
    sn_order = []
    chi_order = []
    for info in pipe["folds"]["post"]:
        bin_order.append(info["nbins"])
        sn_order.append(info["sn"])
        chi_order.append(info["chi"])
    bin_order, sn_order, chi_order = zip(
        *sorted(zip(bin_order, sn_order, chi_order)))
    bin_order = bin_order[::-1]
    sn_order = sn_order[::-1]
    chi_order = chi_order[::-1]
    # now find the one with the most bins that meet the sn and chi conditions
    best_i = None
    bin_lim = bin_sampling_limit(pulsar)
    for i in range(len(bin_order)):
        # only consider profiles where the number of bins used is lower than the bin upper limit
        if bin_order[i] <= bin_lim:
            if sn_order[i] >= pipe["run_ops"]["thresh_sn"] and chi_order[i] >= pipe["run_ops"]["thresh_chi"]:
                best_i = i
                break
    if best_i is None:
        raise NoSuitableProfileError(
            "No profiles fit the threshold parameters")
    logger.info("Adequate profile found with {0} bins".format(
        bin_order[best_i]))
    prof_name = glob.glob(ospj(pipe["run_ops"]["my_dir"], f"*b{bin_order[best_i]}_*{pulsar[1:]}*.bestprof"))[0]
    info = bestprof_info(prof_name)
    pipe["source"]["my_DM"] = info["period"]
    pipe["source"]["my_P"] = info["dm"]
    pipe["run_ops"]["vdif"] = bool(glob.glob(ospj(pipe["run_dir"]["my_dir"], ".vdif")))
    if pipe["run_ops"]["vdif"]:
        pipe["source"]["my_bins"] = 1024
    else:
        pipe["source"]["my_bins"] = info["nbins"]
    return info["nbins"]


def upload_and_move_main(pipe):
    fill_pipe_folds(pipe)
    get_best_fold(pipe)
    dep_id = submit_to_db(pipe)
    if "group" not in pipe["run_ops"]["dir"].split("/"):
        dep_id = move_product_dir(pipe, dep_id=dep_id, dep_type="afterok")
    dpp.resubmit_self(pipe, dep_id=dep_id, dep_type="afterok")
