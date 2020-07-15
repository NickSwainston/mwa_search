#!/usr/bin/env python3

from binfinder import NoSuitableProfileError, bestprof_info
from plotting_toolkit import plot_bestprof

def submit_to_db(pipe):
    """Submits the best fold profile to the pulsar database. Will also submit .ppps, .png and .pfd"""
    bin_list = pipe["folds"]["init"].keys().append(pipe["folds"]["post"].keys())
    commands = []
    commands.append(f"cd {pipe['run_ops']['my_dir']}")
    for bin_count in bin_list:
        ppps = f"*_b{bin_count}*_{pipe['source']['name'][1:]}*.pfd.ps"
        ppps = glob.glob(os.path.join(pipe["run_ops"]["my_dir"], ppps))[0]
        bestprof = f"*_b{bin_count}*_{pipe['source']['name'][1:]}*".format(best_bins)
        bestprof = glob.glob(os.path.join(pipe["run_ops"]["my_dir"], bestprof))[0]
        png = f"*_b{bin_count}*_{pipe['source']['name'][1:]}*.png".format(best_bins)
        png = glob.glob(os.path.join(pipe["run_ops"]["my_dir"], png))[0]
        pfd = f"*_b{bin_count}*_{pipe['source']['name'][1:]}*.pfd".format(best_bins)
        pfd = glob.glob(os.path.join(pipe["run_ops"]["my_dir"], pfd))[0]
        #Make a nice plot
        plotting_toolkit.plot_bestprof(os.path.join(pipe["run_ops"]["my_dir"], bestprof, out_dir=pipe["run_ops"]["my_dir"])
        commands.append(f"echo 'Submitting profile to database with {bin_count} bins'")
        commands.append(f"submit_to_database.py -o {pipe['obs']['id']} --cal_id {pipe['obs']['cal_id']} -p {pipe['source']['name']} --bestprof {bestprof} --ppps {ppps}"

    #submit job
    name = f"Submit_db_{pipe['source']['name']}_{pipe['obs']['id']}"
    batch_dir = os.path.join(comp_config['base_product_dir'], pipe['obs']['id'], "batch")
    dep_id = submit_slurm(name, commands,
                batch_dir=batch_dir,
                slurm_kwargs={"time": "04:00:00"},
                depend=run_params.stokes_dep,
                module_list=[f"mwa_search/{pipe['run_ops']['mwa_search']}", f"mwa_search/{pipe['run_ops']['vcstools']}"],
                submit=True)

    logger.info(f"Submission script on queue for profile: {bestprof}")
    logger.info(f"Job name: {name}")
    logger.info(f"Job ID: {dep_id}")

    return dep_id

def move_product_dir(pipe)

    if run_params.pulsar[-1].isalpha():
        commands.append("cp -ru {0} {1}".format(run_params.pointing_dir, new_pointing_dir))
    else:
        commands.append("mv {0} {1}".format(run_params.pointing_dir, new_pointing_dir))


def fill_pipe_folds(pipe):
    """Fills the pipe with thje bestprof information from all available folds"""
    for bins in pipe["folds"]["init"]:
        bestprof = glob.glob(f"*b{bins}**{pipe['source']['name'][1:]}*.bestprof")[0]
        info = bestprof_info(filename=bestprof)
        pipe["folds"]["init"][bins] = info
    for bins in pipe["folds"]["post"]:
        bestprof = glob.glob(f"*b{bins}**{pipe['source']['name'][1:]}*.bestprof")[0]
        info = bestprof_info(filename=bestprof)
        pipe["folds"]["post"][bins] = info

def get_best_fold(pipe):
    """Finds the profile with the highest bins with acceptable SN and Chi values in a directory"""
    bestprof_names = glob.glob(f"*b*{pulsar[1:]}*.bestprof")
    #throw all of the information from each bestprof into an array
    bin_order = []
    sn_order = []
    chi_order = []
    for info in pipe["folds"]["post"]:
        bin_order.append(info["nbins"])
        sn_order.append(info["sn"])
        chi_order.append(info["chi"])
    bin_order, sn_order, chi_order = zip(*sorted(zip(bin_order, sn_order, chi_order)))
    bin_order = bin_order[::-1]
    sn_order = sn_order[::-1]
    chi_order = chi_order[::-1]
    #now find the one with the most bins that meet the sn and chi conditions
    best_i = None
    bin_lim = bin_sampling_limit(pulsar)
    for i in range(len(bin_order)):
        if bin_order[i]<=bin_lim: #only consider profiles where the number of bins used is lower than the bin upper limit
            if sn_order[i]>=pipe["run_ops"]["thresh_sn"] and chi_order[i]>=pipe["run_ops"]["thresh_chi"]:
                best_i = i
                break
    if best_i is None:
        raise NoSuitableProfileError("No profiles fit the threshold parameters")
    logger.info("Adequate profile found with {0} bins".format(bin_order[best_i]))
    prof_name = glob.glob(f"*b{bin_order[best_i]}_*{pulsar[1:]}*.bestprof")[0]
    info = bestprof_info(prof_name)
    pipe["source"]["my_DM"] = info["period"]
    pipe["source"]["my_P"] = info["dm"]

    return bin_order[best_i]

    def submit_move_main(pipe):
        fill_pipe_folds(pipe)
        get_best_fold(pipe)
