#!/usr/bin/env python3

import os
import glob
import logging
import argparse
import sys
from config_vcs import load_config_file
import psrqpy
import datetime
import numpy as np

import data_processing_pipeline as dpp
from job_submit import submit_slurm
import plotting_toolkit
import find_pulsar_in_obs as fpio
import sn_flux_est as snfe
import stokes_fold
import check_known_pulsars
logger = logging.getLogger(__name__)

#get ATNF db location
try:
    ATNF_LOC = os.environ['PSRCAT_FILE']
except KeyError:
    logger.warn("ATNF database could not be loaded on disk. This may lead to a connection failure")
    ATNF_LOC = None

#load config
comp_config = load_config_file()

#----------------------------------------------------------------------
def move_to_product_dir(pulsar, pointing_dir, obsid):
    """
    Copies the files fitting the input specifications from a prepfold operation to a data product directory
    Parameters:
    -----------
    bins: int
        The number of bins from the prepfold operation
    pulsar: string
        The J name of the pulsar
    pointing_dir: string
        The pointing directory holding the data products
    obsid: int
        The obsid of the observation
    """
    base_dir = comp_config['base_product_dir']
    product_dir = os.path.join(base_dir, obsid, "data_products", pointing_dir)
    all_bins = find_bins_in_dir(pointing_dir)
    bin_limit = bin_sampling_limit(pulsar)
    copy_bins=[]
    if bin_limit<50:
        copy_bins.append(all_bins[0])
    else:
        copy_bins.append( min(all_bins[1:]))
    if bin_limit<100:
        copy_bins.append(50)
    else:
        copy_bins.append(100)

    data_products=[]
    for bins in copy_bins:
        files = glob.glob(os.path.join(pointing_dir, "*{0}*bins*{1}*".format(bins, pulsar[1:])))
        for afile in files:
            data_products.append(afile)

    for product in data_products:
        dpp.copy_data(product, product_dir)

#----------------------------------------------------------------------
def find_fold_times(pulsar, obsid, beg, end, min_z_power=(0.3, 0.1)):
    """
    Finds the fractional time the pulsar is in the beam at some zenith normalized power

    Parameters:
    -----------
    pulsar: string
        Pulsar J name
    obsid: int
        The observation ID
    beg: int
        The beginning of the observation time in gps time
    end: int
        The end of the observation time in gps time
    min_z_power: tuple/list
        OPTIONAL - evaluated the pulsar as 'in the beam' at this normalized zenith power. If None will use [0.3, 0.1] Default: None

    Returns:
    [enter, leave, power]: list
        enter: float
            The time the pulsar enters the beam as a normalized fraction of beg and end. None if pulsar not in beam
        leave: float
            The time the pulsar leaves the beam as a normalized fraction of beg and end. None if pulsar not in beam
        power: float
            The power for which enter and leave are calculated
    """
    if min_z_power is None:
        min_z_power = [0.3, 0.1]
    if not isinstance(min_z_power, list):
        min_z_power = list(min_z_power)

    min_z_power = sorted(min_z_power, reverse=True)
    names_ra_dec = fpio.grab_source_alog(pulsar_list=[pulsar])
    pow_dict, _ = check_known_pulsars.find_pulsars_power(obsid, powers=min_z_power, names_ra_dec=names_ra_dec)
    for power in pow_dict.keys():
        psr_list = pow_dict[power][obsid]
        enter = None
        leave = None
        if psr_list: #if pulsar is in beam for this power coverage
            this_enter, this_leave = snfe.pulsar_beam_coverage(obsid, pulsar, beg=beg, end=end, min_z_power=power)
            if this_enter is not None and this_leave is not None:
                enter = this_enter
                leave = this_leave
                break

    return [enter, leave, power]

#----------------------------------------------------------------------
def add_prepfold_to_commands(run_dir, files="*.fits", pulsar=None, commands=None, prep_ops="", **kwargs):
    """
    Adds prepfold commands to a list

    Parameters:
    -----------
    run_dir: string
        The directory to work in. Typically the pointing directory.
    files: string
        The files to fold on wrt the run directory. Default: '*.fits'
    puslar: string
        OPTIONAL - The J name of the pulsar. If supplied, will use the archived dm and period values for this pulsar. Default: None
    commands: list
        OPTIONAL - A list of commands. Can be empty, this list will be appended to by the function. Default: None
    prep_ops: str
        OPTIONAL - Any prepfold options supplied in string form. Default: ''
    **kwargs:
        Any arguments that can be handed to prepfold. eg. p=0.447
        For any tags, use a blank string. eg. -noclip=""
        Any none values will be ignored
        A suitable dictionary can be generated from make_my_fold_dict()

    Returns:
    --------
    commands: list
        The commands list that was input with the prepfold commands appended
    """
    if commands is None:
        commands = []

    options=""
    for key, val in kwargs.items():
        if val is not None or val is True:
            options += " -{0} {1}".format(key, val)
    options += " {}".format(prep_ops)
    options += " {}".format(files)

    commands.append('cd {0}'.format(run_dir))
    if pulsar:
        commands.append('echo "Folding on known pulsar {}"'.format(pulsar))
        commands.append('psrcat -e {0} > {0}.eph'.format(pulsar))
        commands.append("sed -i '/UNITS           TCB/d' {}.eph".format(pulsar))
        commands.append("prepfold -par {0}.eph {1}".format(pulsar, options))
        commands.append('errorcode=$?')
        commands.append('if [ "$errorcode" != "0" ]; then')
        commands.append('   echo "Folding using the -psr option"')
        commands.append('   prepfold -psr {0} {1}'.format(pulsar, options))
        commands.append('fi')
    else:
        commands.append("prepfold {}".format(options))

    return commands

def make_my_fold_dict(run_params, nbins, initial=None):
    """
    Makes a dictionary that can be sent to add_prepfold_to_commands. The dictionary can be modified or appended to with\
    any other kwargs that prepfold handles. Any key with a blank string is interpreted as a tag.

    Parameters:
    -----------
    run_params: object
        The run_params object defined by data_proces_pipeline
    nbins: int
        The number of bins to fold over
    initial: boolean
        OPTIONAL - True if this is the first run. Will usea wider search range if true. If none, will try to work out. Default: None

    Returns:
    --------
    prepfold_dict: dictionary
        Contains numerous run paramters for prepfold
    """
    if initial is None:
        profs_in_dir = glob.glob(os.path.join(run_params.pointing_dir, "*.pfd.bestprof"))
        if profs_in_dir:
            initial = False
        else:
            initial = True

    #find enter and end times
    enter, leave, _ = find_fold_times(run_params.pulsar, run_params.obsid, run_params.beg, run_params.end, min_z_power=[0.3, 0.1])
    if enter is None or leave is None:
        logger.warn("{} not in beam for given times. Will use entire integration time to fold.".format(run_params.pulsar))
        logger.warn("Used the following parameters:")
        logger.warn("pulsar: {}".format(run_params.pulsar))
        logger.warn("obsid: {}".format(run_params.obsid))
        logger.warn("beg: {}".format(run_params.beg))
        logger.warn("end: {}".format(run_params.end))
        enter=0
        leave=1

    #check for mask
    mask=None
    check_mask = glob.glob(os.path.join(comp_config['base_product_dir'], run_params.obsid, "incoh", "*mask"))
    if check_mask:
        mask = check_mask[0]
    name = "{0}_{1}_bins_{2}".format(run_params.obsid, nbins, run_params.pulsar)

    #make prepfold kwargs
    prepfold_dict = {}
    prepfold_dict["mask"] = mask
    prepfold_dict["o"] = name
    prepfold_dict["start"] = enter
    prepfold_dict["end"] = leave
    prepfold_dict["n"] = nbins
    prepfold_dict["runavg"] = ""
    prepfold_dict["noxwin"] = ""
    prepfold_dict["noclip"] = ""
    prepfold_dict["nsub"] = 256
    prepfold_dict["pstep"] = 1
    prepfold_dict["pdstep"] = 2
    prepfold_dict["dmstep"] = 1
    prepfold_dict["npart"] = 120

    if run_params.dm:
        prepfold_dict["dm"] = run_params.dm
    if run_params.period:
        prepfold_dict["p"] = run_params.period

    #choose the search range basd on initial and nbins
    if initial is None:
        profs_in_dir = glob.glob(os.path.join(run_params.pointing_dir, "*.pfd.bestprof"))
        if profs_in_dir:
            initial = False
        else:
            initial = True
    if nbins == 100:
        prepfold_dict["npfact"] = 1
        prepfold_dict["ndmfact"] = 1
    elif nbins == 50:
        prepfold_dict["npfact"] = 2
        prepfold_dict["ndmfact"] = 1
    elif initial and nbins<300:
        prepfold_dict["npfact"] = 4
        prepfold_dict["ndmfact"] = 3
        prepfold_dict["dmstep"] = 3
        prepfold_dict["npart"] = 40
    else:
        prepfold_dict["npfact"] = 1
        prepfold_dict["ndmfact"] = 1

    if nbins>=300:
        prepfold_dict["nopdsearch"] = ""

    is_bin = is_binary(run_params.pulsar)
    if initial or is_bin:
        prepfold_dict["dm"] = None
        prepfold_dict["p"] = None
        if is_bin:
            logger.info("This is a binary pulsar")
        if initial:
            logger.info("This is the iniital fold")
        logger.info("Fold using pulsar ephemeris: {0}".format(run_params.pulsar))
    else:
        prev_bins = how_many_bins_previous(run_params.pulsar, run_params.pointing_dir)
        info_dict = bestprof_info(glob.glob("{0}/*_{1}*_bins*pfd.bestprof".format(run_params.pointing_dir, prev_bins))[0])
        prepfold_dict["dm"] = info_dict["dm"]
        prepfold_dict["p"] = info_dict["period"]
        logger.info("Will fold using DM: {0} and Period: {1}".format(run_params.dm, run_params.period))
    logger.info("Will fold with npfact: {0} and ndmfact: {1}".format(prepfold_dict["npfact"], prepfold_dict["ndmfact"]))

    return prepfold_dict

#----------------------------------------------------------------------
def submit_prepfold(run_params, nbins, initial=None):
    """
    Submits a prepfold job for the given parameters

    Parameters:
    -----------
    run_params: object
        The run_params object defined by data_processing_pipeline
    nbins: int
        The number of bins to fold on
    initial: boolean
        OPTIONAL - Whether this is the first fold on this pulsar. This affects the search range.\
                    If none, will try to figure it out. Default: None
    """

    if initial is None:
        profs_in_dir = glob.glob(os.path.join(run_params.pointing_dir, "*.pfd.bestprof"))
        if profs_in_dir:
            initial = False
        else:
            initial = True

    if initial or is_binary(run_params.pulsar):
        psr = run_params.pulsar
    else:
        psr = None

    prepfold_dict = make_my_fold_dict(run_params, nbins, initial)
    #make the commands
    commands = []
    commands.append("echo '############### Prepfolding on {} bins ###############'".format(nbins))
    commands = add_prepfold_to_commands(run_params.pointing_dir, files="*.fits", pulsar=psr, commands=commands, **prepfold_dict)

    #Check if prepfold worked:
    commands.append("errorcode=$?")
    commands.append("echo 'errorcode' $errorcode")
    commands.append('if [ "$errorcode" != "0" ]; then')
    commands.append("   echo 'Prepfold operation failure!'")
    commands.append("   exit $errorcode")
    commands.append("fi")

    #binfinder relaunch:
    commands.append("echo '############### Relaunching binfinder script ###############'" )
    bf_relaunch = dpp.binfinder_launch_line(run_params)
    commands.append(bf_relaunch)

    batch_dir = os.path.join(comp_config['base_product_dir'], run_params.obsid, "batch")
    name = "bf_{0}_{1}_{2}_bins".format(run_params.pulsar, run_params.obsid, nbins)

    time = dpp.prepfold_time_alloc(prepfold_dict, run_params.beg, run_params.end)
    if time > 86399.:
        logger.warn("Estimation for prepfold time greater than one day")
        time = 86399
    time = str(datetime.timedelta(seconds = int(time)))

    logger.info("Submitting prepfold and resubmission job:")
    job_id = submit_slurm(name, commands,\
                batch_dir=batch_dir,\
                slurm_kwargs={"time": time},\
                module_list=['mwa_search/{0}'.format(run_params.mwa_search),\
                            'presto/master'],\
                submit=True, vcstools_version="{0}".format(run_params.vcs_tools))

    logger.info("Pointing directory:        {}".format(run_params.pointing_dir))
    logger.info("Pulsar name:               {}".format(run_params.pulsar))
    logger.info("Number of bins to fold on: {}".format(nbins))
    logger.info("Job name:                  {}".format(name))
    logger.info("Time Allocation:           {}".format(time))
    logger.info("Job ID:                    {}".format(job_id))

    return job_id

#----------------------------------------------------------------------
def bestprof_info(filename):
    """
    Finds various information on a .bestprof file

    Parameters:
    filename: string
        The path of the bestprof file
    Returns:
    info_dict: dictionary
        A dictionary consisting of the following:
        obsid: int
            The ID of the observation
        puslar: string
            The J name of the pulsar
        nbins: int
            The number of bins used to fold this profile
        chi: float
            The reduced Chi squared value of the fold
        sn: float
            The signal to noise ratio of the fold
        dm: float
            The pulsar's dispersion measure
        period: float
            The pulsar's period
        period_error: float
            The error in the pulsar's period measurement
    """
    #open the file and read the info into a dictionary
    info_dict = {}
    f = open(filename, "r")
    lines = f.read()
    f.close()
    lines = lines.split("\n")
    #info:
    info_dict["obsid"] = int(lines[0].split()[4].split("_")[0])
    info_dict["pulsar"] = lines[1].split()[3].split("_")[1]
    info_dict["nbins"] = int(lines[9].split()[4])
    info_dict["chi"] = float(lines[12].split()[4])
    info_dict["sn"] = float(lines[13].split()[4][2:])
    info_dict["dm"] = float(lines[14].split()[4])
    info_dict["period"] = float(lines[15].split()[4])/1e3 #in seconds
    info_dict["period_error"] = float(lines[15].split()[6])/1e3
    f.close()
    return info_dict

#----------------------------------------------------------------------
def bin_sampling_limit(pulsar, sampling_rate=1e-4):
    """
    Finds the sampling limit of the input pulsar in units of number of bins

    Parameters:
    -----------
    puslar: string
        The J name of the pulsar to check the sampling limit for
    sampling_rate: float
        OPTIONAL - the sampling rate of the instrument. Default=1e-4 (MWA VCS)

    Returns:
    bin_lim: int
        The highest number of bins that can be logically folded on
    """
    query = psrqpy.QueryATNF(params=["P0"], psrs=[pulsar], loadfromdb=ATNF_LOC).pandas
    period = query["P0"][0]
    bin_lim = int(period/sampling_rate + 1) #the +1 is to round the limit up every time
    logger.debug("Bin limit: {0}".format(bin_lim))
    return bin_lim

def is_binary(pulsar):
    """
    Checks the ATNF database to see if a pulsar is part of a binary system

    Parameters:
    -----------
    pulsar: string
        The J name of the pulsar

    Returns:
    --------
    boolean
        True if the pulsar is a binary. False otherwise
    """
    query = psrqpy.QueryATNF(params=["BINARY"], psrs=[pulsar], loadfromdb=ATNF_LOC).pandas
    if isinstance(query["BINARY"][0], str):
        return True
    else:
        return False

#----------------------------------------------------------------------
def submit_to_db_and_continue(run_params, best_bins):
    """
    Submits the best fold profile to the pulsar database. Will also submit .ppps, .png and .pfd

    Parameters:
    -----------
    run_params: object
        The run_params object defined in data_processing_pipeline
    best_bins: int
    """
    ppps = "*_{}*_bins*.pfd.ps".format(best_bins)
    ppps = glob.glob(os.path.join(run_params.pointing_dir, ppps))[0]
    bestprof = "*_{}*_bins*.pfd.bestprof".format(best_bins)
    bestprof = glob.glob(os.path.join(run_params.pointing_dir, bestprof))[0]
    png = "*_{}*_bins*.png".format(best_bins)
    png = glob.glob(os.path.join(run_params.pointing_dir, png))[0]
    pfd = "*_{}*_bins*.pfd".format(best_bins)
    pfd = glob.glob(os.path.join(run_params.pointing_dir, pfd))[0]

    commands = []
    commands.append("cd {}".format(run_params.pointing_dir))
    commands.append("echo 'Submitting profile to database with {} bins'".format(best_bins))
    commands.append('submit_to_database.py -o {0} --cal_id {1} -p {2} --bestprof {3} --ppps {4}'\
    .format(run_params.obsid, run_params.cal_id, run_params.pulsar, bestprof, ppps))

    #Make a nice plot
    plotting_toolkit.plot_bestprof(os.path.join(run_params.pointing_dir, bestprof),\
                                    out_dir=run_params.pointing_dir)

    bin_lim = bin_sampling_limit(run_params.pulsar)
    if bin_lim>100:
        b_standard=100
    else:
        b_standard=50

    if best_bins != b_standard:
        #do the same for 100/50 bin profiles depending on whether this is an msp or not
        ppps = "*_{}*_bins*.pfd.ps".format(b_standard)
        ppps = glob.glob(os.path.join(run_params.pointing_dir, ppps))[0]
        bestprof = "*_{}*_bins*.pfd.bestprof".format(b_standard)
        bestprof = glob.glob(os.path.join(run_params.pointing_dir, bestprof))[0]
        png = "*_{}*_bins*.png".format(b_standard)
        png = glob.glob(os.path.join(run_params.pointing_dir, png))[0]
        pfd = "*_{}*_bins*.pfd".format(b_standard)
        pfd = glob.glob(os.path.join(run_params.pointing_dir, pfd))[0]

        commands.append("echo 'Submitting profile to database with {} bins'".format(b_standard))
        commands.append('submit_to_database.py -o {0} --cal_id {1} -p {2} --bestprof {3} --ppps {4}'\
                        .format(run_params.obsid, run_params.cal_id, run_params.pulsar, bestprof, ppps))

    if run_params.stokes_dep:
        #submit inverse pfb profile if it exists
        ipfb_archive = os.path.join(run_params.pointing_dir, "{0}_{1}_ipfb_archive.txt".format(run_params.obsid, run_params.pulsar))
        commands.append("echo 'Submitting inverse PFB profile to database'")
        commands.append("submit_to_database.py -o {0} --cal_id {1} -p {2} --ascii {3} --ppps {4} --start {5} --stop {6}"\
                        .format(run_params.obsid, run_params.cal_id, run_params.pulsar, ipfb_archive, ppps,\
                        run_params.beg, run_params.end))

    #Move the pointing directory
    move_loc = os.path.join(comp_config["base_product_dir"], run_params.obsid, "data_products")
    pointing = run_params.pointing_dir.split("/")
    pointing = [i for i in pointing if i != ""]
    new_pointing_dir = os.path.join(move_loc, pointing[-1])
    logger.info("New pointing directory: {}".format(new_pointing_dir))

    #in case the previous command fails. Don't move stuff around
    commands.append("errorcode=$?")
    commands.append("echo 'errorcode' $errorcode")
    commands.append('if [ "$errorcode" != "0" ]; then')
    commands.append("   echo 'Submission Failure!'")
    commands.append("   exit $errorcode")
    commands.append("fi")
    commands.append("echo 'submitted profile to database: {0}'".format(bestprof))
    commands.append("echo 'Moving directory {0} to location {1}'".format(run_params.pointing_dir, move_loc))
    commands.append("mkdir -p {}".format(move_loc))
    if run_params.pulsar[-1].isalpha():
        commands.append("cp -ru {0} {1}".format(run_params.pointing_dir, new_pointing_dir))
    else:
        commands.append("mv {0} {1}".format(run_params.pointing_dir, new_pointing_dir))

    #submit job
    name = "Submit_db_{0}_{1}".format(run_params.pulsar, run_params.obsid)
    batch_dir = os.path.join(comp_config['base_product_dir'], run_params.obsid, "batch")
    logger.info("Submitting submission script for profile: {0}".format(bestprof))
    logger.info("Job name: {}".format(name))

    dep_id = submit_slurm(name, commands,\
                batch_dir=batch_dir,\
                slurm_kwargs={"time": "04:00:00"},\
                depend=run_params.stokes_dep,
                module_list=['mwa_search/{0}'.format(run_params.mwa_search)],\
                submit=True, vcstools_version="{0}".format(run_params.vcs_tools))

    if not run_params.stop:
        #Run stokes fold
        run_params.stokes_bins = best_bins
        launch_line = dpp.stokes_launch_line(run_params, dpp=True, custom_pointing=new_pointing_dir)
        commands=[launch_line]

        name = "dpp_stokes_{0}_{1}".format(run_params.pulsar, run_params.obsid)
        batch_dir = os.path.join(comp_config['base_product_dir'], run_params.obsid, "batch")
        logger.info("Submitting Stokes Fold script")
        logger.info("Job Name: {}".format(name))

        #wait for pfb inversion if it exists
        dep_ids = [dep_id]
        if run_params.stokes_dep:
            dep_ids.append(run_params.stokes_dep)

        submit_slurm(name, commands,\
                    batch_dir=batch_dir,\
                    slurm_kwargs={"time": "00:20:00"},\
                    depend=dep_ids, depend_type="afterany",\
                    module_list=['mwa_search/{0}'.format(run_params.mwa_search)],\
                    submit=True, vcstools_version="{0}".format(run_params.vcs_tools))

#----------------------------------------------------------------------
def get_best_profile_in_dir(pointing_dir, pulsar):
    """
    Finds the profile with the highest bins with acceptable SN and Chi values in a directory

    Parameters:
    -----------
    potining_dir: string
        The poitning directory for the function to look through
    pulsar: string
        The J name of the pulsar

    Returns:
    --------
    prof_name: string
        The name of the best profile found in the directory
    """
    #find all of the relevant bestprof profiles in the pointing directory
    bestprof_names = glob.glob("*bins*{0}*.bestprof".format(pulsar[1:]))
    if len(bestprof_names)==0:
        logger.error("No bestprofs found in directory! Exiting")
        sys.exit(1)

    #throw all of the information from each bestprof into an array
    bin_order = []
    sn_order = []
    chi_order = []
    for prof in bestprof_names:
        prof_info = bestprof_info(prof)
        bin_order.append(prof_info["nbins"])
        sn_order.append(prof_info["sn"])
        chi_order.append(prof_info["chi"])
    bin_order, sn_order, chi_order = zip(*sorted(zip(bin_order, sn_order, chi_order)))
    bin_order = bin_order[::-1]
    sn_order = sn_order[::-1]
    chi_order = chi_order[::-1]

    #now find the one with the most bins that meet the sn and chi conditions
    best_i = None
    bin_lim = bin_sampling_limit(pulsar)
    for i in range(len(bin_order)):
        if bin_order[i]<=bin_lim: #only consider profiles where the number of bins used is lower than the bin upper limit
            if sn_order[i]>10. and chi_order[i]>4.:
                best_i = i
                break

    if best_i is None:
        logger.info("No profiles fit the threshold parameter")
        return None
    else:
        logger.info("Adequate profile found with {0} bins".format(bin_order[best_i]))
        prof_name = glob.glob("*{0}_bins*{1}*.bestprof".format(bin_order[best_i], pulsar[1:]))[0]
        return prof_name

#----------------------------------------------------------------------
def sn_chi_test(bestprof, sn_thresh=10., chi_thresh=4.):
    """
    Checks whether the input signal to noise ratio and chi values are good or bad

    Parameters:
    -----------
    bestprof: string
        The name of the bestprof file to check
    sn_thresh: float
        OPTIONAL - The lower limit on what is an acceptable sn. Default=10.
    ch_thresh: float
        OPTIONAL - The lower limit on what is an acceptable reduced chi squared value. Defaultl=4.

    Returns:
    --------
    test: boolean
        Whether or not the input sn and chi combination is acceptable
    """
    test = False
    info_dict = bestprof_info(bestprof)
    sn = info_dict["sn"]
    chi = info_dict["chi"]
    dm = info_dict["dm"]
    if sn >= sn_thresh and chi >= chi_thresh:
        test = True
    elif sn == 0. and chi >= chi_thresh:
        test = True
    if dm == 0.:
        test = False
        logger.info("This is a satellite")
    return test

#----------------------------------------------------------------------
def submit_multiple_pointings(run_params):
    """
    Submits many pointings on either 50 or 100 bins depending on the pulsar's period. Launches the next step of binfinder

    Parameters:
    -----------
    run_params: object
        The run_params object defined in data_processing_pipeline
    """
    job_ids = []
    #Check number of bins to use
    bin_limit=bin_sampling_limit(run_params.pulsar)
    if bin_limit<100:
        nbins=50
    else:
        nbins=64

    logger.info("Submitting multiple prepfold jobs:")
    logger.info("Pulsar name: {}".format(run_params.pulsar))
    logger.info("Number of bins to fold on: {}".format(nbins))

    #submit a job for each pointing
    for i, pointing in enumerate(run_params.pointing_dir):
        if not glob.glob(os.path.join(pointing, "*.pfd.bestprof")):
            prepfold_dict = make_my_fold_dict(run_params, nbins, True)

            #create slurm job:
            commands =  add_prepfold_to_commands(pointing, files="*.fits", pulsar=run_params.pulsar, **prepfold_dict)
            name = "bf_multi_{0}_{1}_{2}_bins_{3}".format(run_params.pulsar, run_params.obsid, nbins, i)
            batch_dir = os.path.join(comp_config['base_product_dir'], run_params.obsid, "batch")

            time = dpp.prepfold_time_alloc(prepfold_dict, run_params.beg, run_params.end)
            if time > 86399.:
                logger.warn("Estimation for prepfold time greater than one day")
                time = 86399
            time = str(datetime.timedelta(seconds = int(time)))

            logger.info("Submitting pointing: {0}".format(pointing))
            logger.info("Job name: {}".format(name))
            myid = submit_slurm(name, commands,\
                        batch_dir=batch_dir,\
                        slurm_kwargs={"time": time},\
                        module_list=['mwa_search/{0}'.format(run_params.mwa_search),\
                                    'presto/master'],\
                        submit=True, vcstools_version="{0}".format(run_params.vcs_tools))

            logger.info("Pulsar name:               {}".format(run_params.pulsar))
            logger.info("Number of bins to fold on: {}".format(nbins))
            logger.info("Job name:                  {}".format(name))
            logger.info("Time Allocation:           {}".format(time))
            logger.info("Job ID:                    {}".format(myid))
            job_ids.append(myid)
        else:
            logger.info("Pointing {} already has a folded profile. Not Folding".format(pointing))

    launch_line = dpp.binfinder_launch_line(run_params, dpp=False)
    commands = [launch_line]

    name="bf_post_multi_{0}".format(run_params.pulsar)
    batch_dir = os.path.join(comp_config['base_product_dir'], run_params.obsid, "batch")
    logger.info("Submitting post-prepfold binfinder:")
    logger.info("Job name: {}".format(name))
    myid = submit_slurm(name, commands,\
            batch_dir=batch_dir,\
            slurm_kwargs={"time": "00:30:00"},\
            module_list=['mwa_search/{0}'.format(run_params.mwa_search),\
                        "presto/master"],\
            submit=True, depend=job_ids, depend_type="afterany",\
            vcstools_version=run_params.vcs_tools)

#----------------------------------------------------------------------
def find_best_pointing(run_params, nbins):
    """
    Finds the pointing directory with the highest S/N then submits a prepfold job for that pointing

    Parameters:
    -----------
    run_params: object
        The run_params object defined in data_processing_pipeline
    nbins: int
        The number of bins that were folded on that the function will look through
    """

    bestprof_info_list = []
    for pointing in run_params.pointing_dir:
        os.chdir(pointing)
        logger.info("searching directory: {0}".format(pointing))
        bestprof_names = glob.glob("*{0}*_bins*{1}*.bestprof".format(nbins, run_params.pulsar[1:]))
        if len(bestprof_names) == 0:
            logger.warn("{} did not successfully fold".format(pointing))
        else:
            prof_name = bestprof_names[0]
            bestprof_info_list.append(bestprof_info(filename=prof_name))
    if len(bestprof_info_list) == 0:
        logger.error("No pointings have successfully folded! Exiting...")
        sys.exit(1)

    #now we loop through all the info and find the best one
    best_sn = 0.0
    best_i = -1
    for i, info_dict in enumerate(bestprof_info_list):
        if info_dict["chi"]>=4.0 and info_dict["sn"]>best_sn:
            best_sn = info_dict["sn"]
            best_i = i

    if best_i<0 and best_sn<run_params.threshold:
        logger.info("No pulsar found in pointings. Exiting...")
        sys.exit(0)
    else:
        logger.info("Pulsar found in pointings. Running binfinder script on pointing: {0}"\
                    .format(run_params.pointing_dir[best_i]))

    #submit the next iteration
    run_params.set_pointing_dir(run_params.pointing_dir[best_i])
    next_bins = how_many_bins_next(run_params.pulsar, run_params.pointing_dir)
    submit_prepfold(run_params, next_bins)

#----------------------------------------------------------------------
def how_many_bins_next(pulsar, directory):
    """
    Works out how many bins to use for the next fold operation.
    The logic fucntion is based only on the pulsar's period and what folds have been done already.

    Parameters:
    -----------
    pulsar: string
        The J name of the pulsar
    directory: string
        The name of the pointing directory

    Returns:
    --------
    next_bins: int
        The number of bins to fold on next. Returns None if no more folds need to be done
    """
    bins_in_dir = find_bins_in_dir(directory)
    bin_limit = bin_sampling_limit(pulsar)

    #if msp
    if bin_limit<100:
        if 50 not in bins_in_dir:
            return 50
        elif bin_limit not in bins_in_dir:
            return bin_limit
        else:
            return None
    #moderate period pulsar
    elif bin_limit<1024:
        if 64 not in bins_in_dir:
            return 64
        if 100 not in bins_in_dir:
            return 100
        elif bin_limit not in bins_in_dir:
            return bin_limit
        else:
            next_bins =int(min(bins_in_dir[2:])/2)
    #regular period pulsar
    else:
        if 64 not in bins_in_dir:
            return 64
        if 100 not in bins_in_dir:
            return 100
        elif 1024 not in bins_in_dir:
            return 1024
        else:
            next_bins = int(min(bins_in_dir[2:])/2)

    #Check to see if this is an acceptable number
    if next_bins<100:
        next_bins = None

    return next_bins

#----------------------------------------------------------------------
def how_many_bins_previous(pulsar, directory):
    """
    Works backwards from 'how_many_bins_next()' to figure out how many bins would have been folded on previously

    Parameters:
    -----------
    pulsar: string
        The J name of the pulsar
    directory: string
        The name of the pointing directory

    Returns:
    --------
    next_bins: int
        The number of bins to fold on next. Returns None if no folds have been done
    """
    bins_in_dir = find_bins_in_dir(directory)
    bin_limit = bin_sampling_limit(pulsar)
    next_bins = how_many_bins_next(pulsar, directory)

    #next_bins could be none if no more folding needs to be done:
    if not next_bins:
        if bin_limit<100:
            return bin_limit
        elif bin_limit<1024:
            return min(bins_in_dir[2:])
        else:
            return 128

    #if msp
    if bin_limit<100:
        if next_bins==50:
            return None
        elif next_bins == bin_limit:
            return 50

    #moderate period pulsar
    elif bin_limit<1024:
        if next_bins == 64:
            return None
        if next_bins == 100:
            return 64
        elif next_bins == bin_limit:
            return 100
        else:
            return next_bins*2

    #regular period pulsar
    else:
        if next_bins == 64:
            return None
        if next_bins == 100:
            return 64
        elif next_bins == 1024:
            return 100
        else:
            return next_bins*2

    logger.warn("Something has gone wrong when trying to find the previous bins")

#----------------------------------------------------------------------
def find_bins_in_dir(directory):
    """
    Works out what folds have been executed in a given directory already.
    Assumes all files are in the format 'obsid_bins_...pfd.bestprof'

    Parameters:
    -----------
    directory: string
        The directory to look in

    Returns:
    --------
    bins_in_dir: list
        A list of ints of the bin numbers found in the directory in order with duplicates removed
    """
    fold_files=glob.glob(os.path.join(directory, "*.pfd.bestprof"))
    bins_in_dir=[]
    for bestprof in fold_files:
        bins_in_dir.append(int(bestprof.split("/")[-1].split("_")[1]))
    bins_in_dir.sort()
    return sorted(list(set(bins_in_dir)))

#----------------------------------------------------------------------
def work_out_what_to_do(run_params):
    """
    A logic structure that decides what to do next in the binfinder pipeline

    Parameters:
    -----------
    run_params: object
        The run_params object defined by data_proces_pipeline
    """
    os.chdir(run_params.pointing_dir)
    hdr_files = glob.glob("*.hdr".format(run_params.pointing_dir))
    ipfb_archive = glob.glob("*ipfb*.ar".format(run_params.pointing_dir))
    #Multiple pointings?
    if isinstance(run_params.pointing_dir, list):
        logger.debug("More than one pointing to be folded on")
        #Have these directories been folded on?
        any_folded = True
        for pdir in run_params.pointing_dir:
            bins_in_dir = find_bins_in_dir(pdir)
            if len(bins_in_dir) > 0:
                any_folded = True
                break
        if not any_folded:
            #No folds done yet, fold on all pointings
            submit_multiple_pointings(run_params)
        else:
            #Check the pointings and fold on the best one, if any
            find_best_pointing(run_params, nbins=bins_in_dir[0])

    elif isinstance(run_params.pointing_dir, str):
        logger.debug("One pointing to be folded on: {}".format(run_params.pointing_dir))
        #Get some info on where we're at
        bin_limit = bin_sampling_limit(run_params.pulsar)
        bins_in_dir = find_bins_in_dir(run_params.pointing_dir)
        bestprofs_in_dir = glob.glob(os.path.join(run_params.pointing_dir, "*.pfd.bestprof"))
        next_bins = how_many_bins_next(run_params.pulsar, run_params.pointing_dir)

        if hdr_files and not run_params.stokes_dep and not ipfb_archive:
            #.vdif files exist and dspsr job not already submitted
            logger.info(".vdif files available. Will fold on pfb inversion")
            job_id = stokes_fold.submit_inverse_pfb_fold(run_params, stop=True)
            run_params.stokes_dep = job_id

        if len(bins_in_dir)==0:
            logger.debug("No folds have been done. Initial fold")
            submit_prepfold(run_params, next_bins, initial=True)

        elif len(bins_in_dir)==1:
            logger.debug("One fold has been done previously")
            test_dir = os.path.join(run_params.pointing_dir, "*pfd.bestprof")
            test_dir = glob.glob(test_dir)[0]
            test = sn_chi_test(test_dir, sn_thresh=run_params.threshold, chi_thresh=3)
            if test:
                logger.info("Pulsar Detected!")
                logger.debug("Will fold again with {} bins".format(next_bins))
                submit_prepfold(run_params, next_bins, initial=False)
                return
            else:
                logger.info("No pulsar found in initial pointing. Exiting...")
                sys.exit(0)

        elif len(bins_in_dir)==2 and bin_limit>100:
            logger.debug("Two folds have been done previously")
            logger.debug("Will fold again with {} bins".format(next_bins))
            submit_prepfold(run_params, next_bins, initial=False)
            return

        else:
            logger.debug("{} Folds have been done previously".format(len(bins_in_dir)))
            last_fold_bins = how_many_bins_previous(run_params.pulsar, run_params.pointing_dir)

            test_dir = os.path.join(run_params.pointing_dir, "*_{0}*_bins*.bestprof".format(last_fold_bins))
            test_dir = glob.glob(test_dir)[0]
            test = sn_chi_test(test_dir, sn_thresh=run_params.threshold, chi_thresh=3)
            if not test:
                if next_bins is not None:
                    submit_prepfold(run_params, next_bins, initial=False)
                else:
                    logger.info("Minimum bin count hit")
                    if bin_limit<50:
                        sub_bins = bins_in_dir[0] 
                    else:
                        sub_bins = bins_in_dir[2] #the number after the standard 64 and 100 bin trials (usually 128)
                    submit_to_db_and_continue(run_params, sub_bins)
            else:
                info_dict = bestprof_info(test_dir)
                run_params.period = info_dict["period"]
                run_params.dm = info_dict["dm"]
                logger.info("Pulsar detected:   {}".format(run_params.pulsar))
                logger.info("Submitting to database with {0} bins".format(last_fold_bins))
                logger.info("Period:            {}".format(run_params.period))
                logger.info("DM:                {}".format(run_params.dm))
                submit_to_db_and_continue(run_params, last_fold_bins)
                return

#----------------------------------------------------------------------
if __name__ == '__main__':
    #dictionary for choosing log-levels
    loglevels = dict(DEBUG=logging.DEBUG,
                     INFO=logging.INFO,
                     WARNING=logging.WARNING,
                     ERROR = logging.ERROR)

    #Arguments
    parser = argparse.ArgumentParser(description="A script that handles pulsar folding operations",\
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required = parser.add_argument_group("Required Inputs:")
    required.add_argument("-d", "--pointing_dir", action="store", nargs="+", help="Pointing directory(s) that contains the spliced fits files.")
    required.add_argument("-o", "--obsid", type=str, help="The observation ID")
    required.add_argument("-O", "--cal_id", type=str, help="The Obs ID of the calibrator")
    required.add_argument("-p", "--pulsar", type=str, help="The name of the pulsar. eg. J2241-5236")
    required.add_argument("--beg", type=int, help="The beginning of the observation. Will try to find if unsupplied")
    required.add_argument("--end", type=int, help="The end of the observation. Will try to find if unsupplied")

    foldop = parser.add_argument_group("Folding Options:")
    foldop.add_argument("--no_ephem", action="store_true", help="Use this to override the use of an ephemeris for folding the pulsar")
    foldop.add_argument("--dm", type=float, default=None, help="The dispersion measure to fold around")
    foldop.add_argument("--period", type=float, default=None, help="The period to fold around in seconds")
    foldop.add_argument("--prep_ops", type=str, default="", help="Any additional options to use with prepfold in string form\
                        eg. ' -dm 20 -p 0.528'")

    other = parser.add_argument_group("Other Options:")
    other.add_argument("-f", "--freq", type=float, help="The central frequency of the observation in MHz")
    other.add_argument("-t", "--threshold", type=float, default=8.0, help="The signal to noise threshold to stop at. Default = 10.0")
    other.add_argument("-L", "--loglvl", type=str, default="INFO", help="Logger verbosity level. Default: INFO", choices=loglevels.keys())
    other.add_argument("-S", "--stop", action="store_true", help="Use this tag to stop the data processing pipeline when finished binfinding")
    other.add_argument("--dspsr_ops", type=str, default="", help="Any additional options to send to dspsr once binfinder is finished")
    other.add_argument("--stokes_dep", type=int, help="Job ID of a job that needs to be completed before the stokes_fold.py job can begin")
    other.add_argument("--mwa_search", type=str, default="master", help="The version of mwa_search to use.")
    other.add_argument("--vcs_tools", type=str, default="master", help="The version of vcs_tools to use.")

    args = parser.parse_args()
    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    ch.setLevel(loglevels[args.loglvl])
    formatter = logging.Formatter('%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.propagate = False

    #Checking required inputs
    if args.pointing_dir == None:
        logger.error("No pointing directory supplied. Please specify the pointing directory path and rerun")
        sys.exit(1)
    elif args.cal_id == None:
        logger.error("No calibrator ID supplied. Please input a cal ID and rerun")
        sys.exit(1)
    elif args.pulsar == None:
        logger.error("No pulsar name supplied. Please input a pulsar and rerun")
        sys.exit(1)
    if args.end is None or args.beg is None:
        logger.error("Beginning and end times not supplied. Please supply and rerun")
        sys.exit(1)

    rp={}
    rp["pointing_dir"] = args.pointing_dir
    rp["cal_id"] = args.cal_id
    rp["pulsar"] = args.pulsar
    rp["obsid"] = args.obsid
    rp["stop"] = args.stop
    rp["mwa_search"] = args.mwa_search
    rp["vcs_tools"] = args.vcs_tools
    rp["loglvl"] = args.loglvl
    rp["threshold"] = args.threshold
    rp["beg"] = args.beg
    rp["end"] = args.end
    rp["freq"] = args.freq
    rp["dspsr_ops"] = args.dspsr_ops
    rp["prep_ops"] = args.prep_ops
    rp["dm"] = args.dm
    rp["period"] = args.period
    rp["stokes_dep"] = args.stokes_dep
    run_params = dpp.run_params_class(**rp)

    if run_params.freq is None:
        run_params.set_freq_from_metadata(run_params.obsid)

    logger.info("Pointing Dir: {}".format(run_params.pointing_dir))
    work_out_what_to_do(run_params)