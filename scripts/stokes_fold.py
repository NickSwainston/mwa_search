#!/usr/bin/env python3

import logging
import argparse
import sys
import os
from os.path import join as ospj
from os.path import isfile as isfile
import numpy as np
import glob
import config
import psrqpy

from job_submit import submit_slurm
import data_processing_pipeline as dpp
import plotting_toolkit
import binfinder
import rm_synthesis

logger = logging.getLogger(__name__)

#get ATNF db location
try:
    ATNF_LOC = os.environ['PSRCAT_FILE']
except:
    logger.warn("ATNF database could not be loaded on disk. This may lead to a connection failure")
    ATNF_LOC = None

comp_config = config.load_config_file()

#---------------------------------------------------------------
class NotFoundError(Exception):
    """Raise when a value is not found in a file"""
    pass

def plot_everything(run_params):
    """
    Plots polarimetry, RVM fits, chi map and stacked profiles

    Parameters:
    -----------
    run_params: object
        The run_params object from data_processing_pipeline.py
    """
    os.chdir(run_params.pointing_dir)
    filenames_dict = create_filenames(run_params)
    if not isfile(filenames_dict["ascii"]):
        logger.error("Cannot plot without ascii archive")
        return

    #get RM
    if isfile(filenames_dict["rmsynth"]):
        rm_dict     = rm_synthesis.read_rmsynth_out(filenames_dict["rmsynth"])
        rm          = rm_dict["0"]["rm"]
        rm_e        = rm_dict["0"]["rm_e"]
    elif isfile(filenames_dict["rmfit"]):
        rm, rm_e    = find_RM_from_file(filenames_dict["rmfit"])
    if not rm:
        rm, rm_e    = find_RM_from_cat(run_params.pulsar)

    logger.info("Plotting dspsr archive {0} in {1}".format(filenames_dict["ascii"], run_params.pointing_dir))

    logger.info("Plotting polarimetry profile without RVM fit")
    plotting_toolkit.plot_archive_stokes(filenames_dict["ascii"],\
        obsid=run_params.obsid, pulsar=run_params.pulsar, freq=run_params.freq, out_dir=run_params.pointing_dir, rm=rm, rm_e=rm_e)

    #Try too get RVM dictionary + chi map
    try:
        rvm_dict    = read_rvm_fit_file(filenames_dict["rvmfit"])
        chi_map     = read_chi_map(filenames_dict["chimap"])
    except NotFoundError as e:
        rvm_dict    = None
        chi_map     = None

    if rvm_dict: #plot rvm
        if rvm_dict["nbins"]>=5:
            logger.info("Plotting polarimetry profile with RVM fit")
            plotting_toolkit.plot_archive_stokes(filenames_dict["ascii"],\
                obsid=run_params.obsid, pulsar=run_params.pulsar, freq=run_params.freq, out_dir=run_params.pointing_dir,\
                rm=rm, rm_e=rm_e, rvm_fit=rvm_dict)
        else:
            logger.info("Not enough PA points to plot RVM fit")

    if chi_map and rvm_dict: #plot chi_map
        logger.info("Plotting RVM fit chi squared map")
        dof = rvm_dict["dof"]
        chis = np.copy(chi_map["chis"])
        chi_map_name = "{}_RVM_reduced_chi_map.png".format(run_params.file_prefix)
        plotting_toolkit.plot_rvm_chi_map(chi_map["chis"][:], chi_map["alphas"][:], chi_map["betas"][:],\
            name=chi_map_name, my_chi=rvm_dict["redchisq"], my_beta=rvm_dict["beta"], my_alpha=rvm_dict["alpha"])

    #retrieve epn data
    try:
        logger.info("Plotting stacked archival profiles")
        pulsar_dict = plotting_toolkit.get_data_from_epndb(run_params.pulsar)
        #read my data
        pulsar_dict, my_lin_pol = plotting_toolkit.add_ascii_to_dict(pulsar_dict, filenames_dict["ascii"], run_params.freq)
        #ignore any frequencies > 15 000 MHz
        ignore_freqs = []

        for f in pulsar_dict["freq"]:
            if f>15000:
                ignore_freqs.append(f)
        #plot intensity stack
        plotting_toolkit.plot_stack(pulsar_dict["freq"][:], pulsar_dict["Iy"][:], run_params.pulsar,\
                out_dir=run_params.pointing_dir, special_freqs=[run_params.freq], ignore_freqs=ignore_freqs)
        #clip anything without stokes
        pulsar_dict = plotting_toolkit.clip_nopol_epn_data(pulsar_dict)
        #get lin pol - but don't change ours because it could have been generated from psrchive
        lin = plotting_toolkit.lin_pol_from_dict(pulsar_dict)
        for i, f in enumerate(pulsar_dict["freq"]):
            if f==run_params.freq:
                lin[i]=my_lin_pol
        #plot the polarimetry stack
        plotting_toolkit.plot_stack_pol(pulsar_dict["freq"], pulsar_dict["Iy"], lin, pulsar_dict["Vy"], run_params.pulsar,\
                out_dir=run_params.pointing_dir, ignore_freqs=ignore_freqs)
    except plotting_toolkit.NoEPNDBError:
        logger.info("Pulsar not on the EPN database")


def find_RM_from_cat(pulsar):
    """
    Gets rotation measure from prscat query. Returns None if not on catalogue

    Parameters:
    -----------
    pulsar: str
        The J-name of the pulsar

    Returns:
    --------
    rm: float
        The rotation measure
    rm_err: float
        The uncertainty in the rotation measure
    """

    query = psrqpy.QueryATNF(params=["RM"], psrs=[pulsar], loadfromdb=ATNF_LOC).pandas
    rm = query["RM"][0]
    rm_err = query["RM_ERR"][0]

    if np.isnan(rm):
        return None, None
    elif np.isnan(rm_err):
        rm_err = 0.15*rm
    return rm, rm_err

def find_RM_from_file(fname):
    """
    Finds the rotation measure from an input filename as generated by rmfit.
    Returns Nones if rm cold not be generates.

    Parameters:
    -----------
    fname: str
        The path to the file

    Returns:
    --------
    rm: float
        The rotation measure from the file
    rm_err: float
        The uncertainty in the rotation measure
    """
    f = open(fname)
    lines=f.readlines()
    f.close()
    rm=None
    rm_err=None
    for line in lines:
        line = line.split()
        if line[0] == "Best" and line[1] == "RM":
            rm=float(line[3])
            if len(line) >= 5:
                rm_err=float(line[5])
            else:
                logger.warn("Uncertainty for RM not available")
                rm_err=None
            break

    if not rm:
        logger.warn("RM could not be generated from archive file")

    return rm, rm_err

def read_rvm_fit_file(filename):
    """
    Reads a file with the output from psrmodel and returns a dictionary of the results.
    Raises NotFoundError if an expected value is not present in the file

    Parameters:
    -----------
    filename: str
        The pathname of the file with the rvm fit

    Returns:
    --------
    rvm_dict: dictionary
        contains keys:
            nbins: int
                The number of bins used in the rvm fit
            psi_0: float
                The derived psi_0 parameter
            psi_0_e: float
                The uncertainty in psi_0
            beta: float
                The derived beta parameter
            beta_e: float
                The uncertainty in beta
            alpha: float
                The derived alpha parameter
            alpha_e:
                The uncertainty in alpha
            phi_0: float
                The derived phi_0 parameter
            phi_0_e: float
                The uncertainty in phi_0
            redchisq: float
                The reduced chi square of the best fit
            dof: int
                The degrees of freedom of the fit
    """
    keylist = ("nbins", "redchisq", "dof", "psi_0", "psi_0_e", "beta", "beta_e",\
                "alpha", "alpha_e", "phi_0",  "phi_0_e")
    rvm_dict={}
    for key in keylist:
        rvm_dict[key]=None
    f = open(filename)
    lines = f.readlines()
    f.close()
    n_elements = 0
    for i, line in enumerate(lines):
        if line.endswith("bins\n"):
            rvm_dict["nbins"] = int(line.split()[-2])
        elif line[0:6] == "chisq=":
            rvm_dict["redchisq"] = float(line.split()[-1])
            rvm_dict["dof"] = int(line.split()[0].split("=")[-1])
        elif line[0:7] == "psi_0=(":
            psi_0_str = line.split()[0].split("=")[-1].split("(")[-1].split(")")[0].split("+")
            rvm_dict["psi_0"] = float(psi_0_str[0])
            rvm_dict["psi_0_e"] = abs(float(psi_0_str[-1]))
        elif line[0:7] == "beta =(":
            beta_str = line.split()[1].split("(")[-1].split(")")[0].split("+")
            rvm_dict["beta"]  = float(beta_str[0])
        elif line[0:7] == "alpha=(":
            alpha_str = line.split()[0].split("(")[-1].split(")")[0].split("+")
            rvm_dict["alpha"]  = float(alpha_str[0])
        elif line[0:7] == "phi_0=(":
            phi_0_str = line.split()[0].split("(")[-1].split(")")[0].split("+")
            rvm_dict["phi_0"]  = float(phi_0_str[0])
            rvm_dict["phi_0_e"]  = abs(float(phi_0_str[-1]))
        elif line[0:6] == "alpha=":
            n_elements += 1

    rvm_dict["alpha_e"]  = 180/np.sqrt(n_elements)/2
    rvm_dict["beta_e"]  = 180/np.sqrt(n_elements)/2

    for key in keylist:
        if rvm_dict[key] is None:
            raise NotFoundError("{0} not found in file: {1}".format(key, filename))

    return rvm_dict

def read_chi_map(map):
    """
    Reads a chi map of an RVM fit output by psrmodel

    Parameters:
    -----------
    map: str
        The pathname of the map to read

    Returns:
    --------
    chi_map: dictionary
        contains keys:
            alphas: list
                The alpha values in degrees
            betas: list
                The beta values in degrees
            chis: list x list
                The chi values corresponding to the alpha/beta pairs
    """
    f = open(map)
    lines = f.readlines()
    f.close()
    alphas = []
    betas = []
    chis = []
    for i, line in enumerate(lines):
        if not line == "\n":
            chis.append(float(line.split()[2]))
            if len(alphas)==0:
                betas.append(float(line.split()[1]))
        else:
            alphas.append(float(lines[i-1].split()[0]))
    f.close()
    chis = np.reshape(chis, (len(alphas), len(betas)))
    chis = np.transpose(chis)
    chi_map = {"alphas":alphas, "betas":betas, "chis":chis}

    return chi_map

def analytic_pa(phi, alpha, beta, psi_0, phi_0):
    #Inputs should be in radians
    numerator = np.sin(alpha) * np.sin(phi - phi_0)
    denominator = np.sin(beta + alpha) * np.cos(alpha) - np.cos(beta + alpha) * np.sin(alpha) * np.cos(phi - phi_0)
    return np.arctan2(numerator,denominator) + psi_0

def add_rvm_to_commands(run_dir, archive_name, rvmfile="RVM_fit.txt", chimap="chimap.txt", commands=None, res=90):
    """
    Adds the RVM fitting commands to a list

    run_dir: str
        The directory to run the commands in
    archive_name: str
        The name of the archive file to fit
    rvmfile: str
        OPTIONAL - The name of the output RVM fit text file. Default: 'RVM_fit.txt'
    chimap: str
        OPTIONAL - The name of the output chi map file. Default: 'chimap.txt'
    commands: list
        OPTIONAL - A list to append the commands to. Default: None
    res: int
        OPTIONAL - The number of solutions to trial for both alpha and beta. Default: 90

    Returns:
    --------
    commands: list
        A list of commands with the RVM fitting commands appended
    """
    if not commands:
        commands = []
    commands.append("cd {}".format(run_dir))
    commands.append("echo 'Fitting RVM'")
    modelcom = "psrmodel {} -resid -psi-resid -x -use_beta -beta -45:45".format(archive_name)
    modelcom += " -s {0}X{0}".format(res)
    modelcom += " &> {}".format(rvmfile)
    modelcom += " > {}".format(chimap)
    commands.append(modelcom)

    return commands

def add_rm_cor_to_commands(run_dir, RM, archive1, archive2, asciifile, commands=None):
    """
    Adds the commands to correct an archive for the rotation measure of the pulsar

    Parameters:
    -----------
    run_dir: str
        The directory to run the comamnds in
    RM: float
        The rotation measure to correct for
    archive1: str
        The name of the archive file to apply the RM corrections to
    archive2: str
        The name of the RM corrected archive file
    asciifile: str
        The name of the ascii text file to write the archive2 file to
    commands: list
        OPTIONAL - A list to append the commands to. Default: None

    Returns:
    --------
    commands: list
        A list of commands with the rm correction commands appended
    """
    if not commands:
        commands = []
    if not isinstance(RM, float) and not isinstance(RM, int):
        raise ValueError("RM is not a valid value: {}".format(RM))

    #correct for RM
    commands.append("cd {}".format(run_dir))
    commands.append("echo 'Correcting for input rotation measure: {}'".format(RM))
    commands.append("pam -e ar2 -R {0} {1}".format(RM, archive1))

    #Turn the archive into a readable ascii file
    commands.append("echo 'Wiritng result to text file'")
    commands.append("pdv -FTtlZ {0} > {1}".format(archive2, asciifile))

    return commands

def add_rm_fit_to_commands(pulsar, run_dir, archive_name, out_name=None, commands=None):
    """
    Adds the commands to find the rotation measure of an archive

    Parameters:
    -----------
    pulsar: str
        The J name of the pulsar
    run_dir: str
        The directory to run the commands in
    archvie_name: str
        The name of the archive to fit the rotation measure to
    out_name: str
        The name of the output text file. Default: *pulsar*_rmfit.txt
    commands: list
        OPTIONAL - A list to append the commands to. Default: None

    Returns:
    --------
    commands: list
        A list of commands with the rm fitting commands appended
    """
    if not commands:
        commands = []
    if not out_name:
        out_name = ospj(run_dir, "{}_rmfit.txt".format(pulsar))

    commands.append("cd {}".format(run_dir))
    commands.append("echo 'Attempting to find rotation measure.\nOutputting result to {}'".format(out_name))
    commands.append("rmfit {0} -t > {1}".format(archive_name, out_name))

    return commands

def add_rmsynth_to_commands(run_dir, archive_name, label="", write=True, plot=True, keep_QUV=False, commands=None):
    """
    Adds the commands to perform RM synthesis

    Parameters:
    -----------
    run_dir: string
        The location to run the commands
    archive_name: string
        The name of the archive (.ar) to run on
    lebel: string
        A label to apply to the output files. Default: ""
    write: boolean
        OPTIONAL - If True, will write the results of rm_synthesis to file. Default: True
    plot: boolean
        OPTIONAL - If True, will plot the RM synthesis. Default: True
    keep_QUV: boolean
        OPTIONAL - If True, will keep the QUVflux.out file from rmfit. Default: False
    commands: list
        A list of commands to append the rmsynth commands to. Default: None
    """
    if commands is None:
        commands = []

    rms_coms = "rm_synthesis.py"
    rms_coms += " -f {}".format(archive_name)
    if label:
        rms_coms += " --label {}".format(label)
    if write:
        rms_coms += " --write"
    if plot:
        rms_coms += " --plot"
    if keep_QUV:
        rms_coms += " --keep_QUV"
    rms_coms += " --force_single"

    commands.append("cd {}".format(run_dir))
    commands.append("echo 'perfoming RM synthesis'")
    commands.append(rms_coms)

    return commands

def add_pfb_inversion_to_commands(run_dir, pulsar, obsid, archive_out, ascii_out,\
                                nbins=1024, seek=None, total=None, commands=None, tscrunch=100, dm=None, period=None):
    """
    Adds a small dspsr folding pipeline to a list of commands.
    This will fold on each channel using .hdr files, combine all channels and then output the profile as an ascii text file.

    run_dir: string
        The directory to work in. Typically the pointing directory.
    pulsar: string
        The J name of the pulsar
    nbins: int
        OPTIONAL - The number of bins to fold with. Default: 1024
    seek: int
        OPTIONAL - The number of seconds into the obs to start folding. Default: None
    total : int
        OPTIONAL - The total number of seconds of data to fold. Default: None
    commands: list
        OPTIONAL - A list of commands to add the dspsr commands to. Default: None
    tscrunch: int
        OPTIONAL - The number of seconds to timescrunch. Default: 100
    dm: float
        OPTIONAL - The dm to fold around. Default=None
    period: float
        OPTIONAL - The period to fold around. Default=None

    Returns:
    --------
    commands: list
        A list of commands with the dspsr inverse pfb bash commands included
    """

    if commands is None:
        commands = []

    dspsr_coms = "dspsr -U 1000 -A -cont -no_dyn"
    dspsr_coms += " -L {}".format(tscrunch)
    if not dm or not period: #fold with custom dm and period if supplied
        dspsr_coms += " -E {}.eph".format(pulsar)
    dspsr_coms += " -b {}".format(nbins)
    if dm:
        dspsr_coms += " -D {}".format(dm)
    if period:
        dspsr_coms += " -c {}".format(period)
    if seek and total:
        dspsr_coms += " -S {0} -T {1}".format(seek, total)

    psradd_coms = "psradd -R -m time *ar -o {0}".format(archive_out)

    commands.append("cd {0}".format(run_dir))
    commands.append("psrcat -e {0} > {0}.eph".format(pulsar))
    commands.append("echo 'Folding on vdif files'")
    commands.append("j=0")
    commands.append("for i in *.hdr;")
    commands.append("   do {0} -O ipfb_$j $i ;".format(dspsr_coms))
    commands.append("   j=$((j+1))")
    commands.append("done;")
    commands.append("echo 'Combining archives'")
    commands.append(psradd_coms)
    commands.append("echo 'Converting to ascii text'")
    commands.append("pdv -FTtlZ {0} > {1}".format(archive_out, ascii_out))

    return commands

def add_dspsr_fold_to_commands(pulsar, run_dir, nbins,\
                               out_pref=None, commands=None, seek=None, total=None, subint=None, dspsr_ops="", dm=None, period=None):
    """
    Adds a dspsr folding command to a list of commands

    Parameters:
    -----------
    pulsar: str
        The J name of the pulsar
    run_dir: str
        The directory to run the folding operation in
    nbins: int
        The number of bins to fold with
    out_pref: str
        OPTIONAL - The name of the output prefix for the archive file. Default: *pulsar*_archive.ar
    commands: list
        OPTIONAL - A list to add the folding commands to
    seek: int
        OPTIONAL - In seconds, where to begin folding. If none, will not use. Default: None
    total: int
        OPTIONAL - In seconds, the duration to integrate over. If none, will not use. Default: None
    subint: float
        OPTIONAL - In seconds, the length of the sub integrations. Default: 10.
    dspsr_ops: str
        OPTIONAL - A string containing any custom options for dspsr to use. Default: ""
    dm: float
        OPTIONAL - The dm to fold around. Default=None
    period: float
        OPTIONAL - The period to fold around. Default=None

    Returns:
    --------
    Commmands: list
        A list with the dspsr folding commands appended

    """
    if not out_pref:
        out_pref = "{}_archive".format(pulsar)
    if dspsr_ops!='':
        logger.info("dspsr custom options: {}".format(dspsr_ops))
    if not commands:
        commands = []

    dspsr_ops += " {}/*.fits".format(run_dir)
    dspsr_ops += " -O {}".format(ospj(run_dir, out_pref))
    dspsr_ops += " -b {}".format(nbins)
    if dm:
        dspsr_ops += " -D {}".format(dm)
    if period:
        dspsr_ops += " -c {}".format(period)
    if not dm or not period:
        dspsr_ops += " -E {}.eph".format(ospj(run_dir, pulsar))
    if subint:
        dspsr_ops += " -L {}".format(subint)
    if seek:
        dspsr_ops += " -S {}".format(seek)
    if total:
        dspsr_ops += " -T {}".format(total)

    commands.append("cd {}".format(run_dir))
    if not run_params.no_ephem:
        commands.append("psrcat -e {0} > {0}.eph".format(pulsar))
    commands.append("echo 'Running DSPSR folding...'")
    commands.append("dspsr -cont -U 8000 -A -K {}".format(dspsr_ops))

    return commands

def submit_inverse_pfb_fold(run_params, stop=False):
    """
    Submits the inverse pfb folding script and fits RM

    Parameters:
    -----------
    run_params: object
        The run_params object defined in data_processing_pipeline

    Returns:
    --------
    job_id: int
        The id of the submitted job
    """
    #Find beam coverage for known pulsars
    if not run_params.cand:
        enter, leave, _ = binfinder.find_fold_times(run_params.pulsar, run_params.obsid, run_params.beg, run_params.end, min_z_power=[0.3, 0.1])
        obs_int = run_params.end - run_params.beg
        if enter is None or leave is None:
            logger.warn("{} not in beam for given times. Will use entire integration time to fold.".format(run_params.pulsar))
            logger.warn("Used the following parameters:")
            logger.warn("pulsar: {}".format(run_params.pulsar))
            logger.warn("obsid: {}".format(run_params.obsid))
            logger.warn("beg: {}".format(run_params.beg))
            logger.warn("end: {}".format(run_params.end))
            enter_sec = 0
            duration = obs_int
        else:
            duration = (leave - enter) * obs_int
            enter_sec = enter * obs_int
            logger.info("{0} enters beam at {1} and leaves at {2}".format(run_params.pulsar, enter, leave))
            logger.info("Integration time: {}".format(duration))
    else:
        enter_sec = None
        duration = None

    #pfb inversion
    filenames_dict = create_filenames(run_params)
    commands = add_pfb_inversion_to_commands(run_params.pointing_dir, run_params.pulsar, run_params.obsid, filenames_dict["archive1"], filenames_dict["ascii"],\
               seek=enter_sec, total=duration, tscrunch=duration, dm=run_params.dm, period=run_params.period)

    #launch RM fitting
    commands = add_rm_fit_to_commands(run_params.pulsar, run_params.pointing_dir, filenames_dict["archive1"], out_name=filenames_dict["rmfit"], commands=commands)

    #launch RM synthesis
    commands = add_rmsynth_to_commands(run_params.pointing_dir, filenames_dict["archive1"], write=True, plot=True, keep_QUV=False, label=run_params.file_prefix, commands=commands)

    if not stop:
        #Relaunch stokes_fold.py
        launch_line = dpp.stokes_launch_line(run_params)
        commands.append(launch_line)
    elif not run_params.stop:
        launch_line = dpp.stokes_launch_line(run_params)
        commands.append(launch_line)

    batch_dir = ospj(comp_config['base_product_dir'], run_params.obsid, "batch")
    name = "inverse_pfb_{0}".format(run_params.file_prefix)

    logger.info("Submitting inverse pfb job:")
    logger.info("Pointing directory: {}".format(run_params.pointing_dir))
    logger.info("Pulsar name: {}".format(run_params.pulsar))
    logger.info("Job name: {}".format(name))
    job_id = submit_slurm(name, commands,\
                        batch_dir=batch_dir,\
                        slurm_kwargs={"time": "10:00:00"},\
                        module_list=['mwa_search/{0}'.format(run_params.mwa_search),\
                                    "dspsr", "psrchive"],\
                        submit=True, vcstools_version="{0}".format(run_params.vcs_tools))

    return job_id

def submit_dspsr_rmfit(run_params):
    """
    Runs dspsr on fits files and relaunches the stokes fold script

    Parameters:
    -----------
    run_params: object
        The run_params object from data_processing_pipeline.py
    """
    if not run_params.cand:
        enter, leave, _ = binfinder.find_fold_times\
                        (run_params.pulsar, run_params.obsid, run_params.beg, run_params.end, min_z_power=[0.3, 0.1])
        obs_int = run_params.end - run_params.beg
        if enter is None or leave is None:
            logger.warn("{} not in beam for given times. Will use entire integration time to fold.".format(run_params.pulsar))
            logger.warn("Used the following parameters:")
            logger.warn("pulsar: {}".format(run_params.pulsar))
            logger.warn("obsid: {}".format(run_params.obsid))
            logger.warn("beg: {}".format(run_params.beg))
            logger.warn("end: {}".format(run_params.end))
            enter_sec = 0
            duration = obs_int
        else:
            duration = (leave - enter) * obs_int
            enter_sec = enter * obs_int
            logger.info("{0} enters beam at {1} and leaves at {2}".format(run_params.pulsar, enter, leave))
            logger.info("Integration time: {}".format(duration))
    else:
        enter_sec = None
        duration = None

    filenames_dict = create_filenames(run_params)
    #dspsr command
    commands = add_dspsr_fold_to_commands(run_params.pulsar, run_params.pointing_dir, run_params.stokes_bins, out_pref=run_params.file_prefix,\
                                        seek=enter_sec, total=duration, subint=run_params.subint, dspsr_ops=run_params.dspsr_ops,\
                                        dm=run_params.dm, period=run_params.period)
    #rmfit command
    commands = add_rm_fit_to_commands(run_params.pulsar, run_params.pointing_dir, filenames_dict["archive1"], out_name=filenames_dict["rmfit"], commands=commands)

    #rmsynth command
    commands = add_rmsynth_to_commands(run_params.pointing_dir, filenames_dict["archive1"], label=run_params.file_prefix, write=True, plot=True, keep_QUV=False, commands=commands)

    #rerun the script
    if not run_params.stop:
        launch_line = dpp.stokes_launch_line(run_params)
        commands.append(launch_line)

    name = "DSPSR_RMfit_{0}".format(run_params.file_prefix)
    batch_dir = "{}".format(ospj(comp_config['base_product_dir'], run_params.obsid, "batch"))
    job_id = submit_slurm(name, commands,\
                        batch_dir=batch_dir,\
                        slurm_kwargs={"time": "08:00:00"},\
                        module_list=["mwa_search/{0}".format(run_params.mwa_search),\
                                    "dspsr/master", "psrchive/master"],\
                        submit=True, vcstools_version=run_params.vcs_tools, mem="")

    logger.info("Job submitted using\n\
                pointing directory:         {0}\n\
                pulsar:                     {1}"\
                .format(run_params.pointing_dir, run_params.pulsar))

    return job_id

def submit_rm_cor_rvm(run_params):
    """
    Runs the RM correction on the dspsr archive and writes the result to a text file.
    Relaunches the stokes_fold script afterwards for plotting

    Parameters:
    -----------
    run_params: object
        The run_params object from data_processing_pipeline.py
    """
    os.chdir(run_params.pointing_dir)
    filenames_dict = create_filenames(run_params)
    #Look for files in pointing dir
    is_rmsynth  = isfile(filenames_dict["rmsynth"])
    is_rmfit    = isfile(filenames_dict["rmfit"])

    #Correct for RM
    if is_rmsynth:
        logger.info("Using RM synthesis result for correction")
        rm_dict     = rm_synthesis.read_rmsynth_out(filenames_dict["rmsynth"])
        RM          = rm_dict["0"]["rm"]
    elif is_rmfit:
        logger.info("Using rmfit result for correction")
        RM          = find_RM_from_file(filenames_dict["rmfit"])[0]
    if not RM:
        RM          = find_RM_from_cat(run_params.pulsar)[0]
    run_params.RM   = RM
    commands        = add_rm_cor_to_commands(run_params.pointing_dir, run_params.RM, filenames_dict["archive1"], filenames_dict["archive2"], filenames_dict["ascii"])

    #RVM fitting
    commands = add_rvm_to_commands(run_params.pointing_dir, filenames_dict["archive2"], rvmfile=filenames_dict["rvmfit"], chimap=filenames_dict["chimap"],\
               commands=commands, res=run_params.rvmres)

    #relaunch
    if not run_params.stop:
        launch_line = dpp.stokes_launch_line(run_params)
        commands.append(launch_line)

    job_name = "ipfb_RMcor_RVM_{}".format(run_params.file_prefix)
    if not run_params.ipfb:
        job_name = job_name[5:]
    batch_dir = "{0}{1}/batch/".format(comp_config['base_product_dir'], run_params.obsid)
    job_id = submit_slurm(job_name, commands,\
                        batch_dir=batch_dir,\
                        slurm_kwargs={"time": "06:00:00"},\
                        module_list=["mwa_search/{0}".format(run_params.mwa_search),
                                    "psrchive/master"],\
                        submit=True, vcstools_version=run_params.vcs_tools, mem="")

    return job_id

def create_filenames(run_params):
    """
    Creates a dictionary of filenames to use. This is here to ensure the filenames are identical between fucntions and runs
    Note: these are file names only and do not contain the path
    
    Parameters:
    -----------
    run_params: object
        The run_params object defined by data_procesing_pipeline
    
    Returns:
    --------
    filenames_dict: dictionary
        contains the following keys with their respective file names:
        archive1: str
        archive2: str
        ascii: str
        rvmfit: str
        rmfit: str
        rmsytnh: str
        rvmfit: str
        chimap: str
    """
    filenames_dict={}
    if run_params.ipfb:
        filenames_dict["archive1"]  = "{0}_ipfb.ar".format(run_params.file_prefix)
        filenames_dict["archive2"]  = "{0}_ipfb.ar2".format(run_params.file_prefix)
        filenames_dict["ascii"]     = "{0}_ipfb.txt".format(run_params.file_prefix)
        filenames_dict["rvmfit"]    = "{0}_ipfb_RVM_fit.txt".format(run_params.file_prefix)
    else:
        filenames_dict["archive1"]  = "{0}.ar".format(run_params.file_prefix)
        filenames_dict["archive2"]  = "{0}.ar2".format(run_params.file_prefix)
        filenames_dict["ascii"]     = "{0}.txt".format(run_params.file_prefix)
        filenames_dict["rvmfit"]    = "{0}_RVM_fit.txt".format(run_params.file_prefix)
    filenames_dict["rmfit"]         = "{}_rmfit.txt".format(run_params.file_prefix)
    filenames_dict["rmsynth"]       = "{}_RMsynthesis.txt".format(run_params.file_prefix)
    filenames_dict["rvmfit"]        = "{}_RVMfit.txt".format(run_params.file_prefix)   
    filenames_dict["chimap"]        = "{}_chi_map.txt".format(run_params.file_prefix)

    return filenames_dict

def work_out_what_to_do(run_params):
    """
    A logic structure that decides what to do next in the stokes_fold pipeline

    Parameters:
    -----------
    run_params: object
        The run_params object defined by data_procesing_pipeline
    """
    os.chdir(run_params.pointing_dir)
    fits_files_in_dir   = glob.glob("*.fits")
    hdr_files_in_dir    = glob.glob("*.hdr")
    filenames_dict      = create_filenames(run_params)
    is_ar_rm            = isfile(filenames_dict["rmfit"]) and isfile(filenames_dict["archive1"])
    is_ar2_rvm          = isfile(filenames_dict["rvmfit"]) and isfile(filenames_dict["archive2"])

    #Main logic structure
    if hdr_files_in_dir or fits_files_in_dir:
        if not is_ar_rm:
            #Submit the fold and rmfit job
            if run_params.ipfb:
                submit_inverse_pfb_fold(run_params)
            else:
                submit_dspsr_rmfit(run_params)
            return

        elif is_ar_rm and not is_ar2_rvm:
            #Submit the rm correction and RVM fitting job
            submit_rm_cor_rvm(run_params)
            return

        elif is_ar2_rvm:
            #plot
            plot_everything(run_params)
            return

        else:
            logger.error("Something has gone wrong trying to process the .vdif or .fits files :/")
            logger.info("Files found in pointing dir: {}".format(run_params.pointing_dir))
            logger.info("RM fit files: {}".format(rm_fit_files))
            logger.info("RVM fit files: {}".format(rvm_fit_files))
            logger.info("archive files: {}".format(ar_files))
            logger.info("RM corrected archive files: {}".format(ar2_files))
            return

    else:
        logger.error("No valid files in directory: {}".format(ospj(run_params.pointing_dir, "*.fits")))
        logger.debug("glob output: {}".format(fits_files_in_dir))
        return

if __name__ == '__main__':

    loglevels = dict(DEBUG=logging.DEBUG,
                     INFO=logging.INFO,
                     WARNING=logging.WARNING,
                     ERROR = logging.ERROR)

    parser = argparse.ArgumentParser(description="""Folds across stokes IQUV and attempts to find the RM""",\
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    foldop = parser.add_argument_group("Folding Options:")
    foldop.add_argument("-d", "--pointing_dir", type=str, help="Pointing directory that contains the fits files")
    foldop.add_argument("-p", "--pulsar", type=str, default=None, help="The J name of the pulsar.")
    foldop.add_argument("-b", "--nbins", type=int, default=0, help="The number of bins to fold the profile with")
    foldop.add_argument("-s", "--subint", type=float, default=None, help="The length of the integrations in seconds")
    foldop.add_argument("-o", "--obsid", type=str, help="The obsid of the observation")
    foldop.add_argument("--beg", type=int, help="The beginning of the observation time in gps time")
    foldop.add_argument("--end", type=int, help="The end of the observation time in gps time")
    foldop.add_argument("--dm", type=float, default=None, help="The dispersion measure to fold around")
    foldop.add_argument("--period", type=float, default=None, help="The period to fold around in milliseconds")
    foldop.add_argument("--dspsr_ops", type=str, default="", help="Provide as a string in quotes any dspsr command you would like to use for folding.\
                        eg: '-D 50.0 -c 0.50625'")
    foldop.add_argument("--no_ephem", action="store_true", help="Use this tag to override the use of the epehemeris")
    foldop.add_argument("--cand", action="store_true", help="Use this tag if this is not a known pulsar")

    rvmop = parser.add_argument_group("RVM Fitting Options:")
    rvmop.add_argument("--rvmres", type=int, default=90, help="The number of degree samples to try for alpha and beta.")

    otherop = parser.add_argument_group("Other Options:")
    otherop.add_argument("-L", "--loglvl", type=str, default="INFO", help="Logger verbosity level. Default: INFO", choices=loglevels.keys())
    otherop.add_argument("--vcs_tools", type=str, default="master", help="The version of vcstools to use. Default: master")
    otherop.add_argument("--mwa_search", type=str, default="master", help="The version of mwa_search to use. Default: master")
    otherop.add_argument("-S", "--stop", action="store_true", help="Use this tag to stop processing data after the chose mode has finished its intended purpose")
    otherop.add_argument("-f", "--freq", type=float, help="The central frequency of the observation in MHz")

    args = parser.parse_args()

    logger.setLevel(loglevels[args.loglvl])
    ch = logging.StreamHandler()
    ch.setLevel(loglevels[args.loglvl])
    formatter = logging.Formatter('%(asctime)s  %(filename)s  %(name)s  %(lineno)-4d  %(levelname)-9s :: %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.propagate = False


    if args.cand:
        if not args.no_ephem:
            logger.warn("no_ephem tag needs to be used for candidate folds, but is turned off. Overriding.")
            args.no_ephem = True
        rounded_period = round(args.period, 5)
        rounded_dm = round(args.dm, 4)
        args.pulsar="cand".format(rounded_period, rounded_dm)

    else:
        if not args.beg or not args.end:
            logger.error("Beginning/end times not supplied. Please run again and specify times")
            sys.exit(1)

    if not args.pointing_dir:
        logger.error("Pointing directory not supplied. Please run again and specify a pointing directory")
        sys.exit(1)
    if not args.obsid:
        logger.error("Obsid not supplied. Please run again and specify an observation ID")
        sys.exit(1)
    if not args.pulsar:
        logger.error("Pulsar name not supplied. Please run again and specify pulsar name")
        sys.exit(1)
    if args.no_ephem and (args.period is None or args.dm is None):
        logger.warn("Period/DM not explicitly supplied and no ephemeris being used")

    os.chdir(args.pointing_dir)
    ipfb = bool(glob.glob("*hdr"))

    rp={}
    rp["pointing_dir"]      = args.pointing_dir
    rp["pulsar"]            = args.pulsar
    rp["obsid"]             = args.obsid
    rp["stop"]              = args.stop
    rp["mwa_search"]        = args.mwa_search
    rp["vcs_tools"]         = args.vcs_tools
    rp["loglvl"]            = args.loglvl
    rp["stokes_bins"]       = args.nbins
    rp["subint"]            = args.subint
    rp["beg"]               = args.beg
    rp["end"]               = args.end
    rp["freq"]              = args.freq
    rp["dspsr_ops"]         = args.dspsr_ops
    rp["no_ephem"]          = args.no_ephem
    rp["dm"]                = args.dm
    rp["period"]            = args.period
    rp["cand"]              = args.cand
    rp["rvmres"]            = args.rvmres
    rp["ipfb"]              = ipfb
    run_params = dpp.run_params_class(**rp)

    work_out_what_to_do(run_params)