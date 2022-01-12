import logging
import sys
from os.path import join
from glob import glob

from vcstools.prof_utils import subprocess_pdv, get_from_ascii, ProfileLengthError, NoFitError
from vcstools.gfit import gfit
from vcstools.job_submit import submit_slurm
from vcstools.config import load_config_file
from dpp.helper_bestprof import bestprof_fit

comp_config = load_config_file()
logger = logging.getLogger(__name__)


def archive_fit(cfg, archive_path, cliptype="verbose"):
    """Fits a profile to the supplied archive and adds it to cfg. Cliptype options found in prof_utils.py"""
    # Get the profile
    subprocess_pdv(archive_path, outfile=cfg["files"]["archive_ascii"])
    profile, _ = get_from_ascii(cfg["files"]["archive_ascii"])
    # Gaussian fit
    g_fitter = gfit(profile, plot_name=cfg["files"]["gfit_plot"])
    g_fitter.auto_gfit()
    fit = g_fitter.fit_dict
    g_fitter.plot_fit()
    # Find the longest component
    longest_comp = 0
    for comp_name in fit["comp_idx"].keys():
        comp = fit["comp_idx"][comp_name]
        if len(comp) > longest_comp:
            longest_comp = len(comp)
            cfg["source"]["my_component"] = comp_name
    # Add the fit to cfg
    cfg["source"]["gfit"] = fit


def fits_to_archive(fits_dir, ar_name, bins, dm, period, out_dir, memory=4000, total=1.0, seek=0, vdif=False, container="/pawsey/mwa/singularity/dspsr/dspsr.sif"):
    """Returns bash commands to fold on a fits file using dspsr"""
    # Normally I'd just use absolute filepaths but psrchive runs into issues if the filename
    # is too long. So instead we cd into the fits dir and process there, then move everything
    # back to the out_dir
    commands = [f"cd {fits_dir}"]
    container_launch = (f"singularity exec -e {container}") # Open the container
    dspsr_cmd = f"{container_launch} dspsr  -A -K -cont"
    dspsr_cmd += f" -U {memory}"
    dspsr_cmd += f" -b {bins}"
    dspsr_cmd += f" -c {period}"
    dspsr_cmd += f" -D {dm}"
    dspsr_cmd += f" -S {seek}"
    dspsr_cmd += f" -T {total}"
    dspsr_cmd += f" -L {total}" # Timescrunch the whole obs
    if vdif: # TODO: fix this because it can't deal with long filenames
        commands.append("j=0")
        commands.append(f"for i in *.hdr;")
        commands.append(f"   do {dspsr_cmd} -O ipfb_$j $i ;")
        commands.append("   j=$((j+1))")
        commands.append("done;")
        commands.append(f"{container_launch} psradd -R -m time ipfb*.ar -o {ar_name}.ar")
    else:
        dspsr_cmd += f" -O {ar_name}"
        dspsr_cmd += f" *.fits"
        commands.append(dspsr_cmd)
    commands.append(f"mv {ar_name}.ar {out_dir}")
    commands.append("cd -")
    commands = "\n".join(commands) # Convert list to string
    return commands


def archive_to_fits(ar_file, extension="fits", container="/pawsey/mwa/singularity/dspsr/dspsr.sif"):
    """Returns bash commands to turn an arhive file to a fits file"""
    container_launch = f"singularity exec -e {container}"
    pam_cmd = f"{container_launch} pam -a PSRFITS"
    pam_cmd += f" -e {extension}"
    pam_cmd += f" {ar_file}"
    return pam_cmd


def remove_baseline(cfg):
    """Removes baseline RFI from a fits file"""
    my_comp = cfg["source"]["my_component"]
    on_pulse = cfg["source"]["gfit"]["comp_idx"][my_comp]
    debase_cmd = "pmod -debase"
    debase_cmd += f" -onpulse '{on_pulse[0]} {on_pulse[-1]}'" # Measured in bins
    debase_cmd += " -ext debase.gg" # Aligns with what's written in the config - Don't touch unless you know what you're doing
    debase_cmd += f" -device /null" # Don't make an image (it's not a very interesing one anyway)
    debase_cmd += f" {cfg['files']['converted_fits']}"
    return debase_cmd


def ppp_archive_creation(cfg, depends_on=None, depend_type="afterany"):
    """Makes commands for converting ot archive file then converting back to fits"""
    fits_dir = join(cfg["files"]["psr_dir"], cfg["source"]["my_pointing"])
    bins = cfg["source"]["my_bins"] if not cfg["run_ops"]["vdif"] else 1024
    dm = cfg["source"]["my_DM"]
    period = cfg["source"]["my_P"]
    total = cfg["source"]["total"]
    seek = cfg["source"]["seek"]
    # Change to working directory
    commands = [f"cd {cfg['files']['psr_dir']}"]
    # Check for vdif
    vdif_hdrs = join(cfg["source"]["my_pointing"], "*.hdr")
    cfg["run_ops"]["vdif"] = bool(glob(vdif_hdrs))
    # Add folds to commands
    psrchive_container = comp_config['prschive_container']
    archive_base = cfg["files"]["archive"].split(".ar")[0] # Archive without .ar extension
    commands.append(fits_to_archive(fits_dir, archive_base, bins, dm, period, cfg["files"]["psr_dir"],
                    total=total, seek=seek, vdif=cfg["run_ops"]["vdif"], container=psrchive_container))
    # Add ar -> fits conversion to commands
    commands.append(archive_to_fits(cfg["files"]["archive"], container=psrchive_container))
    #Submit_job
    name = f"to_archive_{cfg['source']['name']}_{cfg['obs']['id']}"
    slurm_kwargs = {"time":"12:00:00"} # dspsr folding can take some time
    modules = ["singularity"]
    mem=32768
    jid = submit_slurm(name, commands,
        slurm_kwargs=slurm_kwargs, module_list=modules, mem=mem, batch_dir=cfg["files"]["batch_dir"], depend=depends_on,
        depend_type=depend_type, vcstools_version=cfg["run_ops"]["vcstools"], submit=True)
    logger.info(f"Submitted archive/fits creation job: {name}")
    logger.info(f"job ID: {jid}")
    if depends_on:
        logger.info(f"Job depends on job id(s): {depends_on}")
    return jid, name


def ppp_baseline_removal(cfg, depends_on=None, depend_type="afterany"):
    """Submits a job that removes baseline RFI"""
    # Change to working directory
    commands = [f"cd {cfg['files']['psr_dir']}"]
    # Fit the profile with gaussian
    try:
        archive_fit(cfg, cfg["files"]["archive"], cliptype="verbose")
    except (ProfileLengthError, NoFitError) as e:
        ex_type, _, _ = sys.exc_info() # Raise different messages for different errors
        if ex_type == NoFitError:
            raise NoFitError("A Gaussian fit to this profile could not be made. This profile is likely too noisy")
        elif ex_type == ProfileLengthError:
            raise ProfileLengthError("No VDIF files available and profile is not long enough to fit Gaussian")
    # Add the baseline removal commands
    commands.append(remove_baseline(cfg))
    #Submit_job
    name = f"debase_{cfg['source']['name']}_{cfg['obs']['id']}"
    slurm_kwargs = {"time":"01:00:00"}
    modules = ["singularity", "psrsalsa"]
    mem=32768
    jid = submit_slurm(name, commands,
        slurm_kwargs=slurm_kwargs, module_list=modules, mem=mem, batch_dir=cfg["files"]["batch_dir"], depend=depends_on,
        depend_type=depend_type, vcstools_version=cfg["run_ops"]["vcstools"], submit=True)
    logger.info(f"Submitted archive/fits creation job: {name}")
    logger.info(f"job ID: {jid}")
    if depends_on:
        logger.info(f"Job depends on job id(s): {depends_on}")
    cfg["completed"]["debase"] = True
    return jid, name