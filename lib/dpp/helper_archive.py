import logging
from os.path import join

from vcstools.job_submit import submit_slurm
from vcstools.config import load_config_file
from dpp.helper_bestprof import bestprof_fit

comp_config = load_config_file()
logger = logging.getLogger(__name__)


def fits_to_archive(fits_dir, ar_name, bins, dm, period, memory=4000, total=1.0, seek=0, vdif=False, container="/pawsey/mwa/singularity/dspsr/dspsr.sif"):
    """Returns bash commands to fold on a fits file usinf dspsr"""
    commands = []
    container_launch = (f"singularity exec -e {container}") # Open the container
    dspsr_cmd = f"{container_launch} dspsr  -A -K -cont"
    dspsr_cmd += f" -U {memory}"
    dspsr_cmd += f" -b {bins}"
    dspsr_cmd += f" -c {period}"
    dspsr_cmd += f" -D {dm}"
    dspsr_cmd += f" -S {seek}"
    dspsr_cmd += f" -T {total}"
    dspsr_cmd += f" -L {total}" # Timescrunch the whole obs
    if vdif:
        psradd_coms = f"psradd -R -m time *ar -o {ar_name}"
        commands.append("j=0")
        commands.append("for i in *.hdr;")
        commands.append(f"   do {dspsr_cmd} -O ipfb_$j $i ;")
        commands.append("   j=$((j+1))")
        commands.append("done;")
        commands.append(psradd_coms)
    else:
        dspsr_cmd += f" -O {ar_name}"
        dspsr_cmd += f" {fits_dir}/*.fits"
        commands.append(dspsr_cmd)
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
    debase_cmd += "-ext debase.gg" # Aligns with what's written in the config - Don't touch unless you know what you're doing
    debase_cmd += f" -device /null"
    debase_cmd += f" {cfg['files']['converted_fits']}"
    return debase_cmd


def ppp_file_creation(cfg, depends_on=None, depend_type="afterany"):
    """Makes commands for converting ot archive file, converting back to fits and then removing baseline"""
    fits_dir = join(cfg["files"]["psr_dir"], cfg["source"]["my_pointing"])
    bins = cfg["source"]["my_bins"]
    dm = cfg["source"]["my_DM"]
    period = cfg["source"]["my_P"]
    total = cfg["source"]["total"]
    seek = cfg["source"]["seek"]
    # Fit the profile with gaussian
    bestprof_fit(cfg, cliptype="verbose")
    # Change to working directory
    commands = [f"cd {cfg['files']['psr_dir']}"]
    # Add folds to commands
    psrchive_container = comp_config['prschive_container']
    commands.append(fits_to_archive(fits_dir, cfg["files"]["archive_basename"], bins, dm, period,
                    total=total, seek=seek, vdif=cfg["run_ops"]["vdif"], container=psrchive_container))
    # Add ar -> fits conversion to commands
    commands.append(archive_to_fits(cfg["files"]["archive"], container=psrchive_container))
    # Add the baseline removal commands
    commands.append(remove_baseline(cfg))
    #Submit_job
    name = f"{cfg['obs']['id']}_{cfg['source']['name']}_archive_creation_and_debase"
    slurm_kwargs = {"time":"08:00:00"}
    modules = ["singularity"]
    mem=8192
    jid = submit_slurm(name, commands,
        slurm_kwargs=slurm_kwargs, module_list=modules, mem=mem, batch_dir=cfg["files"]["batch_dir"], depend=depends_on,
        depend_type=depend_type, vcstools_version=cfg["run_ops"]["vcstools"], submit=True)
    logger.info(f"Submitted archive/fits creation job: {name}")
    logger.info(f"job ID: {jid}")
    if depends_on:
        logger.info(f"Job depends on job id(s): {depends_on}")
    return jid, name