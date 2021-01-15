import logging
from os.path import join

from vcstools.job_submit import submit_slurm
from vcstools.config import load_config_file

comp_config = load_config_file()
logger = logging.getLogger(__name__)


def generate_archive_name(cfg):
    return f"{cfg['obs']['id']}_{cfg['source']['my_pointing']}_{cfg['source']['name']}_b{cfg['source']['my_bins']}"


def fits_to_archive(fits_dir, ar_name, bins, dm, period, memory=4000, total=1.0, seek=0, vdif=False):
    """Returns bash commands to fold on a fits file usinf dspsr"""
    commands = []
    commands.append(f"singularity shell {comp_config['prschive_container']}") # Open the container
    dspsr_cmd = "dspsr  -A -K -cont"
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
        commands.append(f"   do {dspsr_coms} -O ipfb_$j $i ;")
        commands.append("   j=$((j+1))")
        commands.append("done;")
        commands.append(psradd_coms)
    else:
        dspsr_com += f" -O {ar_name}"
        dspsr_com += f" {fits_dir}/*.fits"
        commands.append(dspsr_cmd)
    commands.append("exit") # Exit the container
    return commands


def archive_to_fits(ar_file, extension="fits"):
    """Returns bash commands to turn an arhive file to a fits file"""
    comands = []
    commands.append(f"singularity shell {comp_config['prschive_container']}") # Open the container
    pam_cmd = "pam -a PSRFITS"
    pam_cmd += f" -e {extension}"
    pam_cmd += f" {ar_file}"
    commands.append(pam_cmd)
    commands.append("exit") # Exit the container
    return commands


def remove_baseline(cfg):
    """Removes baseline RFI from a fits file"""
    ar_name = generate_archive_name(cfg)
    on_pulse = cfg["source"]["gfit"]["comp_idx"]["0"]
    debase_cmd = "pmod -debase"
    debase_cmd += f" -onpulsef '{on_pulse[0]} {on_pulse[-1]}'"
    debase_cmd += "-ext debase.gg"
    debase_cmd += f" -device /null"
    debase_cmd += f" {ar_name}.fits"
    commands.append(debase_cmd)
    return commands


def ppp_file_creation(cfg, depends_on=None, depend_type="afterany"):
    """Makes commands for converting ot archive file, converting back to fits and then removing baseline"""
    fits_dir = join(cfg["files"]["psr_dir"], cfg["source"]["my_pointing"])
    ar_name = generate_archive_name(cfg)
    ar_file = f"{ar_name}*.ar"
    bins = cfg["source"]["my_bins"]
    dm = cfg["source"]["my_DM"]
    period = cfg["source"]["my_P"]
    total = cfg["source"]["total"]
    seek = cfg["source"]["seek"]
    # Change to working directory
    commands = [f"cd {cfg['files']['psr_dir]}"]
    # Add folds to commands
    commands.append(fits_to_archive(fits_dir, ar_name, bins, dm, period, total=total, seek=seek, vdif=cfg["run_ops"]["vdif"]))
    # Add ar -> fits conversion to commands
    commands.append(archive_to_fits(ar_file))
    # Add the baseline removal commands
    commands.append(remove_baseline(cfg))
    #Submit_job
    name = f"{cfg['obs']['id']}_{cfg['source']['name']}_archive_creation_and_debase"
    slurm_kwargs = {"time":"08:00:00"}
    modules = ["singularity"]
    mem=8192
    jid = submit_slurm(name, cmds,
        slurm_kwargs=slurm_kwargs, module_list=modules, mem=mem, batch_dir=cfg["files"]["batch_dir"], depend=depends_on,
        depend_type=depend_type, vcstools_version=cfg["run_ops"]["vcstools"], submit=True)
    logger.info(f"Submitted archive/fits creation job: {name}")
    logger.info(f"job ID: {jid}")
    if depends_on:
        logger.info(f"Job depends on job id(s): {depends_on}")
    return jid, name