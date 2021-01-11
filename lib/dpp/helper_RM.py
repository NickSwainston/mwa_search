#!/usr/bin/env python3
import logging
from os.path import join

from vcstools.prof_utils import subprocess_pdv, get_from_ascii, auto_gfit
from dpp.helper_archive import generate_archive_name
from rm_synthesis import rm_synth_main

logger = logging.getLogger(__name__)


def archive_prof_fit(archive_file, ascii_outfile, kwargs=None):
    """Fits a profile to an archive. Takes prof_utils.auto_gfit kwargs"""
    # Convert archive to .txt and read it
    subprocess_pdv(archive_file, outfile=ascii_outfile)
    profile, _ = get_from_ascii(ascii_outfile)
    # Gaussian fit
    fit_dict = auto_gfit(profile, **kwargs)
    return fit_dict

def RM_synth(cfg):
    archive_name = generate_archive_name(cfg)
    archive_file = join(cfg["run_ops"]["psr_dir"], f"{archive_name}.ar")
    ascii_outfile = join(cfg["run_ops"]["psr_dir"], f"{archive_name}.txt")
    gfit_kwargs = {"cliptype":"verbose", "period":cfg["source"]["my_P"], "plot_name":f"{archive_name}.png"}
    # Do the fitting
    profile_fit = archive_prof_fit(archive_file, ascii_outfile, kwargs=gfit_kwargs)
    cfg["source"]["gfit"] = profile_fit
    # RM synthesis - Initial
    rms_kwargs = {}
    rms_kwargs["phase_ranges"] = f"{profile_fit['component_1'][0]} {profile_fit['component_1'][1]}"
    rms_kwargs["phi_steps"] = 10000
    rms_kwargs["phi_range"] = (-300, 300)
    rms_kwargs["label"] = f"{archive_name}_RMsynth_initial"
    rms_kwargs["write"] = False
    rms_kwargs["plot"] = False
    rms_kwargs["work_dir"] = cfg["run_ops"]["psr_dir"]
    rm_dict = rm_synth_pipe(rms_kwargs)
    # Find the important things
    rm_inital = rm_dict["0"]["rm"]
    rms_kwargs["phi_range"] = (rm_initial-10, rm_final+10)
    rms_kwargs["label"] = f"{archive_name}_RMsynth_final"
    rms_kwargs["write"] = True
    rms_kwargs["plot"] = True
    # RM synthesis - Fine
    rm_dict = rm_synth_pipe(rms_kwargs)
    rm_final = rm_dict["0"]["rm"]
    rm_final_e = rm_dict["0"]["rm_e"]
    cfg["pol"]["RM"] = rm_final
    cfg["pol"]["RM_e"] = rm_final_e
    logger.info(f"Calculated RM through synthesis: {rm_final} +/- {rm_final_e}")


def RM_cor(cfg, dpeends_on=None, depend_type="afterany"):
    """Submits an RM correction job"""
    # Create PSRSALSA commands
    myfits = join(cfg['run_ops']['psr_dir'], f"{generate_archive_name(cfg)}.fits")
    on_pulse_min = cfg["source"]["gfit"]["component_1"][0]
    on_pulse_max = cfg["source"]["gfit"]["component_1"][1]
    outfile = join(cfg['run_ops']['psr_dir'], f"{generate_archive_name(cfg)}_rmcorrected.txt")
    plotname = join(cfg['run_ops']['psr_dir'], f"{generate_archive_name(cfg)}_rmcorrected.ps")
    rm = cfg["pol"]["RM"]
    ppol_cmd = "ppol -v"
    ppol_cmd += f" -onpulse '{on_pulse_min} {on_pulse_max}'"
    ppol_cmd += f" -header 'rm {rm}'"
    ppol_cmd += " -TSCR -FSCR"
    ppol_cmd += f" -ofile {outifle}"
    ppol_cmd += " -device /null -device2 /ps"
    cmds = []
    cmds.append(f"cd {cfg['run_ops']['psr_dir']}")
    cmds.append(ppol_cmd)
    # Submit the job
    name = f"{cfg['obs']['id']}_{cfg['source']['name']}_RM_correction"
    slurm_kwargs = {"time":"01:00:00"}
    modules = ["psrsalsa"]
    mem=8192
    jid = submit_slurm(name, cmds,
        slurm_kwargs=slurm_kwargs, module_list=modules, mem=mem, batch_dir=cfg["run_ops"]["batch_dir"], depend=depends_on,
        depend_type=depend_type, vcstools_version=cfg["run_ops"]["vcstools"], submit=True)
    logger.info(f"Submitted Rm correctoin creation job: {name}")
    logger.info(f"job ID: {jid}")
    if depends_on:
        logger.info(f"Job depends on job id(s): {depends_on}")
    cfg["RM"] = True
    return jid, name