import logging
from os.path import join, basename

from vcstools.job_submit import submit_slurm
from rm_synthesis import rm_synth_pipe

logger = logging.getLogger(__name__)


def RM_synth(cfg):
    # Get some basic things in order
    my_comp = cfg["source"]["my_component"]
    comp_min = cfg["source"]["gfit"]["comp_idx"][my_comp][0]
    comp_max = cfg["source"]["gfit"]["comp_idx"][my_comp][-1]
    prof_len = len(cfg["source"]["gfit"]["profile"])
    comp_range = comp_max - comp_min
    if comp_range < 7: # Errors arise if phase range is too small
        comp_min -= (7-comp_range)/2
        comp_max += (7-comp_range)/2
    phase_min = comp_min/prof_len
    phase_max = comp_max/prof_len
    # RM synthesis - Initial
    rms_kwargs = {}
    # Fill out the useless stuff
    rms_kwargs["write"] = False
    rms_kwargs["plot_gfit"] = False
    rms_kwargs["plot_rm"] = False
    rms_kwargs["force_single"] = False
    rms_kwargs["keep_QUV"] = False
    # Find the important things
    rms_kwargs["phase_ranges"] = (phase_min, phase_max)
    rms_kwargs["phi_steps"] = 10000
    rms_kwargs["phi_range"] = (-300, 300)
    rms_kwargs["label"] = "RMsynth_initial"
    rms_kwargs["work_dir"] = cfg["files"]["psr_dir"]
    rms_kwargs["archive"] = cfg["files"]["archive"]
    rm_dict, _ = rm_synth_pipe(rms_kwargs)
    rm_initial = rm_dict["0"]["rm"]
    # RM synthesis - Fine
    rms_kwargs["phi_range"] = (rm_initial-10, rm_initial+10)
    rms_kwargs["label"] = "RMsynth_final"
    rms_kwargs["plot_rm"] = True
    rm_dict, _ = rm_synth_pipe(rms_kwargs)
    rm_final = rm_dict["0"]["rm"]
    rm_final_e = rm_dict["0"]["rm_e"]
    cfg["pol"]["RM"] = rm_final
    cfg["pol"]["RM_e"] = rm_final_e
    logger.info(f"Calculated RM through synthesis: {rm_final} +/- {rm_final_e}")


def RM_cor(cfg, depends_on=None, depend_type="afterany"):
    """Submits an RM correction job"""
    # Create PSRSALSA commands
    my_comp = cfg["source"]["my_component"]
    myfits = join(cfg['files']['psr_dir'], cfg["files"]["archive"])
    on_pulse_min = cfg["source"]["gfit"]["comp_idx"][my_comp][0]
    on_pulse_max = cfg["source"]["gfit"]["comp_idx"][my_comp][-1]
    profile_ps = basename(cfg['files']['ppol_profile_ps'])
    polar_prof_ps = basename(cfg['files']['ppol_polar_profile_ps'])
    rm = cfg["pol"]["RM"]
    ppol_cmd = "ppol -v"
    ppol_cmd += f" -onpulse '{on_pulse_min} {on_pulse_max}'"
    ppol_cmd += f" -header 'rm {rm}'"
    ppol_cmd += " -TSCR -FSCR"
    ppol_cmd += f" -ext paswing"
    ppol_cmd += f" -device {profile_ps}/cps" # is the profile
    ppol_cmd += f" -device2 {polar_prof_ps}/cps" # is the polarimetry profile
    ppol_cmd += f" {cfg['files']['debased_fits']}"
    commands = [f"cd {cfg['files']['psr_dir']}"]
    commands.append(ppol_cmd)
    # Submit the job
    name = f"{cfg['obs']['id']}_{cfg['source']['name']}_RM_correction"
    slurm_kwargs = {"time":"01:00:00"}
    modules = ["psrsalsa"]
    mem=8192
    jid = submit_slurm(name, commands,
        slurm_kwargs=slurm_kwargs, module_list=modules, mem=mem, batch_dir=cfg["files"]["batch_dir"], depend=depends_on,
        depend_type=depend_type, vcstools_version=cfg["run_ops"]["vcstools"], submit=True)
    logger.info(f"Submitted Rm correction creation job: {name}")
    logger.info(f"job ID: {jid}")
    if depends_on:
        logger.info(f"Job depends on job id(s): {depends_on}")
    cfg["completed"]["RM"] = True
    return jid, name