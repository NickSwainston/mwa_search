import logging

logger = logging.getLogger(__name__)


def RVM_fit(cfg):
    """Calculates parameters for a PSRSALSA ppolFit job and submits it"""
    alpha = cfg["pol"]["alpha"]
    beta = cfg["pol"]["beta"]
    if not alpha or not beta: # Initial
        trials = 200
        alpha_range = np.array((0, 180))
        beta_range = np.array((-30, 30))
        #Decide the longitude range to fit
        my_comp = cfg["source"]["my_component"]
        component_min = cfg["source"["gfit"]["comp_idx"][my_comp][0] * 360/len(cfg["source"]["gfit"]["profile"])
        component_max = cfg["source"["gfit"]["comp_idx"][my_comp][-1] * 360/len(cfg["source"]["gfit"]["profile"])
        l_cmd = f" -l {component_min} 1"
        maxdl_cmd = f" -maxdl {component_max - component_min}"
        chigrid_file = cfg['file']['chigrid_initial_ps']
        paswing_file = cfg['file']['paswing_initial_ps']
        outfile = cfg['file']['RVM_fit_initial']
    else: # Final
        trials = 400
        alpha_range = np.array((alpha - 20, alpha + 20))
        beta_range = np.array((beta - 10, beta + 10))
        alpha_range.clip(0, 180)
        beta_range.clip(-30, 30) # forcing the range to reasonable values
        l_cmd = f" -l {cfg['pol']['l0'] - 10}"
        maxdl_cmd = " -maxdl 20"
        chigrid_file = cfg['file']['chigrid_final_ps']
        paswing_file = cfg['file']['paswing_final_ps']
        outfile = cfg['file']['RVM_fit_final']
    # Create the job commands
    cmds = [f"cd {cfg['files']['psr_dir']}"]
    ppol_cmd = "ppolFit -showwedge"
    ppol_cmd += f" -g '{trials} {trials}'"
    ppol_cmd += f" -A '{alpha_range[0]} {alpha_range[1]}'" # Alpha range
    ppol_cmd += f" -B '{beta_range[0]} {beta_range[1]}'" # Beta range
    ppol_cmd += l_cmd # longitude start and step size
    ppol_cmd += maxdl_cmd #longitude search range
    ppol_cmd += " -best" # return the best fit values
    ppol_cmd += f" -device1 {chigrid_file}/cps"
    ppol_cmd += f" -device2 {paswing_file}.ps/cps"
    ppol_cmd += f" -device1res '900 900'"
    ppol_cmd += f" *.paswing"
    ppol_cmd += f" > {outfile}" #write to stdout


def read_RVM_fit(cfg):
    pass


ppolFit -g "[NUM_ALPHA_PTS] [NUM_BETA_PTS]" -l "[START_VALUE] [STEP_SIZE]" -showwedge -best [ASCII_OUT_FILE_FROM_PPOL]
#This can be done iteratively and with better results by giving estimates and ranges for a, b, l0, pa0 and pulse width (wmp)


