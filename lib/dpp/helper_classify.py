import logging
from os.path import join, exists

from vcstools.job_submit import submit_slurm
from dpp.helper_files import setup_classify
from dpp.helper_relaunch import relaunch_ppp
from dpp.helper_config import dump_to_yaml


logger = logging.getLogger(__name__)


def add_classify_to_commands(cfg, container="/pawsey/mwa/singularity/lofar_pulsar_ml/lofar_pulsar_ml.sif"):
    """Makes the classify commands"""
    container_launch = f"singularity exec -e {container}"
    cmds = [f"cd {cfg['run_ops']['classify_dir']}"]
    # The container needs all this stuff to run properly for some reason
    #singularity_launch = 'set +u; env - PATH="$PATH"'
    #singularity_launch += ' SINGULARITYENV_TMP="$TMP"'
    #singularity_launch += ' SINGULARITYENV_TMPDIR="$TMPDIR"'
    #singularity_launch += ' SINGULARITYENV_NXF_DEBUG=${NXF_DEBUG:=0}'
    #singularity_launch += ' ${SINGULARITY_BINDPATH:+SINGULARITY_BINDPATH="$SINGULARITY_BINDPATH"}'
    #singularity_launch += ' ${SINGULARITYENV_LD_LIBRARY_PATH:+SINGULARITYENV_LD_LIBRARY_PATH="$SINGULARITYENV_LD_LIBRARY_PATH"}'
    #singularity_launch += ' /pawsey/mwa/singularity/lofar_pulsar_ml/lofar_pulsar_ml.sif'
    #cmds.append(singularity_launch)
    # Run the feature extractor
    cmds.append("REALPATH=`realpath feature_extraction.arff`")
    cmds.append(f"{container_launch} python /usr/local/bin/PulsarFeatureLab.py -d `pwd` -f feature_extraction.arff -t 6 -c 3 --meta --arff")
    #Run the features through the 5 models
    cmds.append("for i in {1..5}; do")
    cmds.append(f"   {container_launch} java -jar /usr/local/bin/LOTAASClassifier.jar -m /home/soft/models/V1.3.1model${{i}}.model -p ${{REALPATH}} -a 1 -d") # ${LOTAAS_MLC_MODEL_DIR}
    cmds.append("   if [ -f '\${REALPATH%arff}positive' ]; then")
    cmds.append("       mv \${REALPATH%arff}positive feature_extraction_m${i}.positive")
    cmds.append("   fi")
    cmds.append("   if [ -f '\${REALPATH%arff}negative' ]; then")
    cmds.append("       mv \${REALPATH%arff}negative feature_extraction_m${i}.negative")
    cmds.append("   fi")
    cmds.append("done")
    return cmds


def submit_classify(cfg):
    """launched a classify job"""
    # Make the commands for the job
    cmds = add_classify_to_commands(cfg)
    # Work out some things for the job
    name = f"{cfg['files']['file_precursor']}_classify"
    slurm_kwargs = {"time":"00:30:00"}
    modules = ["singularity"]
    mem = 8192
    # Submit Job
    jid = submit_slurm(name, cmds,
        slurm_kwargs=slurm_kwargs, module_list=modules, mem=mem, batch_dir=cfg["files"]["batch_dir"], load_vcstools=False, submit=True)
    logger.info(f"Submitted classiy job: {name}")
    logger.info(f"Job ID: {jid}")
    return jid, name


def read_classifications(cfg):
    """Reads the output of the classifier and updates cfg with the information"""
    negfile = join(cfg["files"]["classify_dir"], "feature_extraction.negative")
    posfile = join(cfg["files"]["classify_dir"], "feature_extraction.positive")
    try:
        with open(posfile, "r") as f:
            pos = f.readlines()
    except FileNotFoundError as e:
        if not exists(negfile): # A least one of the pos and neg files should exist
                raise FileNotFoundError(f"Classifier outputs not found in dir: {cfg['run_ops']['classify_dir']}")
        else:
            pos = []
    for pointing in cfg["folds"].keys():
        # Count positive model classifications
        cfg["folds"][pointing]["classifier"] = sum(pointing in s for s in pos)
        logger.debug(f"{pointing} Positive models found: {cfg['folds'][pointing]['classifier']}")


def classify_main(cfg):
    """initiates and launches a classify job"""
    setup_classify(cfg)
    jid, _ = submit_classify(cfg)
    cfg["completed"]["classify"] = True
    return jid