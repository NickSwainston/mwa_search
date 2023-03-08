#! /usr/bin/env python3
"""
Setup for mwa_search
"""
import os
import sys
from setuptools import setup
from subprocess import check_output


def read(fname):
    """Read a file"""
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

# The following two functions were taken from the repo: https://github.com/pyfidelity/setuptools-git-version/blob/master/setuptools_git_version.py


def format_version(version, fmt='{tag}.{commitcount}_{gitsha}'):
    parts = version.split('-')
    if len(parts) < 4:
        return parts[0]
    assert len(parts) in (3, 4)
    dirty = len(parts) == 4
    tag, count, sha = parts[:3]
    if count == '0' and not dirty:
        return tag
    return fmt.format(tag=tag, commitcount=count, gitsha=sha.lstrip('g'))


def get_git_version():
    git_version = check_output(
        'git describe --tags --long --dirty --always'.split()).decode('utf-8').strip()
    return format_version(version=git_version)

# Since we mostly run this on supercomputers it probably isn't correct to
# pip install all these modules
reqs = ['argparse>=1.4.0',
        'numpy>=1.13.3',
        'matplotlib>=2.1.0',
        'astropy>=2.0.2',
        'mwa-vcstools',
        'psrqpy',
        'pyyaml']

mwa_search_version = get_git_version()
#make a temporary version file to be installed then delete it
with open('version.py', 'a') as the_file:
    the_file.write('__version__ = "{}"\n'.format(mwa_search_version))

setup(
    name="mwa_search",
    version=mwa_search_version,
    description="Scripts used to search for pulsars with the Murchison Widefield Array's Voltage Capture System data",
    url="https://github.com/NickSwainston/mwa_search",
    #long_description=read('README.md'),
    python_requires='>=3.6',
    packages=['dpp', 'mwa_search'],
    package_dir={
        'dpp': 'dpp',
        'mwa_search': 'mwa_search',
    },
    package_data={'mwa_search':['data/*.npy']},
    install_requires=reqs,
    scripts=[
        'version.py',
        'scripts/calc_beamformer_benchmarks.py',
        # mwa_search
        'scripts/cold_storage_mover.py',
        'scripts/grid.py',
        'scripts/lfDDplan.py',
        'scripts/LOTAAS_wrapper.py',
        'scripts/search_launch_loop.sh',
        'scripts/rsync_rm_loop.sh',
        'scripts/bestgridpos.py',
        'scripts/find_clustered_and_known_pulsar_candidates.py',
        # dpp
        'scripts/pulsars_in_fov.py',
        'scripts/observation_processing_pipeline.py',
        'scripts/opp_status.py',
        # plotting
        'scripts/plotting/plot_obs_pulsar.py',
        'scripts/plotting/position_sn_heatmap_fwhm.py',
        # nextflow
        'nextflow/mwa_search_pipeline.nf',
        'nextflow/beamform.nf',
        'nextflow/beamform_module.nf',
        'nextflow/pulsar_search.nf',
        'nextflow/pulsar_search_module.nf',
        'nextflow/classifier.nf',
        'nextflow/classifier_module.nf',
        'nextflow/nextflow.config',
        'nextflow/data_processing_pipeline.nf',
        'nextflow/candidate_TOAs.nf',
        'nextflow/find_candidate_position.nf',
        'nextflow/pdmp.nf',
        'nextflow/dspsr_module.nf',
        'nextflow/benchmark_beamformer.nf',
        'nextflow/vcs_download.nf',
    ],
    setup_requires=['pytest-runner'],
    tests_require=['pytest']
)

os.remove('version.py')
