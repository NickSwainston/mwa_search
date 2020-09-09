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
        'psrqpy',
        'pyyaml']

mwa_search_version = get_git_version()
#make a temporary version file to be installed then delete it
with open('version.py', 'a') as the_file:
    the_file.write('__version__ = "{}"\n'.format(mwa_search_version))

setup(name="mwa_search",
      version=mwa_search_version,
      description="Scripts used to search for pulsars with the Murchison Widefield Array's Voltage Capture System data",
      url="https://github.com/NickSwainston/mwa_search",
      #long_description=read('README.md'),
      python_requires='>=3.6',
      packages=['dpp', 'mwa_search'],
      package_dir={'dpp': 'lib/dpp',
                   'mwa_search': 'lib/mwa_search'},
      install_requires=reqs,
      scripts=['version.py',
               # mwa_search
               'scripts/mwa_search/check_known_pulsars.py', 'scripts/mwa_search/cold_storage_mover.py',
               'scripts/mwa_search/grid.py', 'scripts/mwa_search/lfDDplan.py',
               'scripts/mwa_search/mwa_search_pipeline.py',
               'scripts/mwa_search/splice_wrapper.py', 'scripts/mwa_search/LOTAAS_wrapper.py',
               'scripts/mwa_search/search_launch_loop.sh', 'scripts/mwa_search/rsync_rm_loop.sh',
               'scripts/mwa_search/bestgridpos.py', 'scripts/mwa_search/data_processing_pipeline.py',
               # dpp
               'scripts/dpp/stokes_fold.py', 'scripts/dpp/make_pulsar_yaml.py',
               'scripts/dpp/pulsars_in_fov.py', 'scripts/dpp/prepfold_cmd_make.py',
               # plotting
               'scripts/plotting/plot_obs_pulsar.py',
               'scripts/plotting/position_sn_heatmap_fwhm.py',
               # nextflow
               'nextflow/beamform.nf', 'nextflow/beamform_module.nf',
               'nextflow/pulsar_search.nf', 'nextflow/pulsar_search_module.nf',
               'nextflow/classifier.nf', 'nextflow/classifier_module.nf',
               'nextflow/nextflow.config',
               'nextflow/candidate_TOAs.nf', 'nextflow/find_candidate_position.nf',
               'nextflow/pdmp.nf', 'nextflow/dspsr_module.nf',
               'nextflow/mwa_search_pipeline.nf', 'nextflow/data_processing_pipeline.nf'],
      setup_requires=['pytest-runner'],
      tests_require=['pytest']  # , 'nose']
      )

os.remove('version.py')
