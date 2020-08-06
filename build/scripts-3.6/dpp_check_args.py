#!python
import logging
import glob


import yaml_helper
import pipe_helper


logger = logging.getLogger(__name__)


class FileNotExistsError:
    """Raise when a file does not exist"""


def yaml_check_args(kwargs):
    """Makes assertions and changes to the kwargs from data_processing_pipeline.py script"""

    if kwargs["yaml"]:
        kwargs = yaml_helper.from_yaml(kwargs["yaml"])

    else:
        unavailable = []
        if not kwargs["run_dirs"]:
            unavailable.append("run_dirs")
        if not kwargs["obsid"]:
            unavailable.append("obsid")
        if not kwargs["cal_id"]:
            unavailable.append("cal_id")
        if not kwargs["obs_beg"]:
            unavailable.append("obs_beg")
        if not kwargs["obs_end"]:
            unavailable.append("obs_end")

        if not kwargs["cand"]:
            if not kwargs["pulsar"]:
                unavailable.append("pulsar")

        parse_unavailable(unavailable)

        if "initiated_yaml" not in kwargs.keys():
            kwargs = pipe_helper.initiate_pipe(kwargs)

    return kwargs


def parse_unavailable(unavailable):
    """raises an error if there is anything in the input list"""
    if unavailable:
        raise ValueError(
            f"The following parameters need to be provided: {unavailable}")
