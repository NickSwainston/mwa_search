#!/usr/bin/env python3
import yaml_helper

def yaml_check_args(kwargs)
    """Makes assertions and changes to the kwargs from data_processing_pipeline.py script"""

    if kwargs["yaml"]:
        _kwargs = yaml_helper.from_yaml(kwargs["yaml"])
        for key in _kwargs.keys():
            if key not in kwargs.keys():
                kwargs[key] = _kwargs[key]

    unavailable=[]
    if not kwargs["run_dir"]:
        unavailable.append("run_dir")
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
        pipe_helper.initiate_yaml(kwargs)

    return kwargs

def parse_unavailable(unavailable):
    """raises an error if there is anything in the input list"""
    if parse_unavailable:
        raise ValueError(f"The following parameters need to be provided: {unavailable}")
        