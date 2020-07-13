#!/usr/bin/env python3
import yaml
import os

def dump_to_yaml(mydict, filepath):
    with open(filepath, 'w') as f:
        yaml.dump(mydict, f, default_flow_style=False)

def from_yaml(filepath):
    with open(filepath) as f:
        my_dict = yaml.safe_load(f)
    return my_dict