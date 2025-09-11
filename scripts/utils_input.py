"""
functions for loading and validating the config
"""

import argparse
import yaml
from typing import Dict, Any
import logging

def parse_args():
    parser = argparse.ArgumentParser(description="Generate ONT-tailed, barcoded gyrB primers from a config file.")
    parser.add_argument("config_path", nargs="?", default="config.yaml", help="Path to config YAML file (default: config.yaml)")
    return parser.parse_args()

def validate_config(config: Dict[str, Any]):
    #TODO add more validation and required fields
    required_fields = {
        "n_primers": int,
    }
    errors = []
    for field, field_type in required_fields.items():
        if field not in config:
            errors.append(f"Missing required field: {field}")
        elif not isinstance(config[field], field_type):
            errors.append(f"Field '{field}' must be of type {field_type.__name__}")
    if errors:
        for error in errors:
            logging.error(f"Config error: {error}")
        raise SystemExit("Please fix config.yaml and try again")

def load_config(path: str):
    with open(path, "r") as f:
        config = yaml.safe_load(f)
    validate_config(config)
    return config