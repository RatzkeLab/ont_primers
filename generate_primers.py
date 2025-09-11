"""
CLI orchestrator for ONT-tailed, barcoded gyrB primer generation
"""
import logging
from scripts.utils_input import parse_args, load_config
from scripts.utils_generate import generate_barcodes
from scripts.utils_output import write_outputs  

def run():
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = parse_args()
    logging.info(f"Loading config from {args.config_path}")
    config = load_config(args.config_path)
    logging.info("Generating barcodes and primers...")
    accepted_pairs, failures = generate_barcodes(config)
    logging.info(f"Generated {len(accepted_pairs)} valid pairs, {len(failures)} failures.")
    logging.info("Writing outputs...")
    write_outputs(accepted_pairs, failures, config)
    logging.info("Done.")

if __name__ == "__main__":
    run()