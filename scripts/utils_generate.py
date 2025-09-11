"""
Barcode and primer generation logic
"""

import logging
import csv
import random
from .utils_sequence import check_gc, check_edit_distance

def generate_random_barcode(length, seed):
    random.seed(seed)
    return ''.join(random.choices('ACGT', k=length))

def get_provided_barcodes(filepath):
    provided_barcodes = []
    if filepath and isinstance(filepath, str):
        with open(filepath, "r") as f:
            reader = csv.reader(f)
            #skip header if present
            header = next(reader)
            provided_barcodes = [row[1] for row in reader if row]
        logging.info(f"Loaded {len(provided_barcodes)} barcodes from {filepath}")
    else:
        logging.info("No barcode file provided, generating barcodes randomly")    
    return provided_barcodes

def generate_barcodes(config):
    n_primers = config.get("n_primers", 32)
    barcode_length = config.get("barcode_length", 10)
    seed = config.get("random_seed", 42)
    n_attempts = config.get("n_attempts", 1000)
    symmetric = config.get("symmetric_pairs", True)
    use_barcodes_from = config.get("use_barcodes_from", None)

    flanking_seq_fwd = config.get("flanking_seq_fwd", "")
    flanking_seq_rev = config.get("flanking_seq_rev", "")
    template_binding_seq_fwd = config.get("template_binding_seq_fwd", "")
    template_binding_seq_rev = config.get("template_binding_seq_rev", "")

    gc_min = config.get("gc_min", 0)
    gc_max = config.get("gc_max", 1)
    avoid_seqs = [m.upper() for m in config.get("avoid_seqs", [])]
    avoid_at_junction_seqs = [m.upper() for m in config.get("avoid_at_junction_seqs", [])]
    homopolymer_max = config.get("homopolymer_max", float('inf'))
    junction_homopolymer_max = config.get("junction_homopolymer_max", float('inf'))
    
    edit_min = config.get("min_edit_distance", 1)
    max_hetero_stretch = config.get("max_hetero_stretch", float('inf'))

    provided_barcodes = get_provided_barcodes(use_barcodes_from)
    accepted_pairs = []
    failures = []
    
    logging.info(f"Generating {n_primers} primer pairs with barcode_length={barcode_length}")
    for i_attempt in range(n_attempts):
        if len(accepted_pairs) < n_primers:
            if i_attempt < len(provided_barcodes):
                barcode = provided_barcodes[i_attempt]
            else:
                barcode = generate_random_barcode(barcode_length, seed)

            gc_check = check_gc(barcode, gc_min, gc_max)
            if gc_check[0] is False:
                failures.append({"barcode": barcode, "reason": gc_check[1]})
                continue
            edit_check = check_edit_distance(barcode, accepted_pairs, edit_min)
            if edit_check[0] is False:
                failures.append({"barcode": barcode, "reason": edit_check[1]})
                continue
            rev_barcode = str(barcode)[:-3] #TODO add option for asymmetric barcode pairs #FIXME just trimming to satisfy idt requirements
            accepted_pairs.append((barcode, rev_barcode))
    return accepted_pairs, failures
