"""
functions for checking barcodes
"""
import random

def check_gc(barcode, gc_min, gc_max):
    gc_count = sum(1 for b in barcode if b in 'GC')
    gc_content = gc_count / len(barcode)
    if not (gc_min <= gc_content <= gc_max):
        return False, f'GC content {gc_content:.2f} out of range ({gc_min}-{gc_max})'
    return True, None

def check_edit_distance(barcode, accepted_pairs, edit_min):
    for fwd, rev in accepted_pairs:
        dist_fwd = sum(1 for a, b in zip(barcode, fwd) if a != b)
        dist_rev = sum(1 for a, b in zip(barcode, rev) if a != b)
        if dist_fwd < edit_min or dist_rev < edit_min:
            return False, f'edit distance to existing barcode less than {edit_min}'
    return True, None

def get_fwd_primer_seq(barcode, flanking_seq, template_binding_seq):
    return f"{flanking_seq}{barcode}{template_binding_seq}"