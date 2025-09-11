"""
File writing utilities for primers/barcodes
"""
import csv
from typing import List, Tuple, Dict
import string
import os
from datetime import datetime

def write_idt_csv(primer_pairs: List[Tuple[str, str]], outfile: str, scale: str, purification: str) -> None:
    with open(outfile, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(["Name", "Sequence", "Scale", "Purification"])
        for i, pair in enumerate(primer_pairs, 1):
            bid = f"BC{i:02d}"
            fwd_seq, rev_seq = pair
            writer.writerow([f"BC{bid}_F)", fwd_seq, scale, purification])
            writer.writerow([f"BC{bid}_R)", rev_seq, scale, purification])

def write_failures_csv(failures: List[Dict[str, str]], outfile: str) -> None:
    with open(outfile, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(["Barcode_or_Pair", "Reason"])
        for row in failures:
            writer.writerow([row.get("barcode", ""), row.get("reason", "")])

def write_fasta(primer_pairs: List[Tuple[str, str]], outfile: str) -> None:
    with open(outfile, 'w') as fh:
        for i, pair in enumerate(primer_pairs, 1):
            bid = f"BC{i:02d}"
            fh.write(f">{bid}_F|barcode={pair[0]}\n{pair[1]}\n")
            fh.write(f">{bid}_R|barcode={pair[0]}\n{pair[1]}\n")

def well_ids(rows: int, cols: int) -> List[str]:
    letters = string.ascii_uppercase[:rows]
    return [f"{letters[r]}{c}" for r in range(rows) for c in range(1, cols+1)]

def write_plate_csvs(primer_pairs: List[Tuple[str, str]], outfile: str, rows: int, cols: int, scale: str, purification: str, split: bool) -> None:
    primers = []
    for i, pair in enumerate(primer_pairs, 1):
        bid = f"BC{i:02d}"
        primers.append((f"{bid}_F", pair[0]))
        primers.append((f"{bid}_R", pair[1]))
    cap = rows * cols
    if not split and len(primers) > cap:
        raise SystemExit(f"Requested plate has capacity {cap} wells, but {len(primers)} primers were generated. Use --plate-split or reduce --num.")
    def write_one(fname: str, chunk):
        wells = well_ids(rows, cols)
        with open(fname, 'w', newline='') as fh:
            w = csv.writer(fh)
            w.writerow(["Well", "Name", "Sequence", "Scale", "Purification"])
            for i, (nm, seq) in enumerate(chunk):
                if i >= len(wells): break
                w.writerow([wells[i], nm, seq, scale, purification])
    if len(primers) <= cap:
        write_one(outfile, primers)
    else:
        nplates = (len(primers) + cap - 1) // cap
        base, ext = (outfile.rsplit('.', 1) + ["csv"])[:2]
        for p in range(nplates):
            chunk = primers[p*cap:(p+1)*cap]
            fname = f"{base}_{p+1}.{ext}"
            write_one(fname, chunk)

def write_outputs(barcode_pairs, failures, config):
    flanking_seq_fwd = config.get("flanking_seq_fwd", "")
    flanking_seq_rev = config.get("flanking_seq_rev", "")
    template_binding_seq_fwd = config.get("template_binding_seq_fwd", "")
    template_binding_seq_rev = config.get("template_binding_seq_rev", "")
    filename_idt = config.get("filename_idt", "")
    filename_failed = config.get("filename_failed", "")
    filename_fasta = config.get("filename_fasta", "")
    filename_plate = config.get("filename_plate", "")
    plate_rows = config.get("plate_rows", 8)
    plate_cols = config.get("plate_cols", 12)
    plate_split = config.get("plate_split", False)
    scale = config.get("scale", "25nm")
    purification = config.get("purification", "STD")
    timestamp_outdir = config.get("timestamp_outdir", True)
    outdir = config.get("outdir", "../output")

    if timestamp_outdir:
        ts = datetime.now().strftime("%y%m%d_%H%M")
        outdir = os.path.join(outdir, ts)
    os.makedirs(outdir, exist_ok=True)

    primer_pairs = []
    for i in range(len(barcode_pairs)):
        fwd_seq = f"{flanking_seq_fwd}{barcode_pairs[i][0]}{template_binding_seq_fwd}"
        rev_seq = f"{flanking_seq_rev}{barcode_pairs[i][1]}{template_binding_seq_rev}"
        primer_pairs.append((fwd_seq, rev_seq))

    if filename_idt:
        filepath_idt = os.path.join(outdir, filename_idt)
        write_idt_csv(primer_pairs, filepath_idt, scale=scale, purification=purification)
    if filename_failed:
        filepath_failed = os.path.join(outdir, filename_failed)
        write_failures_csv(failures, filepath_failed)
    if filename_fasta:
        filepath_fasta = os.path.join(outdir, filename_fasta)
        write_fasta(primer_pairs, filepath_fasta)
    if filename_plate:
        filepath_plate = os.path.join(outdir, filename_plate)
        write_plate_csvs(primer_pairs, filepath_plate, plate_rows, plate_cols, scale=scale, purification=purification, split=plate_split)