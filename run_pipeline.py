#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_pipeline.py — eDNA Multi-Barcode Ensemble Index Integrated Pipeline
========================================================================
Executes the full preprocessing and ensemble index workflow for multiple
barcodes in a single command.

Step 1  (per barcode):
    1a. merge_otu_taxonomy.py    — attach taxonomy to OTU table
    1b. merge_otu_taxon.py       — collapse OTUs to taxon level
    1c. otutab_normalize.py      — convert counts to relative abundance

Step 2  (cross-barcode):
    2.  eDNA_ensemble_index_pipeline.py — compute eDNA ensemble index

Reference methodology:
    Kelly et al. (2019), Djurhuus et al. (2020)

Usage example
-------------
python run_pipeline.py \\
    --barcodes 16S 18S COI 12S \\
    --otutab-dir data/ \\
    --taxonomy-dir data/ \\
    --workdir intermediate/ \\
    --outdir results/ensemble/ \\
    --matrix-map "B:Biofilm,M:Membrane,S:Sediment" \\
    --exclude-matrix COI:Biofilm 12S:Biofilm \\
    --prefix Yongding_2024
"""

import argparse
import subprocess
import sys
import os
from pathlib import Path


# ─────────────────────────────────────────────────────────────────────────────
# Logging helpers
# ─────────────────────────────────────────────────────────────────────────────

def log(msg: str, level: str = "INFO"):
    tag = {
        "INFO":  "[ INFO ]",
        "STEP":  "[ STEP ]",
        "OK":    "[  OK  ]",
        "WARN":  "[ WARN ]",
        "ERR":   "[ ERR  ]",
    }.get(level, "[      ]")
    print(f"{tag} {msg}", flush=True)


def separator(char="─", width=60):
    log(char * width, "INFO")


# ─────────────────────────────────────────────────────────────────────────────
# Subprocess runner
# ─────────────────────────────────────────────────────────────────────────────

def run_script(script_name: str, args: list, step_label: str):
    """
    Run a pipeline script as a subprocess.
    Exits with code 1 if the script returns a non-zero exit code.
    """
    # Resolve script path relative to this file's location
    script_dir = Path(__file__).parent
    script_path = script_dir / script_name

    if not script_path.exists():
        log(f"Script not found: {script_path}", "ERR")
        log("Ensure all pipeline scripts are in the same directory as run_pipeline.py.", "ERR")
        sys.exit(1)

    cmd = [sys.executable, str(script_path)] + [str(a) for a in args]
    log(f"Running: {' '.join(cmd)}", "INFO")

    result = subprocess.run(cmd, text=True)

    if result.returncode != 0:
        log(f"Script failed at step: {step_label}", "ERR")
        log(f"Exit code: {result.returncode}", "ERR")
        sys.exit(1)

    log(f"Completed: {step_label}", "OK")


# ─────────────────────────────────────────────────────────────────────────────
# Per-barcode preprocessing (Steps 1a–1c)
# ─────────────────────────────────────────────────────────────────────────────

def preprocess_barcode(barcode: str,
                       otutab_dir: Path,
                       taxonomy_dir: Path,
                       workdir: Path) -> Path:
    """
    Run Steps 1a–1c for a single barcode.
    Returns the path to the normalized relative abundance table.

    Expected input files:
        {otutab_dir}/{barcode}_otutab.txt
        {taxonomy_dir}/{barcode}_taxonomy.txt
    """
    separator()
    log(f"Preprocessing barcode: {barcode}", "STEP")
    separator()

    otutab_file   = otutab_dir   / f"{barcode}_otutab.txt"
    taxonomy_file = taxonomy_dir / f"{barcode}_taxonomy.txt"
    merged_file   = workdir / f"{barcode}_merged.txt"
    taxon_file    = workdir / f"{barcode}_taxon.txt"
    norm_file     = workdir / f"{barcode}_normalized.txt"

    # ── Validate inputs ────────────────────────────────────────────────────
    for fpath, label in [(otutab_file, "OTU table"), (taxonomy_file, "taxonomy file")]:
        if not fpath.exists():
            log(f"Missing {label} for barcode {barcode}: {fpath}", "ERR")
            log("Check that input files follow the naming convention "
                "{BARCODE}_otutab.txt / {BARCODE}_taxonomy.txt", "ERR")
            sys.exit(1)

    # ── Step 1a: merge taxonomy ────────────────────────────────────────────
    log(f"[{barcode}] Step 1a: Merging taxonomy with OTU table ...", "STEP")
    run_script(
        "merge_otu_taxonomy.py",
        ["-t", taxonomy_file, "-o", otutab_file, "-out", merged_file],
        step_label=f"{barcode} — merge taxonomy"
    )

    # ── Step 1b: collapse OTUs to taxon level ─────────────────────────────
    log(f"[{barcode}] Step 1b: Collapsing OTUs to taxon level ...", "STEP")
    run_script(
        "merge_otu_taxon.py",
        ["-i", merged_file, "-o", taxon_file],
        step_label=f"{barcode} — collapse to taxon"
    )

    # ── Step 1c: normalize to relative abundance ───────────────────────────
    log(f"[{barcode}] Step 1c: Normalizing to relative abundance ...", "STEP")
    run_script(
        "otutab_normalize.py",
        ["-i", taxon_file, "-o", norm_file, "-m", "relative"],
        step_label=f"{barcode} — normalize"
    )

    log(f"[{barcode}] Preprocessing complete → {norm_file}", "OK")
    return norm_file


# ─────────────────────────────────────────────────────────────────────────────
# Cross-barcode ensemble index (Step 2)
# ─────────────────────────────────────────────────────────────────────────────

def compute_ensemble(norm_files: dict,
                     outdir: Path,
                     matrix_map: str,
                     exclude_matrix: list,
                     prefix: str,
                     output_merged: str):
    """
    Run eDNA_ensemble_index_pipeline.py on all normalized barcode tables.

    Parameters
    ----------
    norm_files     : {barcode: Path}  — normalized table paths per barcode
    outdir         : output directory for ensemble results
    matrix_map     : prefix→matrix mapping string, e.g. "B:Biofilm,M:Membrane,S:Sediment"
    exclude_matrix : list of "BARCODE:Matrix" strings
    prefix         : output filename prefix
    output_merged  : explicit filename for the final ensemble index table
    """
    separator()
    log("Step 2: Computing cross-barcode eDNA ensemble index ...", "STEP")
    separator()

    # Build --inputs argument: "BARCODE:path" pairs
    inputs_args = []
    for bc, fpath in norm_files.items():
        inputs_args.append(f"{bc}:{fpath}")

    script_args = ["--inputs"] + inputs_args
    script_args += ["--matrix-map", matrix_map]
    script_args += ["--outdir", str(outdir)]

    if exclude_matrix:
        script_args += ["--exclude-matrix"] + exclude_matrix

    if prefix:
        script_args += ["--prefix", prefix]

    if output_merged:
        script_args += ["--output-merged", output_merged]

    run_script(
        "eDNA_ensemble_index_pipeline.py",
        script_args,
        step_label="ensemble index computation"
    )


# ─────────────────────────────────────────────────────────────────────────────
# Argument parser
# ─────────────────────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_pipeline.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__
    )

    # ── Required ──────────────────────────────────────────────────────────
    parser.add_argument(
        "--barcodes", nargs="+", required=True,
        metavar="BARCODE",
        help=(
            "List of barcode identifiers to process. "
            "Must match input filename prefixes. "
            "Example: --barcodes 16S 18S COI 12S"
        )
    )
    parser.add_argument(
        "--otutab-dir", required=True, metavar="DIR",
        help=(
            "Directory containing raw OTU tables. "
            "Files must be named {BARCODE}_otutab.txt"
        )
    )
    parser.add_argument(
        "--taxonomy-dir", required=True, metavar="DIR",
        help=(
            "Directory containing taxonomy annotation files. "
            "Files must be named {BARCODE}_taxonomy.txt"
        )
    )
    parser.add_argument(
        "--workdir", required=True, metavar="DIR",
        help="Directory for intermediate per-barcode files (created if absent)."
    )
    parser.add_argument(
        "--outdir", required=True, metavar="DIR",
        help="Directory for final ensemble index outputs (created if absent)."
    )

    # ── Ensemble index options ────────────────────────────────────────────
    parser.add_argument(
        "--matrix-map", default="B:Biofilm,M:Membrane,S:Sediment",
        metavar="STR",
        help=(
            "Sample column prefix → matrix name mapping. "
            "Format: 'PREFIX1:Name1,PREFIX2:Name2,...' "
            "Default: 'B:Biofilm,M:Membrane,S:Sediment'"
        )
    )
    parser.add_argument(
        "--exclude-matrix", nargs="*", default=[], metavar="BARCODE:Matrix",
        help=(
            "Barcode–matrix combinations to exclude from normalization. "
            "Specify as BARCODE:Matrix pairs. "
            "Example: --exclude-matrix COI:Biofilm 12S:Biofilm"
        )
    )

    # ── Output naming ─────────────────────────────────────────────────────
    parser.add_argument(
        "--prefix", default="", metavar="STR",
        help="Uniform prefix added to all output filenames."
    )
    parser.add_argument(
        "--output-merged", default="", metavar="FILENAME",
        help=(
            "Explicit filename for the final ensemble index table. "
            "Overrides --prefix for this file. "
            "Example: --output-merged Yongding_ensemble_final.txt"
        )
    )

    # ── Control flags ─────────────────────────────────────────────────────
    parser.add_argument(
        "--skip-preprocessing", action="store_true",
        help=(
            "Skip Steps 1a–1c and go directly to ensemble index computation. "
            "Expects normalized files already present in --workdir, "
            "named {BARCODE}_normalized.txt."
        )
    )
    parser.add_argument(
        "--only-preprocessing", action="store_true",
        help=(
            "Run Steps 1a–1c only for all barcodes; skip ensemble index computation."
        )
    )

    return parser


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = build_parser()
    args = parser.parse_args()

    # ── Validate mutually exclusive flags ─────────────────────────────────
    if args.skip_preprocessing and args.only_preprocessing:
        log("--skip-preprocessing and --only-preprocessing cannot be used together.", "ERR")
        sys.exit(1)

    # ── Create output directories ─────────────────────────────────────────
    workdir  = Path(args.workdir)
    outdir   = Path(args.outdir)
    workdir.mkdir(parents=True, exist_ok=True)
    outdir.mkdir(parents=True, exist_ok=True)

    otutab_dir   = Path(args.otutab_dir)
    taxonomy_dir = Path(args.taxonomy_dir)

    separator("═")
    log("eDNA Multi-Barcode Ensemble Index Pipeline", "INFO")
    log(f"Barcodes       : {', '.join(args.barcodes)}", "INFO")
    log(f"OTU table dir  : {otutab_dir}", "INFO")
    log(f"Taxonomy dir   : {taxonomy_dir}", "INFO")
    log(f"Working dir    : {workdir}", "INFO")
    log(f"Output dir     : {outdir}", "INFO")
    log(f"Matrix map     : {args.matrix_map}", "INFO")
    log(f"Exclude matrix : {args.exclude_matrix if args.exclude_matrix else 'none'}", "INFO")
    log(f"Prefix         : '{args.prefix}'" if args.prefix else "Prefix         : (none)", "INFO")
    separator("═")

    # ── Per-barcode preprocessing ─────────────────────────────────────────
    norm_files = {}

    if not args.skip_preprocessing:
        for bc in args.barcodes:
            norm_path = preprocess_barcode(bc, otutab_dir, taxonomy_dir, workdir)
            norm_files[bc] = norm_path
    else:
        log("Skipping preprocessing (--skip-preprocessing flag set).", "WARN")
        for bc in args.barcodes:
            norm_path = workdir / f"{bc}_normalized.txt"
            if not norm_path.exists():
                log(f"Expected normalized file not found: {norm_path}", "ERR")
                sys.exit(1)
            norm_files[bc] = norm_path
            log(f"Using existing normalized file: {norm_path}", "INFO")

    # ── Ensemble index computation ────────────────────────────────────────
    if not args.only_preprocessing:
        compute_ensemble(
            norm_files=norm_files,
            outdir=outdir,
            matrix_map=args.matrix_map,
            exclude_matrix=args.exclude_matrix,
            prefix=args.prefix,
            output_merged=args.output_merged
        )
    else:
        log("Skipping ensemble index computation (--only-preprocessing flag set).", "WARN")

    # ── Summary ───────────────────────────────────────────────────────────
    separator("═")
    log("Pipeline finished successfully.", "OK")

    if not args.only_preprocessing:
        if args.output_merged:
            final_file = outdir / args.output_merged
        elif args.prefix:
            final_file = outdir / f"{args.prefix}_ensemble_index_merged.txt"
        else:
            final_file = outdir / "ensemble_index_merged.txt"
        log(f"Final ensemble index → {final_file}", "OK")

    separator("═")


if __name__ == "__main__":
    main()
