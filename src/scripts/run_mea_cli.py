#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "numpy",
#   "pandas>=2.0",
#   "requests",
#   "pyarrow",
#   "tqdm",
#   "natsort",
#   "scipy",
#   "cyclopts"
# ]
# ///
"""Run Kinase Library Motif Enrichment Analysis (MEA) from a rank file."""

from pathlib import Path
from typing import Literal

import sys

import cyclopts
import pandas as pd

from kinase_library.enrichment import mea

app = cyclopts.App(
    name="run_mea",
    help="Run Kinase Library Motif Enrichment Analysis (MEA) from a rank file.",
)


def load_data(input_file: Path) -> pd.DataFrame:
    """Load data from CSV, TSV, Excel, or RNK file."""
    file_ext = input_file.suffix.lower()

    if file_ext == '.csv':
        return pd.read_csv(input_file)
    elif file_ext in ['.tsv', '.txt', '.rnk']:
        return pd.read_csv(input_file, sep='\t')
    elif file_ext in ['.xlsx', '.xls']:
        return pd.read_excel(input_file)
    else:
        print("Unsupported file extension. Trying CSV...")
        return pd.read_csv(input_file)


@app.default
def main(
    input_file: Path,
    *,
    rank_col: str,
    seq_col: str = "Sequence",
    output: Path = Path("mea_results.csv"),
    kin_type: Literal["ser_thr", "tyrosine"] = "ser_thr",
    method: Literal["percentile", "score"] = "percentile",
    threshold: float = 90.0,
    permutations: int = 1000,
    threads: int = 4,
) -> None:
    """Run Kinase Library Motif Enrichment Analysis.

    Parameters
    ----------
    input_file
        Path to input CSV/TSV/Excel file containing sequences and ranks.
    rank_col
        Column name containing the ranking statistic (e.g., t-statistic, fold change).
    seq_col
        Column name containing the peptide sequences.
    output
        Path to save the output CSV file.
    kin_type
        Kinase type to use for analysis.
    method
        Scoring method to use.
    threshold
        Threshold for the scoring method (90 for percentile).
    permutations
        Number of permutations for GSEA.
    threads
        Number of threads to use.
    """
    # Load Data
    print(f"Loading data from {input_file}...")
    data = load_data(input_file)

    if rank_col not in data.columns:
        print(f"Error: Rank column '{rank_col}' not found in input file.")
        print(f"Available columns: {list(data.columns)}")
        sys.exit(1)

    if seq_col not in data.columns:
        print(f"Error: Sequence column '{seq_col}' not found in input file.")
        print(f"Available columns: {list(data.columns)}")
        sys.exit(1)

    # Initialize RankedPhosData
    print("Initializing RankedPhosData...")

    try:
        rpd = mea.RankedPhosData(
            dp_data=data,
            rank_col=rank_col,
            seq_col=seq_col,
            pp=True  # Assuming data derived from phosphoproteomics usually indicates the p-site
        )
    except Exception as e:
        print(f"Error initializing data structure: {e}")
        sys.exit(1)

    # Run MEA
    print(f"Running MEA ({kin_type}, {method} > {threshold})...")
    try:
        results = rpd.mea(
            kin_type=kin_type,
            kl_method=method,
            kl_thresh=threshold,
            permutation_num=permutations,
            threads=threads
        )
    except Exception as e:
        print(f"Error running MEA: {e}")
        sys.exit(1)

    # Save Results
    print(f"Saving results to {output}...")
    results.enrichment_results.to_csv(output)
    print("Done.")


if __name__ == "__main__":
    app()
