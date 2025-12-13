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
"""Scan sequences against Kinase Library motifs to generate kinase-substrate assignments."""

from pathlib import Path
from typing import Literal

import sys

import cyclopts
import pandas as pd

from kinase_library.objects import phosphoproteomics as pps

app = cyclopts.App(
    name="scan_motifs",
    help="Scan sequences against Kinase Library motifs to generate kinase-substrate assignments.",
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
    seq_col: str = "Sequence",
    output: Path = Path("motif_matches.csv"),
    kin_type: Literal["ser_thr", "tyrosine"] = "ser_thr",
    method: Literal["percentile", "score"] = "percentile",
    threshold: float = 90.0,
    term2gene: bool = False,
) -> None:
    """Scan sequences against Kinase Library motifs.

    Parameters
    ----------
    input_file
        Path to input CSV/TSV/Excel file containing sequences.
    seq_col
        Column name containing the peptide sequences.
    output
        Path to save the output CSV file.
    kin_type
        Kinase type to use for scoring.
    method
        Scoring method to use.
    threshold
        Threshold for calling a match (90 for percentile, 1.0 for score).
    term2gene
        Format output as TERM2GENE (term=Kinase, gene=Sequence) for clusterProfiler.
    """
    # Load Data
    print(f"Loading data from {input_file}...")
    df = load_data(input_file)

    if seq_col not in df.columns:
        print(f"Error: Sequence column '{seq_col}' not found in input file.")
        print(f"Available columns: {list(df.columns)}")
        sys.exit(1)

    # Initialize PhosphoProteomics Object
    print("Initializing PhosphoProteomics object...")
    try:
        # Take unique sequences to score
        unique_seqs = df[[seq_col]].drop_duplicates().reset_index(drop=True)

        ppi = pps.PhosphoProteomics(
            data=unique_seqs,
            seq_col=seq_col,
            pp=True  # Assume standard phosphorylation notation or center
        )
    except Exception as e:
        print(f"Error initializing PhosphoProteomics object: {e}")
        sys.exit(1)

    # Calculate Scores/Percentiles
    print(f"Calculating {method}s (Type: {kin_type})...")

    try:
        if method == 'percentile':
            ppi.percentile(kin_type=kin_type, return_values=False)
            results_df = getattr(ppi, f"{kin_type}_percentiles")
        else:
            ppi.score(kin_type=kin_type, return_values=False)
            results_df = getattr(ppi, f"{kin_type}_scores")

    except Exception as e:
        print(f"Error during calculation: {e}")
        sys.exit(1)

    # Apply Threshold to find matches
    print(f"Applying threshold >= {threshold}...")

    # Filter using stack to get a Series, then filter
    stacked = results_df.stack()
    matches = stacked[stacked >= threshold]

    # Convert to DataFrame
    matches_df = matches.reset_index()
    matches_df.columns = [seq_col, 'Kinase', 'Value']

    # Format Output
    if term2gene:
        # TERM2GENE format: term (Kinase) first, then gene (Sequence)
        output_df = matches_df[['Kinase', seq_col]].copy()
        output_df.columns = ['term', 'gene']
    else:
        output_df = matches_df

    # Save
    print(f"Found {len(output_df)} assignments. Saving to {output}...")
    output_df.to_csv(output, index=False)
    print("Done.")


if __name__ == "__main__":
    app()
