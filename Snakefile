"""
Snakefile for running kinase-library CLI tools on example data.

Usage:
    snakemake --cores 4           # Run all targets
    snakemake scan_motifs --cores 1   # Run only motif scanning
    snakemake mea --cores 4       # Run only MEA analysis
"""

# Configuration
INPUT_DIR = "data/example_inputs"
OUTPUT_DIR = "data/example_outputs"
INPUT_ZIP = "data/example_inputs.zip"

# Get all rank files for MEA (after unzip)
CONTRASTS = ["dilstein_InhibitorA_vs_Vehicle", "dilstein_InhibitorB_vs_Vehicle",
             "dilstein_InhibitorC_vs_Vehicle", "dilstein_InhibitorD_vs_Vehicle"]


THRESHOLDS = [90, 95]


# =============================================================================
# Unzip Rule
# =============================================================================

rule unzip_inputs:
    input:
        INPUT_ZIP,
    output:
        directory(INPUT_DIR),
    shell:
        """
        unzip -o {input} -d data/
        """


rule all:
    input:
        # Motif scanning outputs (90th percentile)
        f"{OUTPUT_DIR}/motif_matches.csv",
        f"{OUTPUT_DIR}/term2gene.csv",
        # Motif scanning outputs (95th percentile)
        f"{OUTPUT_DIR}/motif_matches_95.csv",
        f"{OUTPUT_DIR}/term2gene_95.csv",
        # MEA outputs for each contrast (90th and 95th percentile)
        expand(f"{OUTPUT_DIR}/mea_{{contrast}}.csv", contrast=CONTRASTS),
        expand(f"{OUTPUT_DIR}/mea_{{contrast}}_95.csv", contrast=CONTRASTS),


# =============================================================================
# Motif Scanning Rules
# =============================================================================

rule scan_motifs:
    input:
        f"{OUTPUT_DIR}/motif_matches.csv",
        f"{OUTPUT_DIR}/term2gene.csv",
        f"{OUTPUT_DIR}/motif_matches_95.csv",
        f"{OUTPUT_DIR}/term2gene_95.csv",


rule scan_motifs_default:
    input:
        seqs=f"{INPUT_DIR}/seqwindows_unique.tsv",
        _unzip=INPUT_DIR,
    output:
        f"{OUTPUT_DIR}/motif_matches.csv",
    shell:
        """
        uv run python src/scripts/scan_motifs_cli.py \
            {input.seqs} \
            --output {output} \
            --kin-type ser_thr \
            --method percentile \
            --threshold 90
        """


rule scan_motifs_term2gene:
    input:
        seqs=f"{INPUT_DIR}/seqwindows_unique.tsv",
        _unzip=INPUT_DIR,
    output:
        f"{OUTPUT_DIR}/term2gene.csv",
    shell:
        """
        uv run python src/scripts/scan_motifs_cli.py \
            {input.seqs} \
            --output {output} \
            --kin-type ser_thr \
            --method percentile \
            --threshold 90 \
            --term2gene
        """


rule scan_motifs_default_95:
    input:
        seqs=f"{INPUT_DIR}/seqwindows_unique.tsv",
        _unzip=INPUT_DIR,
    output:
        f"{OUTPUT_DIR}/motif_matches_95.csv",
    shell:
        """
        uv run python src/scripts/scan_motifs_cli.py \
            {input.seqs} \
            --output {output} \
            --kin-type ser_thr \
            --method percentile \
            --threshold 95
        """


rule scan_motifs_term2gene_95:
    input:
        seqs=f"{INPUT_DIR}/seqwindows_unique.tsv",
        _unzip=INPUT_DIR,
    output:
        f"{OUTPUT_DIR}/term2gene_95.csv",
    shell:
        """
        uv run python src/scripts/scan_motifs_cli.py \
            {input.seqs} \
            --output {output} \
            --kin-type ser_thr \
            --method percentile \
            --threshold 95 \
            --term2gene
        """


# =============================================================================
# MEA (Motif Enrichment Analysis) Rules
# =============================================================================

rule mea:
    input:
        expand(f"{OUTPUT_DIR}/mea_{{contrast}}.csv", contrast=CONTRASTS),
        expand(f"{OUTPUT_DIR}/mea_{{contrast}}_95.csv", contrast=CONTRASTS),


rule run_mea:
    input:
        rnk=f"{INPUT_DIR}/{{contrast}}.rnk",
        _unzip=INPUT_DIR,
    output:
        f"{OUTPUT_DIR}/mea_{{contrast}}.csv",
    threads: 4
    shell:
        """
        uv run python src/scripts/run_mea_cli.py \
            {input.rnk} \
            --rank-col Score \
            --output {output} \
            --kin-type ser_thr \
            --method percentile \
            --threshold 90 \
            --permutations 1000 \
            --threads {threads}
        """


rule run_mea_95:
    input:
        rnk=f"{INPUT_DIR}/{{contrast}}.rnk",
        _unzip=INPUT_DIR,
    output:
        f"{OUTPUT_DIR}/mea_{{contrast}}_95.csv",
    threads: 4
    shell:
        """
        uv run python src/scripts/run_mea_cli.py \
            {input.rnk} \
            --rank-col Score \
            --output {output} \
            --kin-type ser_thr \
            --method percentile \
            --threshold 95 \
            --permutations 1000 \
            --threads {threads}
        """


# =============================================================================
# Utility Rules
# =============================================================================

rule clean:
    shell:
        f"rm -rf {OUTPUT_DIR}/*.csv"


rule help:
    shell:
        """
        echo "Available targets:"
        echo "  all          - Run all analyses (default)"
        echo "  scan_motifs  - Run motif scanning on unique sequences"
        echo "  mea          - Run MEA on all contrast rank files"
        echo "  clean        - Remove all output files"
        echo ""
        echo "Example usage:"
        echo "  snakemake --cores 4"
        echo "  snakemake mea --cores 4"
        """
