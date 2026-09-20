# Kinase Library

Kinase Library is a Python package for phosphoproteomics analyses based on kinase substrate preferences. It supports kinase prediction, binary enrichment, differential phosphorylation analysis, motif enrichment analysis (MEA), and visualization.

## Installation

```bash
pip install kinase-library
```

The package provides two command-line programs:

```bash
scan-motifs --help
run-mea --help
```

`scan-motifs` assigns phosphosite sequences to kinases from motif scores or percentiles. `run-mea` performs native GSEApy-based motif enrichment analysis and can write both a result table and the portable shared GSEA JSON representation.

## Python package

The Python package exposes the original Kinase Library analysis objects and functions:

```python
import kinase_library as kl
```

The repository contains worked notebooks for substrate scoring, binary enrichment, differential phosphorylation, and MEA.

## Project links

- [Kinase Library web tool](https://kinase-library.phosphosite.org)
- [Source repository](https://github.com/wolski/kinase-library)
- [PyPI package](https://pypi.org/project/kinase-library/)
