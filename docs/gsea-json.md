# Shared GSEA JSON

Native MEA results can be serialized to the portable shared GSEA JSON document used by the PTM pipeline. The document stores ordinary JSON values only: complete result rows, the ordered ranking, kinase membership, leading-edge members, running-score traces, identifiers, and analysis parameters.

## Python API

After running MEA, write one contrast with `write_gsea_result_json`:

```python
results.write_gsea_result_json("mea_result.json", "treated_vs_control")
```

Use `to_gsea_result_document` when the caller needs the document in memory:

```python
document = results.to_gsea_result_document("treated_vs_control")
```

The document identifies its method as `MEA`, its category as `KinaseLib`, and its backend as `gseapy`. It uses the same `gsea_result` block as other GSEA producers, so downstream readers do not need a Kinase Library-specific serialization.

## Command line

`run-mea` writes the table selected by `--output` and optionally writes the JSON document selected by `--json-output`:

```bash
run-mea ranked-sites.tsv \
  --rank-col statistic \
  --seq-col Sequence \
  --output mea_treated_vs_control.csv \
  --json-output mea_treated_vs_control.json \
  --contrast treated_vs_control
```

The stored running-score arrays come from the same GSEApy result used for the MEA statistics. They can therefore reproduce enrichment score versus ranked-list plots without rerunning the enrichment analysis.
