# configs/

Path and dataset configuration for reproducing the xGATE manuscript analyses.

## Usage

```bash
cp configs/paths.example.yaml    configs/paths.yaml
cp configs/datasets.example.yaml configs/datasets.yaml
# then edit the copies for your machine
```

`configs/paths.yaml` and `configs/datasets.yaml` are **git-ignored** — they hold
machine-specific absolute paths and should never be committed. Only the
`*.example.yaml` templates are tracked.

## What each file holds

| File | Purpose |
|------|---------|
| `paths.example.yaml` | Repo root, data root, results/figure output roots, optional Python environment prefix for the wrapper. |
| `datasets.example.yaml` | Per-dataset accessions and expected input file names. |

## Current limitation (author review)

The analysis scripts and notebooks still contain a **placeholder absolute root**
`/path/to/xGATE` (introduced when private cluster paths were removed for public
release). Until full config-driven path injection is added, either:

- edit `/path/to/xGATE` in the scripts you run to point at your clone, or
- set the data root to match and run from the repository root.

This is tracked in `docs/files_needing_author_review.md`.
