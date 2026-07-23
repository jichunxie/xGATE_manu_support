# Figure 2 - pancreas: topology, marginal expression, embeddings

Pancreatic beta-cell analysis showing that co-expression **graph topology** (not
edge density or marginal expression alone) separates active from inactive
pathways.

**Canonical data source:** the *updated* pancreas dataset, public accession
**GSE148073** (the committed adjacency matrices under `data/` derive from it).
Run the notebooks in your local Jupyter environment; run Figure 2b in R (Seurat).

## Manuscript items

| Panel | Backend script/notebook | Required inputs | Expected output |
|-------|-------------------------|-----------------|-----------------|
| 2a co-expression graphs | `analysis/figure2_pancreas/fig2a_coexpr_graphs.ipynb` | `data/adj_matrix_pancreas_ctrl_final.csv`, pathway gene sets | graph panels for `Figure2.jpeg` |
| 2b marginal-expression boxplots | `analysis/figure2_pancreas/fig2b_Mean_Exp.Rmd` (render to HTML locally) | `data/pancreas_human.rds` (GSE148073 Seurat object), KEGG pathway genes | boxplot JPGs for `Figure2.jpeg` |
| 2c MDS embedding space | `analysis/figure2_pancreas/fig2c_embedding_comparison.ipynb` | `Pancreas/Embeddings/*_xGATE_pancreas_embeddings.csv` plus FeatherGraph/Graph2Vec/NetLSD embeddings | MDS scatter for `Figure2.jpeg` |

## Figure 2b

`analysis/figure2_pancreas/fig2b_Mean_Exp.Rmd` is the R/Seurat notebook that produces the
marginal-expression boxplots. Render it locally to read the analysis and expected output:

```bash
Rscript -e 'rmarkdown::render("analysis/figure2_pancreas/fig2b_Mean_Exp.Rmd")'
```

(The rendered HTML is a generated output and is not committed — see `.gitignore`.)

It loads the GSE148073 pancreas Seurat object (`data/pancreas_human.rds`), selects control
`type B pancreatic cell`s, normalizes with `SCTransform(vst.flavor="v2")`, pulls pathway genes
from KEGG via `KEGGREST` / `org.Hs.eg.db`, and draws one boxplot per pathway with one point per
pathway gene. The pathways match Figure 2a (e.g. `AMPK signaling pathway`,
`Bacterial invasion of epithelial cells`).

Requires R with `Seurat`, `dplyr`, `tidyverse`, `KEGGREST`, and `org.Hs.eg.db`
(see `docs/r_environment.md`). Render with:

```bash
Rscript -e 'rmarkdown::render("analysis/figure2_pancreas/fig2b_Mean_Exp.Rmd")'
```

after placing the GSE148073 Seurat object at `data/pancreas_human.rds`. Output boxplot JPGs are
written to `figures/` (git-ignored) and assembled into `Figure2.jpeg` manually.

## Other Figure 2 panels

```bash
jupyter lab
# open:
# analysis/figure2_pancreas/fig2a_coexpr_graphs.ipynb
# analysis/figure2_pancreas/fig2c_embedding_comparison.ipynb
```

## Notes

- Embedding CSVs are large and git-ignored; regenerate from the notebook.
