# Chromatin alterations of cis-regulatory domains in human microglia in Alzheimer’s disease

- [Project Proposal](https://docs.google.com/document/d/14lGdC1wouelZSBIiexuzfM_xTdoa0PFkcUsDgmrrtlI/edit?tab=t.yje1p2c6e8i#heading=h.neg6n6u2fk6i)
- [Project Slides](https://docs.google.com/presentation/d/1DeebPCTVOYL7_xXkz5jO_dYZ5dXAhq5-tiSk8S1OhLw/edit?slide=id.h2c3801deca726a03_0_38#slide=id.h2c3801deca726a03_0_38)

Roussos Lab Fall 2026 Rotation Project

## Results

Here are all the compiled quarto notebooks detailing my results

1. [SV Results](./compiled.notebooks/01_Surrogate.Variable.Analysis.html)
1. [CRD Merging + Annotation](./compiled.notebooks/02_CRD.Investigation.html)
1. [Differential CRDs](./compiled.notebooks/03_Differential.CRD.Analysis.html)
1. [CRD Functional Enricment](./compiled.notebooks/04_CRD.Functional.Enrichment.Analysis.html)
1. [GWAS Results](./compiled.notebooks/05_GWAS.Analysis.html)
1. [CRD + GWAS + Functional Integration](./compiled.notebooks/06_CRD.GWAS.Integration.html)

## Reproducibility

I manage R dependencies with [rv](https://a2-ai.github.io/rv-docs/), I have included the `rv.lock` and `rproject.toml` here to list dependencies and allow for easy install.
```bash
$ rv sync  # install all dependencies specified in the .toml file
```
Run scripts to generate results
```bash
# Generate SVs and variancePartition() analysis of SVs
$ Rscript ./scripts/generate.SV.results.R
# Generate CRDs with SVs various regressed out
$ Rscript ./scripts/run.decorate.R
# Load CRDs, merge CRDs across meanClusterSizes + sample sets 
# Also annotated whether CRDs are AD/Control-specific or not
$ Rscript ./scripts/merge.and.annotate.CRDs.R
```
Compile quarto notebooks
```bash
$ find notebooks/*.qmd -type f -exec quarto render {} --to html \;
```

