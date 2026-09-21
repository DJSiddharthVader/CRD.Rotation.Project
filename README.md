# Chromatin alterations of cis-regulatory domains in human microglia in Alzheimer’s disease

[Project Proposal](https://docs.google.com/document/d/14lGdC1wouelZSBIiexuzfM_xTdoa0PFkcUsDgmrrtlI/edit?tab=t.yje1p2c6e8i#heading=h.neg6n6u2fk6i)

Roussos Lab Fall 2026 Rotation Project

## Reproducibility

I manage R dependencies with [rv](https://a2-ai.github.io/rv-docs/), I have included the `rv.lock` and `rproject.toml` here to list dependencies and allow for easy install.
```bash
rv sync  # install all dependencies specified in the .toml file
```

## TODO

- [ ] Subset samples of peak residuals for downstream analyses
  - AD Only, Controls Only, No Others, All Samples
- [ ] correlation of SVs estimated on different sample sets
  - Controls vs All on Control samples
  - AD vs All on AD samples
- [ ] Generate SVs from peak residuals
- [ ] Implement code to residualize a set of SVs out of the peak residuals matrix
- [ ] Compare how LEF changes when including different numbers of SVs
  - Elbow plot of LEF vs #SVs included
- [ ] Compare SV ~ metadata correlations
  - All SVs vs AD Phenotypes + Age + Covariates
  - use `canCorPairs()` 
