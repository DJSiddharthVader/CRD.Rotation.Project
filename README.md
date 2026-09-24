# Chromatin alterations of cis-regulatory domains in human microglia in Alzheimer’s disease

[Project Proposal](https://docs.google.com/document/d/14lGdC1wouelZSBIiexuzfM_xTdoa0PFkcUsDgmrrtlI/edit?tab=t.yje1p2c6e8i#heading=h.neg6n6u2fk6i)

Roussos Lab Fall 2026 Rotation Project

## Reproducibility

I manage R dependencies with [rv](https://a2-ai.github.io/rv-docs/), I have included the `rv.lock` and `rproject.toml` here to list dependencies and allow for easy install.
```bash
rv sync  # install all dependencies specified in the .toml file
```

## TODO

- [x] Subset samples of peak residuals for downstream analyses
  - AD Only, Controls Only, No Others, All Samples
- [x] correlation of SVs estimated on different sample sets
  - Controls vs All on Control samples
  - AD vs All on AD samples
- [x] Generate SVs from peak residuals
- [ ] compare SV correlations across n.SVs i.e. same  matrix + CPM + sample set 
  - merge all SVs into a single matrix w/ n.SV_SV{1-N} for mega heatmap
    - or just do free space facets?
- [x] Compare SV ~ metadata correlations
  - All SVs vs AD Phenotypes + Age + Covariates
  - use `canCorPairs()` 
- [x] Implement code to residualize a set of SVs out of the peak residuals matrix
- [ ] SV varianceParrition analysis
  - [ ] get feature variances before residualizing SVs
  - [ ] get feature variances after residualizing SVs
  - [ ] compute pre/post SV-residualizing variance difference per feature
  - [ ] plot SV variances 
- [ ] Compare how LEF changes when including different numbers of SVs
  - Elbow plot of LEF vs #SVs included
  - for each number of SVs generated (n.SV), test incrementally adding SVs to the residualization
    - i.e. elbow plot for each n.SV + all other residual/SV sample sets
    - use different summary stats of LEFs across clusters per re-residualized matrix
- [ ] Compare cluster features as well
  - [ ] number of clusters
  - [ ] size distribution
  - [ ] LEF distribution
