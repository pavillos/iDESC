
<!-- README.md is generated from README.Rmd. Please edit that file -->

# iDESC

<!-- badges: start -->
<!-- badges: end -->

The goal of iDESC is to identify differentially expressed genes between groups in single cell RNA-seq data with several subjects per group, using a negative binomial mixed model with the subject as random effect.

This is a fork of [yl883/iDESC](https://github.com/yl883/iDESC). Version 0.3.0 changes the model and the output (see `NEWS.md`); version 0.2.0 is commit `2be58b9`.

## Installation

``` r
devtools::install_github("pavillos/iDESC")
# version 0.2.0, to reproduce earlier analyses
devtools::install_github("pavillos/iDESC", ref = "2be58b9")
```

## Example

``` r
library(iDESC)
data(IPF_example)
mat=IPF_example$mat
meta=IPF_example$meta
sequencing_depth=IPF_example$meta$sequencing_depth
str(meta)
result=iDESC(mat,meta,subject_var="subject",group_var="disease",
             norm_opt="User",user_sf = sequencing_depth,cores=4)

# one row per gene: estimates, Wald tests, likelihood-ratio test of the group effect
head(result$model_results)
# what happened to every input gene
table(result$gene_status$status)
```

`iDESC()` fits, for every gene, `y ~ group + covariates + offset(log(norm_sf)) + (1 | subject)` with `glmmTMB` (`zi_model = "none"`; `zi_model = "constant"` adds a zero-inflation intercept) and reports

* the likelihood-ratio test of the group effect: `LRT_stat`, `Pval_LRT_group`, `FDR_LRT_group`;
* the Wald test of every coefficient (`Beta_*`, `Pval_Beta_*`, `FDR_Beta_*`) and, for three or more groups, of the remaining pairwise contrasts (`Beta_<B>_vs_<A>`);
* a `Status` per gene and the table `gene_status`.

With few subjects per group the asymptotic p-values are too small. By default (`permutation = "auto"`) permutation p-values are added when the smallest group has fewer than 10 subjects: the group label is permuted between subjects and the columns `Pval_perm_group`, `FDR_perm_group`, `Pval_perm_Beta_*`, `FDR_perm_Beta_*` and `N_perm` are added.

``` r
result=iDESC(mat,meta,subject_var="subject",group_var="disease",
             norm_opt="User",user_sf = sequencing_depth,cores=4,
             permutation="yes",n_permutations=999,seed=1)
```
