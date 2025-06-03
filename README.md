# `causalLFO`: R package for Causal Inference for Latent Factor-Modeled Outcomes

This R package provides all algorithms discussed in the paper "Causal Inference for Latent Factor-Modeled Outcomes". Code to reproduce results from our paper can be found in the [jennalandy/causalLFO_PAPER](https://github.com/jennalandy/causalLFO_PAPER/tree/master) repository.

## Installation

```{r}
remotes::install_github("jennalandy/causalLFO")

library(NMF)
library(causalLFO)
```

`NMF::nmf()` internally uses `setupLibPaths("NMF")`, which calls `path.package("NMF")`. This requires the NMF package to be attached, not just imported, so the user must library `NMF` as well as `causalLFO`.

Please install `NMF` if you have not yet done so. `NMF` requires the `Biobase` package, which may have to be installed separately from `Bioconductor`.

## Algorithms

This section follows the notation of the paper "Causal Inference for Latent Factor-Modeled Outcomes". Please refer to sections 2 and 3 for details.

The primary algorithm provided by this package is the `impute_and_stabilize` algorithm proposed in the paper "Causal Inference for Latent Factor-Modeled Outcomes". We first estimate unobserved potential outcomes $\mathbf Y(1-\mathbf T)$ with imputations $\tilde{\mathbf Y}_{1-\mathbf T}$. We fit the factor model $\hat{\boldsymbol \lambda}$ on a matrix of the original sample size $\tilde{\mathbf Y}_{\mathbf 0}$ as a combination of observed (for $T_i = 0$) and imputed (for $T_i = 1$) values to estimate $\ell_{\text{IS}, 0}(Y_i, T_i)$. A non-negative linear model is then used on the remaining $\tilde{\mathbf Y}_{1}$ with fixed $\hat{\boldsymbol \lambda}$ to measure $\ell_{\text{IS}, 1}(Y_i, T_i)$. Pairwise comparisons are used to estimate ATE with a mean of individual treatment effects estimator.

We also provide two baselines for comparison. First, `all_data` uses the full dataset $\mathbf Y$ to perform matrix decomposition to estimate $\hat{\boldsymbol \lambda}$ and $\ell_{\text{AD}, T_i}(Y_i, T_i)$, then again the full dataset to estimate causal effects with difference of means on $\ell_{\text{AD}, T_i}(Y_i, T_i)$. Second, we consider the `random_split` approach suggested by prior literature: identify a random subset of 50% of indices $S$, use $\mathbf Y_S$ as input to matrix decomposition to estimate $\hat{\boldsymbol \lambda}$, use a non-negative linear model on remaining data $\mathbf Y_{/S}$ to estimate latent outcomes $\ell_{\text{RS}, T_i}(Y_i, T_i)$ for $i \notin S$, and finally estimating the ATE with difference of means on those estimates.

Finally, we provide two ablations of the novel Impute and Stabilize algorithm to identify the relative benefits of each component. In the `impute` ablation, we perform imputation as before, but we perform matrix decomposition on the observed data $\mathbf Y$ to estimate $\ell_{\text{I}, T_i}(Y_i, T_i)$ for observed $T_i$ and non-negative linear model on the imputed $\tilde{\mathbf Y}_{1-\mathbf T}$ to estimate $\ell_{\text{I}, 1-T_i}(Y_i, T_i)$ for unobserved $T_i$. In the `stabilize` ablation, we perform matrix decomposition only on the untreated subset of observed data $\mathbf Y_{\{i:T_i = 0\}}$ and non-negative linear model on the treated subset $\mathbf Y_{\{i:T_i = 1\}}$, and a simple difference of means is used on the estimated latent outcomes $\ell_{\text{S}, T_i}(Y_i, T_i)$.

## To Do

-   algorithm roxygen documentation

-   example / quickstart
