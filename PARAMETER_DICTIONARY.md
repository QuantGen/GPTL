# GPTL Parameter Dictionary

This document defines the canonical parameter names used throughout GPTL. New
functions should reuse these names whenever they represent the same quantity.
Existing function names and parameter names are not changed by this document.

## Naming rules

1. **Mathematical parameters keep their mathematical notation.** Use the exact
   letters, including capitalization, when the capitalization carries meaning:
   `X`, `y`, `XX`, `Xy`, `yy`, `B`, `b`, `R2`, `C`, and `RSS`.
2. **Single-word parameters are lower case:** `thin`, `verbose`, `lambda`,
   `seed`, and `tol`.
3. **Multi-word parameters use lower camel case.** The first word is lower
   case and each following word starts with a capital letter: `nIter`,
   `priorProb`, `learningRate`, and `returnPath`.
4. **Parameters with a mathematical subscript or component are separated by a
   period:** `df0.E`, `S0.E`, and `df0.b`. The period separates the base
   quantity from its component; the component capitalization should follow its
   mathematical meaning.

Do not introduce alternate spellings such as `niter`, `learning_rate`,
`prior_prob`, or `return_path` for these established parameters.

## Core data and sufficient statistics

| Parameter | Meaning | Typical form | Used by |
|---|---|---|---|
| `X` | Target-population genotype/design matrix. Rows are individuals and columns are variants. | Numeric matrix, `n` × `p` | `GDXy`, `GDXy2`, `GD.CV`, `GD.CV.ES`, `GD.ES.CV.SS` |
| `y` | Target-population phenotype/response vector. | Numeric vector of length `n` | Functions that take individual-level data |
| `XX` | Cross-product matrix `X'X`. | `p` × `p` matrix or sparse matrix | `GD`, `PR`, `BMM`, `getCor` |
| `Xy` | Cross-product vector `X'y`. | Numeric vector or one-column matrix of length `p` | `GD`, `PR`, `BMM`, `getCor` |
| `yy` | Phenotype cross-product `y'y`. | Numeric scalar | `getCor` |
| `B` | Matrix of prior effect estimates. Each row is a variant and each column is a prior source or mixture component. | `p` × `k` matrix | `BMM`, `BMM_Block`, `getSS`, `getCor` |
| `b` | Current or initial target-population effect estimates. | Numeric vector or one-column matrix of length `p` | `GD`, `PR`, gradient-descent helpers |

Variant IDs should be stored as row names and column names where applicable.
`XX`, `Xy`, and `B` are aligned by these IDs before fitting.

## General model and optimization parameters

| Parameter | Meaning | Typical form | Used by |
|---|---|---|---|
| `n` | Target-population sample size used by the model. | Positive integer or effective sample size | `BMM`, `GD.Auto.ErrVar`, `getSS` output |
| `my` | Mean of the target phenotype. | Numeric scalar | `BMM`, `BMM_Block` |
| `vy` | Variance of the target phenotype. | Positive numeric scalar | `BMM`, `BMM_Block` |
| `nIter` | Number of optimization or sampling iterations. | Positive integer | `GD`, `BMM`, `PR`, `GDXy2` |
| `maxIter` | Maximum number of iterations before an algorithm stops. | Positive integer | Automatic stopping and coordinate-descent functions |
| `learningRate` | Step size used by gradient descent. | Positive numeric scalar | `GD`, `GDXy`, `GD.CV.ES`, automatic GD helpers |
| `lambda` | Shrinkage or ridge-penalty parameter. | Numeric scalar or vector | `GD`, `PR`, GD helpers |
| `nLambda` | Number of penalty values to generate when `lambda` is not supplied. | Positive integer | `PR` |
| `alpha` | Elastic-net mixing parameter: `0` is ridge and `1` is lasso. | Numeric value in `[0, 1]` | `PR` |
| `convThreshold` | Convergence tolerance for iterative coordinate descent. | Small positive numeric scalar | `PR` |
| `returnPath` | Whether to return the full sequence of effect estimates across iterations. | Logical | `GD`, `PR`, `GDXy2` |
| `verbose` | Whether to print progress or diagnostic messages. | Logical | Most fitting functions |
| `tol` | Numerical convergence tolerance for a lower-level solver. | Small positive numeric scalar | `fitLSYS` |

## Bayesian mixture parameters

| Parameter | Meaning | Typical form | Used by |
|---|---|---|---|
| `burnIn` | Number of initial MCMC samples discarded before posterior summaries are accumulated. | Non-negative integer | `BMM`, `BMM_Block` |
| `thin` | Number of MCMC iterations between retained samples. | Positive integer | `BMM`, `BMM_Block` |
| `R2` | Prior expected proportion of phenotype variance explained by the regression. | Numeric value between `0` and `1` | `BMM`, `BMM_Block` |
| `nComp` | Number of prior or mixture components. | Positive integer | `BMM`, `BMM_Block` |
| `K` | Inverse of `nComp`, used in the prior variance calculation. | Numeric scalar or vector | `BMM`, `BMM_Block` |
| `df0.E` | Degrees of freedom for the scaled inverse-chi-squared prior on residual variance. | Positive numeric scalar | `BMM`, `BMM_Block` |
| `S0.E` | Scale parameter for the residual-variance prior. | Positive numeric scalar | `BMM`, `BMM_Block` |
| `df0.b` | Degrees of freedom for the prior on effect variance for each mixture component. | Numeric vector of length `nComp` | `BMM`, `BMM_Block` |
| `priorProb` | Prior probabilities of the mixture components. | Numeric vector of length `nComp` | `BMM`, `BMM_Block` |
| `priorCounts` | Prior pseudo-counts for the mixture probabilities. | Numeric vector of length `nComp` | `BMM`, `BMM_Block` |
| `fixVarE` | Whether to fix the residual variance instead of sampling it. | Logical | `BMM`, `BMM_Block` |
| `fixVarB` | Whether to fix the effect variance for each mixture component. | Logical vector of length `nComp` | `BMM`, `BMM_Block` |

## Preprocessing, validation, and cross-validation parameters

| Parameter | Meaning | Typical form | Used by |
|---|---|---|---|
| `centerX` | Whether to center the columns of `X`. | Logical | `GDXy`, `GD.CV.ES`, `GD.ES.CV.SS` |
| `scaleX` | Whether to scale the columns of `X`. | Logical | `GDXy`, `GD.CV.ES`, `GD.ES.CV.SS` |
| `earlyStop` | Whether to stop gradient descent based on an improvement criterion before `maxIter` or `nIter`. | Logical | `GDXy`, `GDXy2` |
| `pctChangeRSS` | Relative change in residual sum of squares used as an early-stopping threshold. | Non-negative numeric scalar | `GDXy`, `GDXy2`, `GD.Auto.RSS` |
| `nTst` | Number of observations assigned to the test set in repeated cross-validation. | Positive integer | `GD.CV` |
| `nRep` | Number of repeated train/test splits. | Positive integer | `GD.CV` |
| `seed` | Random-number seed used to make resampling reproducible. | Integer or `NULL` | `GD.CV` |
| `trn` | Integer indices identifying the training observations. | Integer vector | `GD.ES.CV.SS` |
| `X_trn` | Training design matrix. | Numeric matrix | `GD.CV.ES` |
| `y_trn` | Training phenotype vector. | Numeric vector | `GD.CV.ES` |
| `X_tst` | Test design matrix. | Numeric matrix | `GD.CV.ES` |
| `y_tst` | Test phenotype vector. | Numeric vector | `GD.CV.ES` |

## Lower-level computational parameters

These names are used by wrappers around the compiled routines and are mainly
relevant to advanced users.

| Parameter | Meaning |
|---|---|
| `C` | Coefficient or system matrix passed to a lower-level solver. |
| `rhs` | Right-hand-side vector or matrix for a linear system. |
| `active` | Indices of the active effects or variables. |
| `RSS` | Residual sum of squares, either current or updated. |
| `B0` | Prior effect matrix passed to the effect sampler. |
| `varE` | Residual/error variance. |
| `varB` | Effect variance, generally one value per prior component. |
| `d` | Mixture-component assignment for each variant. |
| `PROB` | Matrix of non-negative sampling probabilities. |
| `times` | Number of bootstrap replicates. |
| `x`, `y` | Generic paired vectors used by correlation utilities; in model-fitting functions, prefer the more specific `X` and `y`. |

## Adding a new parameter

Before introducing a new parameter name:

1. Check this dictionary for an existing name with the same meaning.
2. Reuse the existing spelling and capitalization if the meaning is the same.
3. Use mathematical notation for mathematical quantities and lower camel case
   for multi-word non-mathematical quantities.
4. Update this dictionary and the relevant `.Rd` documentation in the same
   change.
