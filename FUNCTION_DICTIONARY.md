# GPTL Function Dictionary

This document is a map of the functions currently present in GPTL. It covers
the active package API, lower-level helpers, developing functions, archived
experiments, and the compiled C routines. The source path is included so that
the implementation can be inspected directly.

## Status conventions

- **Active/public**: part of the current package workflow and/or exported in
  `NAMESPACE`.
- **Active/internal**: used by active functions but not intended as the main
  user-facing workflow.
- **Developing**: under development and not part of the stable package
  layout.
- **Archived**: retained for historical reference or old analyses. New code
  should generally use the active equivalent where one exists.

The archived functions are documented for discoverability, not as a promise
that their old interfaces or results remain compatible with the active API.

## Active/public functions

| Function | Purpose | Main inputs / outputs | Source |
|---|---|---|---|
| `GD` | Fits target-population effects by gradient descent from sufficient statistics. | `XX`, `Xy`, optional prior effects `b`; returns effects or an iteration path. | [`R/GD.R`](R/GD.R) |
| `GDXy` | Convenience wrapper that accepts individual-level `X` and `y`, constructs `XX` and `Xy`, and calls `GD` or automatic RSS stopping. | `X`, `y`, preprocessing flags, and GD arguments; returns the corresponding GD result. | [`R/GDXy.R`](R/GDXy.R) |
| `GDXy2` | Direct individual-level gradient-descent implementation with optional early stopping. | `X`, `y`, optional `b`; returns final effects or an effect path. | [`R/GDXy2.R`](R/GDXy2.R) |
| `PR` | Fits penalized regression with shrinkage toward prior effects; `alpha=0` is ridge and `alpha=1` is lasso. | `XX`, `Xy`, `b`, `lambda`, `alpha`; returns penalized-effect paths and convergence information. | [`R/PR.R`](R/PR.R) |
| `BMM` | Gibbs sampler for Bayesian transfer learning with a finite mixture prior. | `XX`, `Xy`, prior matrix `B`, `my`, `vy`, `n`, and MCMC/prior controls; returns posterior effects and variance summaries. | [`R/BMM.R`](R/BMM.R) |
| `GD.CV` | Repeated random train/test evaluation of gradient descent using prediction correlation. | Individual-level `X`, `y`, split controls, and GD arguments; returns a correlation matrix over repetitions and iterations. | [`R/GD.CV.R`](R/GD.CV.R) |
| `GD.CV.ES` | Gradient descent with early stopping based on prediction correlation on a supplied test set. | `X_trn`, `y_trn`, `X_tst`, `y_tst`; returns the selected effect estimate and stopping iteration. | [`R/GD.CV.ES.R`](R/GD.CV.ES.R) |
| `GD.ES.CV.SS` | Early-stopping gradient descent using one train/test split while building sufficient statistics from the training data. | `X`, `y`, training indices `trn`, preprocessing and GD controls; returns fitted effects and stopping information. | [`R/GD.ES.CV.SS.R`](R/GD.ES.CV.SS.R) |
| `getSS` | Aligns LD, GWAS, and optional prior data by variant ID and constructs sufficient statistics. | `ld`, `gwas`, optional `B`; returns `XX`, `Xy`, effective `n`, and optionally aligned `B`. | [`R/getSS.R`](R/getSS.R) |
| `getCor` | Computes prediction correlation from sufficient statistics and one or more effect vectors. | `XX`, `Xy`, `yy`, `B`; returns one correlation per effect vector. | [`R/getCor.R`](R/getCor.R) |
| `fitLSYS` | Wrapper around the compiled Gauss-Seidel solver for a linear system. | `C`, `rhs`, initial `b`, active indices, `RSS`, `maxIter`, `tol`; returns updated `b` and `RSS`. | [`R/fitLSYS.R`](R/fitLSYS.R) |
| `sample_effects` | Wrapper for sampling marker effects in the Bayesian mixture implementation. | Effect-system inputs, mixture assignments, prior effects, and variance parameters. | [`R/sample_effects.R`](R/sample_effects.R) |
| `rMultinom` | Samples categorical mixture assignments from a matrix of non-negative probabilities. | `PROB`; returns one sampled component per column. | [`R/rMultinom.R`](R/rMultinom.R) |

## Active/internal functions

These functions are defined in active source files but are implementation
helpers rather than primary workflows.

| Function | Purpose | Source |
|---|---|---|
| `sample_effects_new` | Updated dense effect sampler that also maintains `RSS`. | [`R/sample_effects.R`](R/sample_effects.R) |
| `sample_effects_new_sparse` | Sparse-matrix version of `sample_effects_new`. | [`R/sample_effects.R`](R/sample_effects.R) |
| `rDirichlet` | Draws a probability vector from a Dirichlet distribution using gamma draws. | [`R/BMM.R`](R/BMM.R) |
| `which.first` | Returns the first matching index in a logical vector. | [`R/BMM.R`](R/BMM.R) |
| `sampleComp` | Older R implementation for sampling mixture components from probability rows. | [`R/BMM.R`](R/BMM.R) |
| `corBootstrap` | Estimates mean correlation over bootstrap resamples. | [`R/corBootstrap.R`](R/corBootstrap.R) |

## Developing functions

These functions are kept separately because their stopping criteria are still
being developed.

| Function | Purpose | Source |
|---|---|---|
| `GD.Auto.RSS` | Gradient descent with automatic stopping based on relative RSS change. | [`developing/GD.Auto.RSS.R`](developing/GD.Auto.RSS.R) |
| `GD.Auto.ErrVar` | Gradient descent with stopping diagnostics based on estimated error variance and its change. | [`developing/GD.Auto.ErrVar.R`](developing/GD.Auto.ErrVar.R) |

## Archived functions

### Bayesian mixture-model experiments

| Function | Historical purpose | Source |
|---|---|---|
| `BMM_Block` | Fits `BMM` separately over LD blocks and combines the results. | [`archive/BMM_Block.R`](archive/BMM_Block.R) |
| `get_block_ids` | Identifies contiguous blocks in a sparse matrix. | [`archive/BMM_Block.R`](archive/BMM_Block.R) |
| `BMM_SCALES` | Experimental Bayesian mixture model with prior-specific scales. | [`archive/BMM_SCALES.R`](archive/BMM_SCALES.R) |
| `linearIndex` | Converts a matrix row/column position to a linear index. | [`archive/BMM_SCALES.R`](archive/BMM_SCALES.R) |
| `getVarB` | Computes component-specific effect variances using prior scales and assignments. | [`archive/BMM_SCALES.R`](archive/BMM_SCALES.R) |
| `BMM.ld` | Earlier BMM workflow that reconstructs sufficient statistics directly from LD and GWAS inputs. | [`archive/BMM_archive.R`](archive/BMM_archive.R) |
| `BMM1` | Earlier sufficient-statistics BMM implementation. | [`archive/BMM_archive.R`](archive/BMM_archive.R) |
| `BMM_old` | Older Gibbs-sampling implementation using `C`, `rhs`, and `B0`. | [`archive/BMM_old.R`](archive/BMM_old.R) |
| `BMM_old2` | Intermediate Gibbs-sampling implementation with thinning and mixture controls. | [`archive/BMM_old2.R`](archive/BMM_old2.R) |
| `rDirichlet` | Archived Dirichlet sampler used by older BMM implementations. | [`archive/BMM_archive.R`](archive/BMM_archive.R), [`archive/BMM_old.R`](archive/BMM_old.R), [`archive/BMM_old2.R`](archive/BMM_old2.R) |
| `which.first` | Archived helper for selecting the first matching component. | [`archive/BMM_archive.R`](archive/BMM_archive.R), [`archive/BMM_old2.R`](archive/BMM_old2.R) |
| `sampleComp` | Archived R mixture-component sampler. | [`archive/BMM_archive.R`](archive/BMM_archive.R), [`archive/BMM_old2.R`](archive/BMM_old2.R) |

### Gradient-descent experiments

| Function | Historical purpose | Source |
|---|---|---|
| `GD.CV` | Earlier fold-based cross-validation implementation with configurable accuracy function and lambda values. | [`archive/GD.CV.R`](archive/GD.CV.R) |
| `GD.Full` | Earlier full gradient-descent implementation using sufficient statistics. | [`archive/GD.Full.R`](archive/GD.Full.R) |
| `GD.ld` | Earlier gradient descent workflow starting from LD and GWAS data. | [`archive/GD_archive.R`](archive/GD_archive.R) |
| `GD.R` | Earlier dense gradient-descent implementation. | [`archive/GD_archive.R`](archive/GD_archive.R) |
| `GD0` | Earlier gradient-descent variant with explicit prior and shrinkage controls. | [`archive/GD_archive.R`](archive/GD_archive.R) |
| `GD1` | Earlier gradient-descent variant. | [`archive/GD_archive.R`](archive/GD_archive.R) |
| `GD_sparse` | Earlier sparse gradient-descent variant. | [`archive/GD_archive.R`](archive/GD_archive.R) |

### Penalized-regression and solver experiments

| Function | Historical purpose | Source |
|---|---|---|
| `PR_TEST` | Experimental penalized-regression implementation. | [`archive/PR_TEST.R`](archive/PR_TEST.R) |
| `PR_lambdah2` | Earlier penalized-regression implementation over a lambda grid. | [`archive/PR_archive.R`](archive/PR_archive.R) |
| `PR.ld` | Earlier penalized-regression workflow starting from LD and GWAS data. | [`archive/PR_archive.R`](archive/PR_archive.R) |
| `PR1` | Earlier sufficient-statistics penalized-regression implementation. | [`archive/PR_archive.R`](archive/PR_archive.R) |
| `LASSO.CD1` | Early coordinate-descent lasso implementation. | [`archive/LASSO.R`](archive/LASSO.R) |
| `LASSO` | Earlier lasso fitting implementation. | [`archive/LASSO.R`](archive/LASSO.R) |
| `LASSO.CD` | Coordinate-descent lasso implementation. | [`archive/LASSO.R`](archive/LASSO.R) |
| `LASSO.GS2` | Gauss-Seidel-style lasso implementation. | [`archive/LASSO.R`](archive/LASSO.R) |
| `OLS` | Ordinary least-squares solver based on sufficient statistics. | [`archive/LASSO.R`](archive/LASSO.R) |
| `RIDGE.CD` | Coordinate-descent ridge implementation. | [`archive/RIDGE_CD.R`](archive/RIDGE_CD.R) |
| `RR` | Earlier ridge-regression wrapper with optional shrinkage toward `b0`. | [`archive/RR.R`](archive/RR.R) |
| `WSS` | Pools weighted sufficient statistics from multiple studies or groups. | [`archive/WSS.R`](archive/WSS.R) |

### Miscellaneous utilities

| Function | Historical purpose | Source |
|---|---|---|
| `nextBlock` | Finds the next block boundary while allowing a specified number of gaps. | [`misc/getBlock.R`](misc/getBlock.R) |
| `findBlocks` | Finds candidate blocks from an LD or relationship matrix using a threshold and maximum gap. | [`misc/getBlock.R`](misc/getBlock.R) |
| `mergeBlocks` | Merges or filters candidate blocks according to a minimum block size. | [`misc/getBlock.R`](misc/getBlock.R) |

## Archived/developing parameter variants

Some archived functions use older names for quantities that have active names
in the parameter dictionary. The most important examples are:

| Historical name | Current comparable name | Note |
|---|---|---|
| `learning_rate` | `learningRate` | Older gradient-descent functions use underscore notation. |
| `conv_threshold` | `convThreshold` | Older penalized-regression functions use underscore notation. |
| `C` / `rhs` | `XX` / `Xy` | Older BMM and solver code names the sufficient statistics as a generic system and right-hand side. |
| `B0` | `B` | Older BMM code uses `B0` for the prior-effect matrix. |
| `lambda0` | — | Legacy parameter controlling additional shrinkage toward `b0`; it is not part of the main active API. |
| `nFolds` | `nTst` / `trn` | Archived cross-validation uses folds, while active workflows use holdout counts or explicit training indices. |

## Compiled routines

The following routines are implemented in `src/` and called through R
wrappers. They are not ordinary user-facing R functions.

| Routine | Purpose | Source |
|---|---|---|
| `GRAD_DESC` | Dense gradient-descent update. | [`src/GRAD_DESC.c`](src/GRAD_DESC.c) |
| `GRAD_DESC_sparse` | Sparse gradient-descent update. | [`src/GRAD_DESC.c`](src/GRAD_DESC.c) |
| `GRAD_DESC_Xy` | Direct `X`/`y` gradient-descent update using residuals and column sums of squares. | [`src/GRAD_DESC.c`](src/GRAD_DESC.c) |
| `ElasticNet` | Dense elastic-net coordinate-descent update. | [`src/ElasticNet.c`](src/ElasticNet.c) |
| `ElasticNet_sparse` | Sparse elastic-net coordinate-descent update. | [`src/ElasticNet.c`](src/ElasticNet.c) |
| `LASSO_CD` | Lasso coordinate-descent update. | [`src/LASSO_CD.c`](src/LASSO_CD.c) |
| `RIDGE_CD` | Ridge coordinate-descent update. | [`src/RIDGE_CD.c`](src/RIDGE_CD.c) |
| `fitLSYS` | Gauss-Seidel linear-system solver. | [`src/fitLSYS.c`](src/fitLSYS.c) |
| `rMultinomial` | Multinomial/categorical sampler used by `rMultinom`. | [`src/sample_multinomial.c`](src/sample_multinomial.c) |
| `sample_effects` | Original dense Bayesian effect sampler. | [`src/sample_effects.c`](src/sample_effects.c) |
| `sample_effects_new` | Dense Bayesian effect sampler that updates residual sum of squares. | [`src/sample_effects_new.c`](src/sample_effects_new.c) |
| `sample_effects_new_sparse` | Sparse Bayesian effect sampler. | [`src/sample_effects_new.c`](src/sample_effects_new.c) |

## Function documentation workflow

When adding or changing a function:

1. Keep the function in the directory that matches its status: `R/` for
   active code, `developing/` for work in progress, and `archive/` for
   historical code.
2. Add or update its entry here, including its purpose, source path, and
   status.
3. Keep parameter names consistent with
   [`PARAMETERS_DICTIONARY.md`](PARAMETERS_DICTIONARY.md).
4. For active/public functions, update the matching `.Rd` documentation and
   `NAMESPACE` as appropriate.
