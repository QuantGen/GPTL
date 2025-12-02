### GPTL using LD reference panels and GWAS results

The following example illustrate how GPTL software works when using an LD reference panel (internal or external) and GWAS results.

**1. Data Preparation**

We use a toy data set `wheatSumStats` of LD matrix `wheat_LD`, GWAS results `wheat_GWAS`, prior effects `wheat_PRIOR`, and individual validation/testing data `wheat_VLD.X, wheat_VLD.y, wheat_TST.X, wheat_TST.y`. These statistics were generated based on the wheat data set collected from CIMMYT's Global Wheat Program, including 599 wheat lines genotype (1279 variants) and phenotype (average grain yield).

```r
library(GPTL)
data(wheatSumStats)
library(Matrix)

wheat_LD[1:3,1:3]
#> 3 x 3 sparse Matrix of class "dgCMatrix"
#>            wPt.1171   c.312549   c.306034
#> wPt.1171 1.00000000 0.34408250 0.06038764
#> c.312549 0.34408250 1.00000000 0.03293636
#> c.306034 0.06038764 0.03293636 1.00000000

head(wheat_GWAS, 3)
#>                 id        beta        se   n allele_freq
#>  wPt.1171 wPt.1171 -0.03196163 0.2448336 167   0.4461078
#>  c.312549 c.312549 -0.13563632 0.1825260 167   0.3892216
#>  c.306034 c.306034 -0.41371141 0.1599184 167   0.3413174

str(wheat_PRIOR)
#>   Named num [1:1279] 0.01324 0.03401 -0.00575 0.00472 0.00775 ...
#>   - attr(*, "names")= chr [1:1279] "wPt.1171" "c.312549" "c.306034" "c.346957" ...
```

- *wheat_LD* is a sparse LD reference matrix in "dgCMatrix" class, with row and column names as variant ID.
- *wheat_GWAS* is a GWAS result data frame, where columns of **beta** (estimated effects), **n** (sample sizes), **allele_freq** (allele frequency), and variant ID as row names are required.
- *wheat_PRIOR* is a vector of prior effects estimated in the source population, with names as variant ID.

We use *getSS()* function to compute the sufficient statistics (**X'X** and **X'y**) based on the LD matrix, GWAS results (and prior effects, to algin the variants).

```R
SS=getSS(ld=wheat_LD,gwas=wheat_GWAS,B=wheat_PRIOR)
#>  There were 1279 variant in common between the LD reference panel, the GWAS, and the prior.
str(SS)
#>  List of 4
#>   $ XX:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>    .. ..@ i       : int [1:30477] 0 1 2 3 0 1 2 3 0 1 ...
#>    .. ..@ p       : int [1:1280] 0 4 8 12 16 20 24 28 32 38 ...
#>    .. ..@ Dim     : int [1:2] 1279 1279
#>    .. ..@ Dimnames:List of 2
#>    .. .. ..$ : chr [1:1279] "wPt.1171" "c.312549" "c.306034" "c.346957" ...
#>    .. .. ..$ : chr [1:1279] "wPt.1171" "c.312549" "c.306034" "c.346957" ...
#>    .. ..@ x       : num [1:30477] 82.53 27.85 4.75 4.88 27.85 ...
#>    .. ..@ factors : list()
#>   $ Xy: num [1:1279, 1] -2.638 -10.77 -31.066 -21.927 0.166 ...
#>    ..- attr(*, "dimnames")=List of 2
#>    .. ..$ : chr [1:1279] "wPt.1171" "c.312549" "c.306034" "c.346957" ...
#>    .. ..$ : NULL
#>   $ n : num 167
#>   $ B :'data.frame':	1279 obs. of  1 variable:
#>    ..$ B: num [1:1279] 0.01324 0.03401 -0.00575 0.00472 0.00775 ...
```

**2. PGS Estimation Using GPTL**

- #### Transfer Learning using Gradient Descent with Early Stopping (*TL-GDES*)

*GD()* function takes as input the above derived sufficient statistics and prior effects. The function returns regression coefficient values over the gradient descent cycles.

```R
fm_GDES=GD(XX=SS$XX, Xy=SS$Xy, b=SS$B, learningRate=1/100, nIter=100, returnPath=T)
dim(fm_GDES)
#> [1] 1279  100
```

We evaluate the prediction accuracy in the calibration set to select the optimal shrinkage parameter (nIter).

```R
Cor_GDES=cor(wheat_VLD.X %*% fm_GDES, wheat_VLD.y)
plot(Cor_GDES, xlab='iteration', ylab='Prediction Corr.', pch=20)
opt_nIter=which.max(Cor_GDES)
```

<p align="left">
    <img src="https://github.com/QuantGen/GPTL/blob/main/man/plots/E2_GDES.png" alt="Description" width="400">
</p>

We then re-estimate the PGS effects using the training set, with the optimal shrinkage parameter, and evaluate the final prediction accuracy in the testing set.

```R
fm_GDES_final=GD(XX=SS$XX, Xy=SS$Xy, b=SS$B, learningRate=1/100, nIter=opt_nIter, returnPath=F)
cor(wheat_TST.X %*% fm_GDES_final, wheat_TST.y)
#> [1] 0.606333
```

- #### Transfer Learning using penalized regressions (*TL-PR*)

*PR()* function takes as inputs the above derived sufficient statistics and prior effects, plus, potentially, values for $\lambda$ and $\alpha$ (if these are not provided, by default $\alpha$ is set equal to zero and the model is fitted over a grid of values of $\lambda$). The function returns estimates for a grid of values of $\lambda$, enabling users to select the optimal model based on cross-validation.

```R
fm_PR=PR(XX=SS$XX, Xy=SS$Xy, b=SS$B, alpha=0, nLambda=100, convThreshold=1e-4,
         maxIter=1000, returnPath=FALSE)
str(fm_PR)
#> List of 4
#>  $ B        : num [1:1279, 1:100] 0.01322 0.03397 -0.00585 0.00464 0.00775 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:1279] "wPt.1171" "c.312549" "c.306034" "c.346957" ...
#>   .. ..$ : chr [1:100] "lambda_617943.6592" "lambda_563047.2476" "lambda_513027.682" "lambda_467451.7167" ...
#>  $ lambda   : num [1:100] 617944 563047 513028 467452 425925 ...
#>  $ alpha    : num 0
#>  $ conv_iter: num [1:100] 3 3 3 3 3 3 3 3 3 3 ...
```

We evaluate the prediction accuracy in the calibration set to select the optimal shrinkage parameter (lambda).

```R
Cor_PR=cor(wheat_VLD.X %*% fm_PR$B, wheat_VLD.y)
plot(x=log(fm_PR$lambda), y=Cor_PR, xlab='log(lambda)', ylab='Prediction Corr.', pch=20)
opt_lambda=fm_PR$lambda[which.max(Cor_PR)]
```

<p align="left">
    <img src="https://github.com/QuantGen/GPTL/blob/main/man/plots/E2_PR.png" alt="Description" width="400">
</p>

We then re-estimate the PGS effects using the training set, with the optimal shrinkage parameter, and evaluate the final prediction accuracy in the testing set.

```R
fm_PR_final=PR(XX=SS$XX, Xy=SS$Xy, b=SS$B, alpha=0, lambda=opt_lambda, convThreshold=1e-4,
            maxIter=1000, returnPath=FALSE)
cor(wheat_TST.X %*% fm_PR_final$B, wheat_TST.y)
#> [1] 0.5645487
```

- #### Transfer Learning using Bayesian model with an informative finite mixture prior (*TL-BMM*)

*BMM()* function takes as inputs the above derived sufficient statistics prior, a matrix (B) whose columns contain the priors (one row per variant, one column per prior source of information), and parameters that control the algorithm. The function returns posterior means and posterior SDs for variant effects and other unknown parameters (including posterior ‘inclusion’ probabilities that link each variant effect to each of the components of the mixture).

When using sparse **X'X** as inputs for *BMM()*, variance parameters cannot be properly updated due to ignoring the LD between LD blocks. Therefore, in this circumstance, we suggest fixing the error and effects variances (by setting function parameters *fixVarE=FALSE, fixVarB=rep(FALSE,nComp)*). The error and effects variances will be computed based on *R2* parameters and we can profile *R2* over a grid of values (e.g., {0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5}). The optimal *R2* value can be selected in a calibration set.

```R
R2s=c(0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)
B_BMM=matrix(nrow=nrow(SS$XX), ncol=length(R2s))
for (i in 1:length(R2s)) {
    fm_BMM=BMM(XX=SS$XX, Xy=SS$Xy, B=cbind(SS$B,0), n=SS$n, R2=R2s[i],
           nIter=12000, burnIn=2000, thin=5, fixVarE=TRUE, fixVarB=rep(TRUE,2), verbose=FALSE)
    B_BMM[,i]=fm_BMM$b
}
```

We evaluate the prediction accuracy in the calibration set to select the optimal *R2* parameter.

```R
Cor_BMM=cor(wheat_VLD.X %*% B_BMM, wheat_VLD.y)
plot(x=R2s, y=Cor_BMM, xlab='R2', ylab='Prediction Corr.', pch=20)
opt_R2=R2s[which.max(Cor_BMM)]
```

<p align="left">
    <img src="https://github.com/QuantGen/GPTL/blob/main/man/plots/E2_BMM.png" alt="Description" width="400">
</p>

We then re-estimate the PGS effects using the training set, with the optimal *R2* parameter, and evaluate the final prediction accuracy in the testing set.

```R
fm_BMM_final=BMM(XX=SS$XX, Xy=SS$Xy, B=cbind(SS$B,0), n=SS$n, R2=opt_R2,
           nIter=12000, burnIn=2000, thin=5, fixVarE=TRUE, fixVarB=rep(TRUE,2), verbose=FALSE)
cor(wheat_TST.X %*% fm_BMM_final$b, wheat_TST.y)
#> [1] 0.5363532
```

[Back to Homepage](https://github.com/QuantGen/GPTL/blob/main/README.md)


