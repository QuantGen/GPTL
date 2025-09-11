### GPTL using LD reference panels and GWAS results

The following example illustrate how GPTL software works when using an LD reference panel (internal or external) and GWAS results.

**1. Data Preparation**

We use a toy data set of LD matrix, GWAS results, and prior effects, including 1947 variants (first 10 LD blocks in chromosome 1 of All of Us African American LD reference panel). The whole LD reference panels and GWAS results are at [Link](https://doi.org/10.5281/zenodo.16923734) and [Link](https://doi.org/10.5281/zenodo.17087604).

```R
library(GPTL)
library(Matrix)
data(toyData)

LD[1:3,1:3]
#> 3 x 3 sparse Matrix of class "dgCMatrix"
#>              JHU_1.737262  rs114339855 JHU_1.761190
#> JHU_1.737262  1.000000000 0.0010638559 0.0015435403
#> rs114339855   0.001063856 1.0000000000 0.0003776402
#> JHU_1.761190  0.001543540 0.0003776402 1.0000000000

head(gwas, 3)
#>                        id chr a1 a0        beta        se     n allele_freq
#> JHU_1.737262 JHU_1.737262   1  A  G  0.06759564 0.1091687 39302  0.06583625
#> rs114339855   rs114339855   1  G  T -0.18147248 0.1860653 39302  0.02037881
#> JHU_1.761190 JHU_1.761190   1  T  C  0.16778221 0.1603125 39302  0.02787886

str(prior)
#>  Named num [1:1947] 0 0.000592 0 0 0 ...
#>  - attr(*, "names")= chr [1:1947] "JHU_1.737262" "rs114339855" "JHU_1.761190" "JHU_1.761763" ...
```

Here *LD* is a sparse LD reference matrix in "dgCMatrix" class, with row and column names as variant ID. *gwas* is a GWAS result data frame, where columns of **beta** (estimated effects), **se** (standard errors), **n** (sample sizes), **allele_freq** (allele frequency), and variant ID as row names are required. *prior* is a vector of prior effects estimated in the source population, with names as variant ID.

We use *getSS()* function to compute the sufficient statistics (**X'X** and **X'y**) based on the LD matrix, GWAS results (and prior effects, to algin the variants).

```R
SS=getSS(ld=LD,gwas=gwas,B=prior)
#>  There were 1719 variant in common between the LD reference panel, the GWAS and, the prior.
str(SS)
#> List of 4
#>  $ XX:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   .. ..@ i       : int [1:604305] 0 1 2 3 4 0 1 2 3 4 ...
#>   .. ..@ p       : int [1:1720] 0 5 10 15 20 25 259 493 727 961 ...
#>   .. ..@ Dim     : int [1:2] 1719 1719
#>   .. ..@ Dimnames:List of 2
#>   .. .. ..$ : chr [1:1719] "JHU_1.737262" "rs114339855" "JHU_1.761190" "JHU_1.761763" ...
#>   .. .. ..$ : chr [1:1719] "JHU_1.737262" "rs114339855" "JHU_1.761190" "JHU_1.761763" ...
#>   .. ..@ x       : num [1:604305] 4834.29 2.93 4.95 11.57 13.34 ...
#>   .. ..@ factors : list()
#>  $ Xy: Named num [1:1719] 327 -285 357 115 301 ...
#>   ..- attr(*, "names")= chr [1:1719] "JHU_1.737262" "rs114339855" "JHU_1.761190" "JHU_1.761763" ...
#>  $ n : num 39302
#>  $ B :'data.frame':	1719 obs. of  1 variable:
#>   ..$ B: num [1:1719] 0 0.000592 0 0 0 ...
```

**2. PGS Estimation Using GPTL**




- #### Transfer Learning using Gradient Descent with Early Stopping (*TL-GDES*)

GD() function takes as input the sufficient statistics derived from the target population and a vector of initial values (prior). The function returns regression coefficient values over the gradient descent cycles.

```R
fm_GDES=GD(XX=XXt_train, Xy=Xyt_train, b=prior, learning_rate=1/50, nIter=100, returnPath=T)
dim(fm_GDES)
#> [1] 1279  100
```

We evaluate the prediction accuracy in the calibration set to select the optimal shrinkage parameter (nIter).

```R
Cor_GDES=getCor(XXt_cali, Xyt_cali, yyt_cali, fm_GDES)
plot(Cor_GDES, xlab='iteration', ylab='Prediction Corr.', pch=20)
opt_nIter=which.max(Cor_GDES)
```

<p align="left">
    <img src="https://github.com/QuantGen/GPTL/blob/main/man/plots/GDES_plot.png" alt="Description" width="400">
</p>

We then re-estimate the PGS effects using both the training and calibration sets, with the optimal shrinkage parameter, and evaluate the final prediction accuracy in the testing set.

```R
fm_GDES_final=GD(XX=XXt_train+XXt_cali, Xy=Xyt_train+Xyt_cali, b=prior, learning_rate=1/50, nIter=opt_nIter, returnPath=F)
getCor(XXt_test, Xyt_test, yyt_test, fm_GDES_final)
#> [1] 0.4334093
```

- #### Transfer Learning using penalized regressions (*TL-PR*)

PR() function takes as inputs the sufficient statistics derived from the target population and a vector of initial values (prior), plus, potentially, values for $\lambda$ and $\alpha$ (if these are not provided, by default $\alpha$ is set equal to zero and the model is fitted over a grid of values of $\lambda$). The function returns estimates for a grid of values of $\lambda$, enabling users to select the optimal model based on cross-validation.

```R
fm_PR=PR(XX=XXt_train, Xy=Xyt_train, b=prior, alpha=0, nLambda=100, conv_threshold=1e-4,
         maxIter=1000, returnPath=FALSE)
str(fm_PR)
#> List of 4
#>  $ B        : num [1:1279, 1:100] 0.01565 0.015802 -0.005457 0.002243 0.000514 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:1279] "wPt.0538" "wPt.8463" "wPt.6348" "wPt.9992" ...
#>   .. ..$ : chr [1:100] "lambda_158389.9608" "lambda_144319.0332" "lambda_131498.1281" "lambda_119816.1968" ...
#>  $ lambda   : num [1:100] 158390 144319 131498 119816 109172 ...
#>  $ alpha    : num 0
#>  $ conv_iter: num [1:100] 3 3 3 3 3 3 3 3 3 3 ...
```

We evaluate the prediction accuracy in the calibration set to select the optimal shrinkage parameter (lambda).

```R
Cor_PR=getCor(XXt_cali, Xyt_cali, yyt_cali, fm_PR$B)
plot(x=log(fm_PR$lambda), y=Cor_PR, xlab='log(lambda)', ylab='Prediction Corr.', pch=20)
opt_lambda=fm_PR$lambda[which.max(Cor_PR)]
```

<p align="left">
    <img src="https://github.com/QuantGen/GPTL/blob/main/man/plots/PR_plot.png" alt="Description" width="400">
</p>

We then re-estimate the PGS effects using both the training and calibration sets, with the optimal shrinkage parameter, and evaluate the final prediction accuracy in the testing set.

```R
fm_PR_final=PR(XX=XXt_train+XXt_cali, Xy=Xyt_train+Xyt_cali, b=prior, alpha=0, lambda=opt_lambda, conv_threshold=1e-4,
            maxIter=1000, returnPath=FALSE)
getCor(XXt_test, Xyt_test, yyt_test, fm_PR_final$B)
#> [1] 0.6141271
```

- #### Transfer Learning using Bayesian model with an informative finite mixture prior (*TL-BMM*)

BMM() function takes as inputs the sufficient statistics derived from the target population, a matrix (B) whose columns contain the priors (one row per variant, one column per prior source of information), and parameters that control the algorithm. The function returns posterior means and posterior SDs for variant effects and other unknown parameters (including posterior ‘inclusion’ probabilities that link each variant effect to each of the components of the mixture).

BMM() only requires a single run of the algorithm because regularization parameters and variant effects are jointly inferred from the posterior distribution. Thus, this method does not require calibrating regularization parameters. We estimate the PGS effects using both the training and calibration sets, and evaluate the final prediction accuracy in the testing set.

```R
fm_BMM=BMM(XX=XXt_train+XXt_cali, Xy=Xyt_train+Xyt_cali, my=mean(c(yt_train,yt_cali)), vy=var(c(yt_train,yt_cali)n), B=cbind(prior,0), n=nrow(Xt_train)+nrow(Xt_cali),
           nIter=12000, burnIn=2000, thin=5, verbose=FALSE)
str(fm_BMM)
#> List of 7
#>  $ b           : Named num [1:1279] -0.020117 0.017591 0.019884 0.003021 -0.000865 ...
#>   ..- attr(*, "names")= chr [1:1279] "wPt.0538" "wPt.8463" "wPt.6348" "wPt.9992" ...
#>  $ POST.PROB   : num [1:1279, 1:2] 0.436 0.522 0.479 0.505 0.504 ...
#>  $ postMeanVarB: num [1:2] 0.0014 0.00235
#>  $ postProb    : num [1:2] 0.5 0.5
#>  $ samplesVarB : num [1:12000, 1:2] 0.000726 0.000718 0.000772 0.000746 0.000772 ...
#>  $ samplesB    : num [1:12000, 1:1279] 0.0204 -0.0154 -0.0104 -0.0226 0.0205 ...
#>  $ samplesVarE : num [1:12000] 0.593 0.537 0.524 0.585 0.458 ...
getCor(XXt_test, Xyt_test, yyt_test, fm_BMM$b)
#> [1] 0.635143
```

[Back to Homepage](https://github.com/QuantGen/GPTL/blob/main/README.md)












###### Not Using


Alternatively, one can provide a LD reference panel (i.e., a matrix of correlations between the variants) and GWAS results (including variants allele frequencies, estimated effects, and SEs) if the sufficient statistics (**X'X** and **X'y**) or individual genotype and phenotype (**X** and **y**) for the target population are not available.

```R
str(ld)
#> num [1:1279, 1:1279] 1 0.12 0.102 0.256 0.241 ...
#> - attr(*, "dimnames")=List of 2
#>  ..$ : chr [1:1279] "wPt.0538" "wPt.8463" "wPt.6348" "wPt.9992" ...
#>  ..$ : chr [1:1279] "wPt.0538" "wPt.8463" "wPt.6348" "wPt.9992" ...

head(gwas, 3)
#>                beta       se         n          allele_freq
#> wPt.0538       0.00235    0.01269    15000      0.55476
#> wPt.8463       -0.01228   0.01206    15000      0.44582
#> wPt.6348       0.00989    0.01162    15000      0.63647
```

In this study we generated LD reference panels for African American and Hispanic ancestry populations using the All of Us (CDRv7, Controlled Tier) data, available at [Link](https://doi.org/10.5281/zenodo.16923734).

getSS() function takes as input the LD reference panel and GWAS results. The function returns the sufficient statistics (**X'X** and **X'y**).

```R
SS=getSS(ld, gwas)
str(SS)
#> List of 3
#>  $ XX:Formal class 'dgeMatrix' [package "Matrix"] with 4 slots
#>   .. ..@ Dim     : int [1:2] 1279 1279
#>   .. ..@ Dimnames:List of 2
#>   .. .. ..$ : chr [1:1279] "wPt.0538" "wPt.8463" "wPt.6348" "wPt.9992" ...
#>   .. .. ..$ : chr [1:1279] "wPt.0538" "wPt.8463" "wPt.6348" "wPt.9992" ...
#>   .. ..@ x       : num [1:1635841] 446.16 3.73 25.69 -69.47 -3.84 ...
#>   .. ..@ factors : list()
#>  $ Xy: Named num [1:1279] -14.1 -20.8 -54.7 38.9 -56.8 ...
#>   ..- attr(*, "names")= chr [1:1279] "wPt.0538" "wPt.8463" "wPt.6348" "wPt.9992" ...
#>  $ n : num 15000
```




#>  $ samplesVarE : num [1:12000] 0.799 0.667 0.664 0.672 0.721 ...
```

