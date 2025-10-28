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

- *LD* is a sparse LD reference matrix in "dgCMatrix" class, with row and column names as variant ID.
- *gwas* is a GWAS result data frame, where columns of **beta** (estimated effects), **se** (standard errors), **n** (sample sizes), **allele_freq** (allele frequency), and variant ID as row names are required.
- *prior* is a vector of prior effects estimated in the source population, with names as variant ID.

We use *getSS()* function to compute the sufficient statistics (**X'X** and **X'y**) based on the LD matrix, GWAS results (and prior effects, to algin the variants).

```R
SS=getSS(ld=LD,gwas=gwas,B=prior)
#>  There were 1719 variant in common between the LD reference panel, the GWAS and the prior.
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

The PGS estimation processes are similar to [Example_Individual_Data](https://github.com/QuantGen/GPTL/blob/main/man/Example_Individual_Data.md) (except for *BMM()*). Therefore, here we only illustrate how GPTL software works, without showing the details of calibration-testing processes.

- #### Transfer Learning using Gradient Descent with Early Stopping (*TL-GDES*)

*GD()* function takes as input the above derived sufficient statistics and prior effects. The function returns regression coefficient values over the gradient descent cycles.

```R
fm_GDES=GD(XX=SS$XX, Xy=SS$Xy, b=SS$B, learningRate=1/50, nIter=100, returnPath=T)
dim(fm_GDES)
#> [1] 1719  100
```

- #### Transfer Learning using penalized regressions (*TL-PR*)

*PR()* function takes as inputs the above derived sufficient statistics and prior effects, plus, potentially, values for $\lambda$ and $\alpha$ (if these are not provided, by default $\alpha$ is set equal to zero and the model is fitted over a grid of values of $\lambda$). The function returns estimates for a grid of values of $\lambda$, enabling users to select the optimal model based on cross-validation.

```R
fm_PR=PR(XX=SS$XX, Xy=SS$Xy, b=SS$B, alpha=0, nLambda=100, conv_threshold=1e-4,
         maxIter=1000, returnPath=FALSE)
str(fm_PR)
#> List of 4
#>  $ B        : num [1:1719, 1:100] 7.56e-06 5.86e-04 8.27e-06 2.66e-06 6.96e-06 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:1719] "JHU_1.737262" "rs114339855" "JHU_1.761190" "JHU_1.761763" ...
#>   .. ..$ : chr [1:100] "lambda_86445040.9413" "lambda_78765501.7566" "lambda_71768191.6674" "lambda_65392503.3211" ...
#>  $ lambda   : num [1:100] 86445041 78765502 71768192 65392503 59583214 ...
#>  $ alpha    : num 0
#>  $ conv_iter: num [1:100] 2 2 2 2 3 3 3 3 3 3 ...
```

- #### Transfer Learning using Bayesian model with an informative finite mixture prior (*TL-BMM*)

*BMM()* function takes as inputs the above derived sufficient statistics prior, a matrix (B) whose columns contain the priors (one row per variant, one column per prior source of information), and parameters that control the algorithm. The function returns posterior means and posterior SDs for variant effects and other unknown parameters (including posterior ‘inclusion’ probabilities that link each variant effect to each of the components of the mixture).

When using sparse **X'X** as inputs for *BMM()*, variance parameters cannot be properly updated due to ignoring the LD between LD blocks. Therefore, in this circumstance, we suggest fixing the error and effects variances (by setting function parameters *fixVarE=FALSE, fixVarB=rep(FALSE,nComp)*). The error and effects variances will be computed based on *R2* parameters and we can profile *R2* over a grid of values (e.g., {0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5}). The optimal *R2* value can be selected in a calibration set (not shown here).

```R
fm_BMM=BMM(XX=SS$XX, Xy=SS$Xy, my=my, vy=vy, B=cbind(SS$B,0), n=N,
           nIter=12000, burnIn=2000, thin=5, fixVarE=TRUE, fixVarB=rep(TRUE,2), verbose=FALSE)
str(fm_BMM)
#> List of 7
#>  $ b           : Named num [1:1719] 0.0454 -0.0768 0.0818 0.021 0.041 ...
#>   ..- attr(*, "names")= chr [1:1719] "JHU_1.737262" "rs114339855" "JHU_1.761190" "JHU_1.761763" ...
#>  $ POST.PROB   : num [1:1719, 1:2] 0.481 0.513 0.498 0.492 0.488 ...
#>  $ postMeanVarB: num [1:2] 0.0177 0.0177
#>  $ postProb    : num [1:2] 0.501 0.499
#>  $ samplesVarB : num [1:12000, 1:2] 0.0177 0.0177 0.0177 0.0177 0.0177 ...
#>  $ samplesB    : num [1:12000, 1:1719] -0.00363 0.00521 0.03518 0.02083 0.15977 ...
#>  $ samplesVarE : num [1:12000] 40 40 40 40 40 ...
```

[Back to Homepage](https://github.com/QuantGen/GPTL/blob/main/README.md)


