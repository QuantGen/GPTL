### GPTL using LD reference panels and GWAS results

The following example illustrate how GPTL software works when using an LD reference panel (internal or external) and GWAS results.

**1. Data Preparation**

We use a toy data set `data(wheatSumStats)` of LD matrix, GWAS results, prior effects, and individual validation/testing data. These statistics were generated based on the wheat data set collected from CIMMYT's Global Wheat Program, including 599 wheat lines genotype (1279 variants) and phenotype (average grain yield).

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
#>  wPt.0538 wPt.0538 -0.44041492 0.1480735 167   0.2604790
#>  wPt.8463 wPt.8463  0.02897992 0.2338454 167   0.4401198
#>  wPt.6348 wPt.6348  0.11440060 0.1528133 167   0.2814371

str(wheat_PRIOR)
#>  Named num [1:1279] 0.01508 0.01603 -0.00451 0.00243 0.00337 ...
#>  - attr(*, "names")= chr [1:1279] "wPt.0538" "wPt.8463" "wPt.6348" "wPt.9992" ...
```

- *wheat_LD* is a sparse LD reference matrix in "dgCMatrix" class, with row and column names as variant ID.
- *wheat_GWAS* is a GWAS result data frame, where columns of **beta** (estimated effects), **se** (standard errors), **n** (sample sizes), **allele_freq** (allele frequency), and variant ID as row names are required.
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
#>    ..$ B: num [1:1279] 0.01451 0.03475 -0.00429 0.00121 0.009 ...
```

**2. PGS Estimation Using GPTL**

The PGS estimation processes are similar to [Example_Individual_Data](https://github.com/QuantGen/GPTL/blob/main/man/Example_Individual_Data.md) (except for *BMM()*). Therefore, here we only illustrate how GPTL software works, without showing the details of calibration-testing processes.

- #### Transfer Learning using Gradient Descent with Early Stopping (*TL-GDES*)

*GD()* function takes as input the above derived sufficient statistics and prior effects. The function returns regression coefficient values over the gradient descent cycles.

```R
fm_GDES=GD(XX=SS$XX, Xy=SS$Xy, b=SS$B, learningRate=1/50, nIter=100, returnPath=T)
dim(fm_GDES)
#> [1] 1279  100
```

- #### Transfer Learning using penalized regressions (*TL-PR*)

*PR()* function takes as inputs the above derived sufficient statistics and prior effects, plus, potentially, values for $\lambda$ and $\alpha$ (if these are not provided, by default $\alpha$ is set equal to zero and the model is fitted over a grid of values of $\lambda$). The function returns estimates for a grid of values of $\lambda$, enabling users to select the optimal model based on cross-validation.

```R
fm_PR=PR(XX=SS$XX, Xy=SS$Xy, b=SS$B, alpha=0, nLambda=100, convThreshold=1e-4,
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


