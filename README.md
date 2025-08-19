GPTL: Genomic Prediction Using Transfer Learning
================================================

GPTL is an R package that implements several methods for Genomic Prediction using Transfer Learning, including
 -  *TL-GDES*: A Gradient Descent Algorithm with Early Stopping.
 -  *TL-PR*: A Penalized regression with srhinkage towards a prior mean (e.g., estimates derived from another data set).
 -  *TL-BMM*: A Bayesian model with a finite mixture prior that allows TL from multiple prior sources of information.

Installation
------------

To install from Github:

```R
# install.packages("remotes")
remotes::install_github("QuantGen/GPTL")
```

This should take less than 30 seconds for a typical computer when including installation of all non-base dependencies.

This package should be compatible with Windows, Mac, and Linux operating systems and has been tested on Windows 7 & 10, macOS Sequoia & Sonoma, and Linux CentOS 7.

Examples
--------

### Loading the package

```R
library(GPTL)
```

### Loading example dataset

Here we use the [wheat](https://doi.org/10.1104/pp.105.063438) data set collected from CIMMYT's Global Wheat Program. We use the 599 wheat lines genotype (1279 variants) and phenotype (average grain yield) data to illustrate the methods in GPTL.

```R
library(BGLR)
data(wheat)
y=wheat.Y[,1]
X=scale(wheat.X, center=TRUE, scale=TRUE)

CLUSTER=kmeans(X,2)
table(CLUSTER$cluster)
#>   1   2 
#> 346 253 
```

We use samples in cluster 1 as the source data set (where information is transferred) and samples in cluster 2 as the target data set (where the PGS will be used). 

```R
X_s=scale(wheat.X[CLUSTER$cluster == 1,], center=TRUE, scale=FALSE)
y_s=y[CLUSTER$cluster == 1]

X_t=scale(wheat.X[CLUSTER$cluster == 2,], center=TRUE, scale=FALSE)
y_t=y[CLUSTER$cluster == 2]
```

We compute the sufficient statistics (**X'X** and **X'y**) for the target data set.

```R
XX_t=crossprod(X_t)
Xy_t=crossprod(X_t, y_t)
```

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

### Estimating prior effects from the source data set

We estimated effects using a Bayesian shrinkage estimation method (a Bayesian model with a Gaussian prior centered at zero, model ‘BRR’ in the **BGLR** R-package).

```R
ETA=list(list(X=X_s, model="BRR"))
fm=BGLR(y=y_s, response_type = "gaussian", ETA = ETA, nIter = 12000,
       burnIn = 2000, verbose = FALSE)
prior=fm$ETA[[1]]$b
names(prior)=colnames(X_s)
```

Alternatively, if only sufficient statistics (**X'X** and **X'y**) for the source data set are provided, one can use *BLRCross()* function in the **BGLR** R-package.

### Transfer Learning using Gradient Descent with Early Stopping (*TL-GDES*)

GD.SS() function takes as input the sufficient statistics derived from the target population and a vector of initial values (prior). The function returns regression coefficient values over the gradient descent cycles.

```R
fm_GDES=GD.SS(XX=XX_t, Xy=Xy_t, b=prior, nIter=100, returnPath=T, learning_rate=1/50)
dim(fm_GDES)
#> [1] 1279  100
```

Alternatively, GD() function takes as input a LD reference panel and GWAS results from the target population, and a vector of initial values (prior).

```R
fm_GDES2=GD(ld=ld, gwas=gwas, b=prior, nIter=100, returnPath=T, learning_rate=1/50)
```

### Transfer Learning using penalized regressions (*TL-PR*)

PR.SS() function takes as inputs the sufficient statistics derived from the target population and a vector of initial values (prior), plus, potentially, values for $\lambda$ and $\alpha$ (if these are not provided, by default $\alpha$ is set equal to zero and the model is fitted over a grid of values of $\lambda$). The function returns estimates for a grid of values of $\lambda$ and $\alpha$, enabling users to select the optimal model based on cross-validation.

```R
fm_PR=PR.SS(XX=XX_t, Xy=Xy_t, b=prior, alpha=0, nLambda=100, conv_threshold=1e-4,
         maxIter=1000, returnPath=FALSE)
str(fm_PR)
#> List of 4
#>  $ B        : num [1:1279, 1:100] 0.00634 0.01794 -0.00267 0.00389 0.00399 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:1279] "wPt.0538" "wPt.8463" "wPt.6348" "wPt.9992" ...
#>   .. ..$ : chr [1:100] "lambda_3886.4305" "lambda_3711.9921" "lambda_3545.3049" "lambda_3386.0244" ...
#>  $ lambda   : num [1:100] 3886 3712 3545 3386 3234 ...
#>  $ alpha    : num 0
#>  $ conv_iter: num [1:100] 6 6 6 6 7 7 7 7 7 7 ...
```

Alternatively, PR() function takes as input a LD reference panel and GWAS results from the target population, a vector of initial values (prior), and regularization parameters.

```R
fm_PR2=PR(ld=ld, gwas=gwas, b=prior, alpha=0, nLambda=100, conv_threshold=1e-4,
         maxIter=1000, returnPath=FALSE)
```

### Transfer Learning using Bayesian model with an informative finite mixture prior (*TL-BMM*)

BMM.SS() function takes as inputs the sufficient statistics derived from the target population, a matrix (B) whose columns contain the priors (one row per variant, one column per prior source of information), and parameters that control the algorithm. The function returns posterior means and posterior SDs for variant effects and other unknown parameters (including posterior ‘inclusion’ probabilities that link each variant effect to each of the components of the mixture).

```R
fm_BMM=BMM.SS(XX=XX_t, Xy=Xy_t, my=mean(y_t), vy=var(y_t), nIter=12000, burnIn=2000, thin=5,
               verbose=FALSE, B=cbind(prior,0), n=nrow(X_t))
str(fm_BMM)
#> List of 7
#>  $ b           : num [1:1279] -0.01643 0.01787 0.02097 0.00445 -0.00112 ...
#>  $ POST.PROB   : num [1:1279, 1:2] 0.471 0.499 0.493 0.493 0.508 ...
#>  $ postMeanVarB: num [1:2] 0.00261 0.00207
#>  $ postProb    : num [1:2] 0.5 0.5
#>  $ samplesVarB : num [1:12000, 1:2] 0.000655 0.000675 0.000647 0.000699 0.000708 ...
#>  $ samplesB    : num [1:12000, 1:1279] 0.002978 0.038193 -0.015793 -0.003802 0.000238 ...
#>  $ samplesVarE : num [1:12000] 0.643 0.582 0.549 0.629 0.604 ...
```

Alternatively, BMM() function takes as input a LD reference panel and GWAS results from the target population, a matrix of initial values (prior), and parameters that control the algorithm.

```R
fm_BMM2=BMM(ld=ld, gwas=gwas, B=cbind(prior,0), my=mean(y_t), vy=var(y_t), nIter=12000, burnIn=2000, thin=5,
               verbose=FALSE, n=nrow(X_t))
```

System Requirements
-------------------

Depends: R (>= 3.5.0)
