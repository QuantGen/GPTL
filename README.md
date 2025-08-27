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

### Estimating prior effects from the source data set

We estimated effects using a Bayesian shrinkage estimation method (a Bayesian model with a Gaussian prior centered at zero, model ‘BRR’ in the **BGLR** R-package).

```R
ETA=list(list(X=X_s, model="BRR"))
fm=BGLR(y=y_s, ETA = ETA, response_type = "gaussian", nIter = 12000, burnIn = 2000, verbose = FALSE)
prior=fm$ETA[[1]]$b
names(prior)=colnames(X_s)
```

Alternatively, if only sufficient statistics (**X'X** and **X'y**) for the source data set are provided, one can use *BLRCross()* function in the **BGLR** R-package.

### Transfer Learning using Gradient Descent with Early Stopping (*TL-GDES*)

GD() function takes as input the sufficient statistics derived from the target population and a vector of initial values (prior). The function returns regression coefficient values over the gradient descent cycles.

```R
fm_GDES=GD(XX=XX_t, Xy=Xy_t, b=prior, learning_rate=1/50, nIter=100, returnPath=T)
dim(fm_GDES)
#> [1] 1279  100
```

### Transfer Learning using penalized regressions (*TL-PR*)

PR.SS() function takes as inputs the sufficient statistics derived from the target population and a vector of initial values (prior), plus, potentially, values for $\lambda$ and $\alpha$ (if these are not provided, by default $\alpha$ is set equal to zero and the model is fitted over a grid of values of $\lambda$). The function returns estimates for a grid of values of $\lambda$ and $\alpha$, enabling users to select the optimal model based on cross-validation.

```R
fm_PR=PR(XX=XX_t, Xy=Xy_t, b=prior, alpha=0, nLambda=100, conv_threshold=1e-4,
         maxIter=1000, returnPath=FALSE)
str(fm_PR)
#> List of 4
#>  $ B        : num [1:1279, 1:100] 0.002611 0.003227 -0.00547 -0.00045 -0.000568 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:1279] "wPt.0538" "wPt.8463" "wPt.6348" "wPt.9992" ...
#>   .. ..$ : chr [1:100] "lambda_4614.6071" "lambda_4407.4852" "lambda_4209.5668" "lambda_4020.443" ...
#>  $ lambda   : num [1:100] 4615 4407 4210 4020 3840 ...
#>  $ alpha    : num 0
#>  $ conv_iter: num [1:100] 10 10 10 10 11 11 12 12 12 13 ...
```

### Transfer Learning using Bayesian model with an informative finite mixture prior (*TL-BMM*)

BMM.SS() function takes as inputs the sufficient statistics derived from the target population, a matrix (B) whose columns contain the priors (one row per variant, one column per prior source of information), and parameters that control the algorithm. The function returns posterior means and posterior SDs for variant effects and other unknown parameters (including posterior ‘inclusion’ probabilities that link each variant effect to each of the components of the mixture).

```R
fm_BMM=BMM(XX=XX_t, Xy=Xy_t, my=mean(y_t), vy=var(y_t), B=cbind(prior,0), n=nrow(X_t),
           nIter=12000, burnIn=2000, thin=5, verbose=FALSE)
str(fm_BMM)
#> List of 7
#>  $ b           : Named num [1:1279] -0.00131 0.01974 0.00557 0.00333 -0.00142 ...
#>   ..- attr(*, "names")= chr [1:1279] "wPt.0538" "wPt.8463" "wPt.6348" "wPt.9992" ...
#>  $ POST.PROB   : num [1:1279, 1:2] 0.429 0.537 0.482 0.482 0.49 ...
#>  $ postMeanVarB: num [1:2] 0.00159 0.00194
#>  $ postProb    : num [1:2] 0.499 0.501
#>  $ samplesVarB : num [1:12000, 1:2] 0.000703 0.000793 0.000812 0.0009 0.000918 ...
#>  $ samplesB    : num [1:12000, 1:1279] 0.00358 -0.07363 -0.01596 -0.02647 -0.06138 ...
#>  $ samplesVarE : num [1:12000] 0.799 0.667 0.664 0.672 0.721 ...
```

System Requirements
-------------------

Depends: R (>= 3.5.0)
