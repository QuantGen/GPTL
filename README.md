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

Examples
--------

### Loading the package

```R
library(GPTL)
```

### Loading example dataset

```R
library(BGLR)
data(wheat)
dim(wheat.X)
#> [1] 599 1279
dim(wheat.Y)
#> [1] 599 4
X=scale(wheat.X, center=TRUE, scale=TRUE)
CLUSTER=kmeans(X,2)
table(CLUSTER$cluster)
#>   1   2 
#> 346 253 
```

We use samples in cluster 1 as the source data set (where information is transferred) and samples in cluster 2 as the target data set (where the PGS will be used). We compute the sufficient statistics (**X'X** and **X'y**) for each set.

```R
X_s=scale(wheat.X[CLUSTER$cluster == 1,], center=TRUE, scale=FALSE)
y_s=wheat.Y[CLUSTER$cluster == 1,1]
C_s=crossprod(X_s)
r_s=crossprod(X_s, y_s)

X_t=scale(wheat.X[CLUSTER$cluster == 2,], center=TRUE, scale=FALSE)
y_t=wheat.Y[CLUSTER$cluster == 2,1]
C_t=crossprod(X_t)
r_t=crossprod(X_t, y_t)
```

### Estimating prior effects from the source data set

We estimated effects using a Bayesian shrinkage estimation method (a Bayesian model with a Gaussian prior centered at zero, model ‘BRR’ in the **BGLR** R-package)

```R
idPriors=rep(1,ncol(C_s))
priors=list(list(model='BRR'))
fm<-BLRCross(n=nrow(X_s),vy=var(y_s),my=mean(y_s),XX=C_s,Xy=r_s,priors=priors,
             idPriors=idPriors, nIter=12000,burnIn=2000,verbose=TRUE)
prior=fm$ETA[[1]]$b
```

### Transfer Learning using Gradient Descent with Early Stopping (*TL-GDES*)

GD() function takes as input the sufficient statistics derived from the target population and a vector of initial values (prior). The function returns regression coefficient values over the GD cycles.

```R
fm_GDES=GD(C_t,r_t,b=prior,nIter=100,returnPath=T,learning_rate=1/50)
dim(fm_GDES)
#> [1] 1279  100
```

### Transfer Learning using penalized regressions (*TL-PR*)

PR() function takes as inputs the sufficient statistics plus, potentially, values for $\lambda$ and $\alpha$ (if these are not provided, by default $\alpha$ is set equal to zero and the model is fitted over a grid of values of $\lambda$). The function returns estimates for a grid of values of $\lambda$ and $\alpha$, enabling users to select the optimal model based on cross-validation.

```R
fm_PR=PR(C_t, r_t, b0=prior, alpha=0, nLambda=100, conv_threshold=1e-4,
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

### Transfer Learning using Bayesian model with an informative finite mixture prior (*TL-BMM*)

BMM() function takes as inputs the sufficient statistics from the target population, a matrix (B) whose columns contain the prior means (one row per SNP, one column per prior source of information), and parameters that control the algorithm. The function returns posterior means (and posterior SDs) for SNP effects and other unknown parameters (including posterior ‘inclusion’ probabilities that link each SNP effect to each of the components of the mixture).

```R
fm_BMM=BMM_new(C=C_t, rhs=r_t, my=mean(y_t), vy=var(y_t), nIter=12000, burnIn=2000, thin=5,
               verbose=FALSE, B0=cbind(prior,0), n=nrow(X_t))
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
