GPTL: Genomic Prediction Using Transfer Learning
================================================

GPTL is an R package that implements several methods for Genomic Prediction Using Transfer Learning, including
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

We use samples in cluster 1 as the source data set (where information is transferred) and samples in cluster 2 as the target data set (where the PGS will be used).  For the target data set, 20% of the samples are randomly selected and retained as a testing set. We compute the sufficient statistics (**X'X** and **X'y**) for each set.

```R
X_s=scale(wheat.X[CLUSTER$cluster == 1,], center=TRUE, scale=FALSE)
y_s=wheat.Y[CLUSTER$cluster == 1,1]
C_s=crossprod(X_s)
r_s=crossprod(X_s, y_s)

X_t=scale(wheat.X[CLUSTER$cluster == 2,], center=TRUE, scale=FALSE)
y_t=wheat.Y[CLUSTER$cluster == 2,1]
set.seed(217)
TST=sample(1:nrow(X_t), round(nrow(X_t)*0.2))
TRN=(1:nrow(X_t))[-TST]

C_t_TRN=crossprod(X_t[TRN,])
r_t_TRN=crossprod(X_t[TRN,], y_t[TRN])
C_t_TST=crossprod(X_t[TST,])
r_t_TST=crossprod(X_t[TST,], y_t[TST])
```

### Estimating prior effects from the source data set

We estimated effects using a Bayesian shrinkage estimation method (a Bayesian model with a Gaussian prior centered at zero, model ‘BRR’ in the **BGLR** R-package)

```R
idPriors=rep(1,ncol(C_s))
priors=list(list(model='BRR'))
fm<-BLRCross(n=nrow(X_s),vy=var(y_s),my=mean(y_s),XX=C_s,Xy=r_s,priors=priors,idPriors=idPriors,
             nIter=12000,burnIn=2000,verbose=TRUE)
prior=fm$ETA[[1]]$b
```

### Transfer Learning using Gradient Descent with Early Stopping (*TL-GDES*)

GD() function takes as input the sufficient statistics derived from the target population and a vector of initial values (prior). The function returns regression coefficient values over the GD cycles.

```R
fm_GDES=GD(C_t_TRN,r_t_TRN,b=prior,nIter=100,returnPath=T,learning_rate=1/50)
dim(fm_GDES)
#> [1] 1279  100
```

### Transfer Learning using penalized regressions (*TL-PR*)

PR() function takes as inputs the sufficient statistics plus, potentially, values for $\lambda$ and $\alpha$ (if these are not provided, by default $\alpha$ is set equal to zero and the model is fitted over a grid of values of $\lambda$). The function returns estimates for a grid of values of $\lambda$ and $\alpha$.

```R
fm_PR=PR(C_t_TRN, r_t_TRN, b0=prior, alpha=0, nLambda=100, conv_threshold=1e-4, maxIter=1000, returnPath=FALSE)
str(fm_PR)
#> List of 4
#>  $ B        : num [1:1279, 1:100] 0.01016 0.02012 0.0016 0.00335 0.00539 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:1279] "wPt.0538" "wPt.8463" "wPt.6348" "wPt.9992" ...
#>   .. ..$ : chr [1:100] "lambda_3096.509" "lambda_2957.5254" "lambda_2824.7175" "lambda_2697.811" ...
#>  $ lambda   : num [1:100] 3097 2958 2825 2698 2577 ...
#>  $ alpha    : num 0
#>  $ conv_iter: num [1:100] 6 6 6 7 7 7 7 7 7 8 ...
```

### Transfer Learning using Bayesian model with an informative finite mixture prior (*TL-BMM*)

PR() function takes as inputs the sufficient statistics plus, potentially, values for $\lambda$ and $\alpha$ (if these are not provided, by default $\alpha$ is set equal to zero and the model is fitted over a grid of values of $\lambda$). The function returns estimates for a grid of values of $\lambda$ and $\alpha$.

BMM() function takes as inputs the sufficient statistics from the target population, a matrix (B) whose columns contain the prior means (one row per SNP, one column per prior source of information), and parameters that control the algorithm. The function returns posterior means (and posterior SDs) for SNP effects and other unknown parameters (including posterior ‘inclusion’ probabilities that link each SNP effect to each of the components of the mixture).

```R
fm_BMM=BMM_new(C=C_t_TRN, rhs=r_t_TRN, my=mean(y_t[TRN]), vy=var(y_t[TRN]), nIter=12000, burnIn=2000, thin=5, verbose=FALSE, B0=cbind(prior,0), n=nrow(X_t[TRN,]))
str(fm_BMM)
#> List of 7
#>  $ b           : num [1:1279] -5.61e-03 1.70e-02 2.51e-02 -5.04e-05 3.41e-03 ...
#>  $ POST.PROB   : num [1:1279, 1:2] 0.469 0.536 0.458 0.477 0.511 ...
#>  $ postMeanVarB: num [1:2] 0.00172 0.00221
#>  $ postProb    : num [1:2] 0.499 0.501
#>  $ samplesVarB : num [1:12000, 1:2] 0.000603 0.000589 0.000536 0.000591 0.000619 ...
#>  $ samplesB    : num [1:12000, 1:1279] 0.03207 0.00817 0.03011 -0.00755 -0.01431 ...
#>  $ samplesVarE : num [1:12000] 0.639 0.552 0.722 0.769 0.575 ...
```
