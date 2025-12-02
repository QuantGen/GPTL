### GPTL using individual genetype and phenotype data

The following example illustrate how GPTL software works when one has access to individual genetype and phenotype data. Here we use the [simulated human genotype and phenotype data sets](https://github.com/QuantGen/GPTL/blob/main/man/Data_Simulation.md) (N: Source-10,000, Target-4,000), including 2450 variants.

**1. Data Preparation**

```R
library(GPTL)
data(Ind_SimData)
library(BGLR)
```

We compute the sufficient statistics (X'X and X'y) for the calibrating and testing sets in the target population.

```R
XXt_trn=crossprod(scale(Xt_trn,center=TRUE,scale=FALSE));Xyt_trn=crossprod(scale(Xt_trn,center=TRUE,scale=FALSE), yt_trn)
XXt_cal=crossprod(scale(Xt_cal,center=TRUE,scale=FALSE));Xyt_cal=crossprod(scale(Xt_cal,center=TRUE,scale=FALSE), yt_cal);yyt_cal=crossprod(yt_cal)
XXt_tst=crossprod(scale(Xt_tst,center=TRUE,scale=FALSE));Xyt_tst=crossprod(scale(Xt_tst,center=TRUE,scale=FALSE), yt_tst);yyt_tst=crossprod(yt_tst)

```

**2. Prior Estimation**

We estimated prior effects from the source data set using a Bayesian shrinkage estimation method (a Bayesian model with a Gaussian prior centered at zero, model ‘BRR’ in the **BGLR** R-package). Alternatively, if only sufficient statistics (**X'X** and **X'y**) for the source data set are provided, one can use *BLRCross()* function in the **BGLR** R-package.

```R
LP=list(snps=list(X=scale(Xs,center=TRUE,scale=FALSE),model='BayesC'))
fm=BLRXy(y=ys,ETA=LP,verbose=FALSE,nIter=6000,burnIn=1000,verbose=FALSE)
prior=fm$ETA[[1]]$b
names(prior)=colnames(Xs)
```

**3. PGS Estimation Using GPTL**

- #### Transfer Learning using Gradient Descent with Early Stopping (*TL-GDES*)

*GD()* function takes as input the sufficient statistics derived from the target population and a vector of initial values (prior). The function returns regression coefficient values over the gradient descent cycles.

```R
fm_GDES=GD(XX=XXt_trn, Xy=Xyt_trn, b=prior, learningRate=1/200, nIter=100, returnPath=T)
dim(fm_GDES)
#> [1] 2450  100
```

We evaluate the prediction accuracy in the calibrating set to select the optimal number of gradient descent cycles (nIter).

```R
Cor_GDES=getCor(XXt_cal, Xyt_cal, yyt_cal, fm_GDES)
plot(Cor_GDES, xlab='iteration', ylab='Prediction Corr.', pch=20, type='o')
opt_nIter=which.max(Cor_GDES)
```

<p align="left">
    <img src="https://github.com/QuantGen/GPTL/blob/main/man/plots/E1_GDES.png" alt="Description" width="400">
</p>

We then re-estimate the PGS effects using both the training and calibrating sets, with the optimal shrinkage parameter, and evaluate the final prediction accuracy in the testing set.

```R
fm_GDES_final=GD(XX=XXt_trn+XXt_cal, Xy=Xyt_trn+Xyt_cal, b=prior, learningRate=1/200, nIter=opt_nIter, returnPath=F)
getCor(XXt_tst, Xyt_tst, yyt_tst, fm_GDES_final)
#> [1] 0.1973647
```

- #### Transfer Learning using penalized regressions (*TL-PR*)

*PR()* function takes as inputs the sufficient statistics derived from the target population and a vector of initial values (prior), plus, potentially, values for $\lambda$ and $\alpha$ (if these are not provided, by default $\alpha$ is set equal to zero, i.e., L2 penalty, and the model is fitted over a grid of values of $\lambda$). The function returns estimates for a grid of values of $\lambda$, enabling users to select the optimal model based on cross-validation.

```R
fm_PR=PR(XX=XXt_trn, Xy=Xyt_trn, b=prior, alpha=0, nLambda=100, convThreshold=1e-4,
         maxIter=1000, returnPath=FALSE)
str(fm_PR)
#> List of 4
#>  $ B        : num [1:2450, 1:100] -0.000234 -0.000301 0.000232 -0.000919 -0.000123 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:2450] "SNP_1" "SNP_2" "SNP_3" "SNP_4" ...
#>   .. ..$ : chr [1:100] "lambda_10013023.2585" "lambda_9123493.8693" "lambda_8312987.8193" "lambda_7574484.893" ...
#>  $ lambda   : num [1:100] 10013023 9123494 8312988 7574485 6901589 ...
#>  $ alpha    : num 0
#>  $ conv_iter: num [1:100] 2 2 2 2 2 2 2 2 2 2 ...
```

Next, we evaluate the prediction accuracy in the calibrating set to select the optimal shrinkage parameter (lambda).

```R
Cor_PR=getCor(XXt_cal, Xyt_cal, yyt_cal, fm_PR$B)
plot(x=log(fm_PR$lambda), y=Cor_PR, xlab='log(lambda)', ylab='Prediction Corr.', pch=20, type='o')
opt_lambda=fm_PR$lambda[which.max(Cor_PR)]
```

<p align="left">
    <img src="https://github.com/QuantGen/GPTL/blob/main/man/plots/E1_PR.png" alt="Description" width="400">
</p>

We then re-estimate the PGS effects using the training set, with the optimal shrinkage parameter, and evaluate the final prediction accuracy in the testing set.

```R
fm_PR_final=PR(XX=XXt_train, Xy=Xyt_train, b=prior, alpha=0, lambda=opt_lambda, convThreshold=1e-4,
            maxIter=1000, returnPath=FALSE)
getCor(XXt_test, Xyt_test, yyt_test, fm_PR_final$B)
#> [1] 0.628841
```

- #### Transfer Learning using Bayesian model with an informative finite mixture prior (*TL-BMM*)

*BMM()* function takes as inputs the sufficient statistics derived from the target population, a matrix (B) whose columns contain the priors (one row per variant, one column per prior source of information), and parameters that control the algorithm. The function returns posterior means and posterior SDs for variant effects and other unknown parameters (including posterior ‘inclusion’ probabilities that link each variant effect to each of the components of the mixture).

*BMM()* only requires a single run of the algorithm (when the input **X'X** is dense) because regularization parameters and variant effects are jointly inferred from the posterior distribution. Thus, this method does not require calibrating regularization parameters. We estimate the PGS effects using the training set, and evaluate the final prediction accuracy in the testing set.

```R
fm_BMM=BMM(XX=XXt_train, Xy=Xyt_train, my=mean(yt_train), vy=var(yt_train), B=cbind(prior,0), n=nrow(Xt_train),
           nIter=12000, burnIn=2000, thin=5, verbose=FALSE)
str(fm_BMM)
#> List of 7
#>  $ b           : Named num [1:1279] -0.00942 0.01315 0.00791 0.00285 0.00434 ...
#>   ..- attr(*, "names")= chr [1:1279] "wPt.0538" "wPt.8463" "wPt.6348" "wPt.9992" ...
#>  $ POST.PROB   : num [1:1279, 1:2] 0.445 0.497 0.482 0.499 0.493 ...
#>  $ postMeanVarB: num [1:2] 0.00142 0.00149
#>  $ postProb    : num [1:2] 0.5 0.5
#>  $ samplesVarB : num [1:12000, 1:2] 0.000559 0.000546 0.00057 0.000528 0.000522 ...
#>  $ samplesB    : num [1:12000, 1:1279] 0.04545 0.01208 -0.04309 0.00146 0.01017 ...
#>  $ samplesVarE : num [1:12000] 0.643 0.623 0.933 0.647 0.693 ...
getCor(XXt_test, Xyt_test, yyt_test, fm_BMM$b)
#> [1] 0.6166595
