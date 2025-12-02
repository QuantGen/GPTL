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
fm=BLRXy(y=ys,ETA=LP,verbose=FALSE,nIter=6000,burnIn=1000)
prior=fm$ETA[[1]]$b
names(prior)=colnames(Xs)
```

**3. PGS Estimation Using GPTL**

- #### Transfer Learning using Gradient Descent with Early Stopping (*TL-GDES*)

*GD()* function takes as input the sufficient statistics derived from the target population and a vector of initial values (prior). The function returns regression coefficient values over the gradient descent cycles.

```R
fm_GDES=GD(XX=XXt_trn, Xy=Xyt_trn, b=prior, learningRate=1/1000, nIter=100, returnPath=T)
dim(fm_GDES)
#> [1] 2450  100
```

We evaluate the prediction accuracy in the calibrating set to select the optimal number of gradient descent cycles (nIter).

```R
Cor_GDES=getCor(XXt_cal, Xyt_cal, yyt_cal, fm_GDES)
opt_nIter=which.max(Cor_GDES)
plot(Cor_GDES, xlab='iteration', ylab='Prediction Corr.', pch=20, type='o');abline(v=opt_nIter, lty=2)
```

<p align="left">
    <img src="https://github.com/QuantGen/GPTL/blob/main/man/plots/E1_GDES.png" alt="Description" width="400">
</p>

We then re-estimate the PGS effects using both the training and calibrating sets, with the optimal shrinkage parameter, and evaluate the final prediction accuracy in the testing set.

```R
fm_GDES_final=GD(XX=XXt_trn+XXt_cal, Xy=Xyt_trn+Xyt_cal, b=prior, learningRate=1/1000, nIter=opt_nIter, returnPath=F)
getCor(XXt_tst, Xyt_tst, yyt_tst, fm_GDES_final)
#> [1] 0.2472281
```

- #### Transfer Learning using penalized regressions (*TL-PR*)

*PR()* function takes as inputs the sufficient statistics derived from the target population and a vector of initial values (prior), plus, potentially, values for $\lambda$ and $\alpha$ (if these are not provided, by default $\alpha$ is set equal to zero, i.e., L2 penalty, and the model is fitted over a grid of values of $\lambda$). The function returns estimates for a grid of values of $\lambda$, enabling users to select the optimal model based on cross-validation.

```R
fm_PR=PR(XX=XXt_trn, Xy=Xyt_trn, b=prior, alpha=0, nLambda=100, convThreshold=1e-4,
         maxIter=1000, returnPath=FALSE)
str(fm_PR)
#> List of 4
#>  $ B        : num [1:2450, 1:100] 8.09e-05 1.21e-04 -8.16e-05 3.76e-04 -3.12e-04 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:2450] "SNP_1" "SNP_2" "SNP_3" "SNP_4" ...
#>   .. ..$ : chr [1:100] "lambda_10079316.6912" "lambda_9183897.9761" "lambda_8368025.7918" "lambda_7624633.4437" ...
#>  $ lambda   : num [1:100] 10079317 9183898 8368026 7624633 6947282 ...
#>  $ alpha    : num 0
#>  $ conv_iter: num [1:100] 2 2 2 2 2 2 2 2 2 2 ...
```

Next, we evaluate the prediction accuracy in the calibrating set to select the optimal shrinkage parameter (lambda).

```R
Cor_PR=getCor(XXt_cal, Xyt_cal, yyt_cal, fm_PR$B)
opt_lambda=fm_PR$lambda[which.max(Cor_PR)]
plot(x=log(fm_PR$lambda), y=Cor_PR, xlab='log(lambda)', ylab='Prediction Corr.', pch=20, type='o');abline(v=log(opt_lambda), lty=2)
```

<p align="left">
    <img src="https://github.com/QuantGen/GPTL/blob/main/man/plots/E1_PR.png" alt="Description" width="400">
</p>

We then re-estimate the PGS effects using both the training and calibrating sets, with the optimal shrinkage parameter, and evaluate the final prediction accuracy in the testing set.

```R
fm_PR_final=PR(XX=XXt_trn+XXt_cal, Xy=Xyt_trn+Xyt_cal, b=prior, alpha=0, lambda=opt_lambda, convThreshold=1e-4,
            maxIter=1000, returnPath=FALSE)
getCor(XXt_tst, Xyt_tst, yyt_tst, fm_PR_final$B)
#> [1] 0.2472293
```

- #### Transfer Learning using Bayesian model with an informative finite mixture prior (*TL-BMM*)

*BMM()* function takes as inputs the sufficient statistics derived from the target population, a matrix (B) whose columns contain the priors (one row per variant, one column per prior source of information), and parameters that control the algorithm. The function returns posterior means and posterior SDs for variant effects and other unknown parameters (including posterior ‘inclusion’ probabilities that link each variant effect to each of the components of the mixture).

*BMM()* only requires a single run of the algorithm (when the input **X'X** is dense) because regularization parameters and variant effects are jointly inferred from the posterior distribution. Thus, this method does not require calibrating regularization parameters. We estimate the PGS effects using both the training and calibrating sets, and evaluate the final prediction accuracy in the testing set.

```R
fm_BMM=BMM(XX=XXt_trn+XXt_cal, Xy=Xyt_trn+Xyt_cal, my=mean(c(yt_trn, yt_cal)), vy=var(c(yt_trn, yt_cal)), B=cbind(prior,0),
           n=nrow(Xt_trn)+nrow(Xt_cal), nIter=6000, burnIn=1000, thin=5, verbose=FALSE)
str(fm_BMM)
#> List of 7
#>  $ b           : Named num [1:2450] 0.003573 0.000378 0.000174 -0.003342 0.00017 ...
#>   ..- attr(*, "names")= chr [1:2450] "SNP_1" "SNP_2" "SNP_3" "SNP_4" ...
#>  $ POST.PROB   : num [1:2450, 1:2] 0.513 0.475 0.5 0.477 0.53 0.482 0.496 0.512 0.513 0.499 ...
#>  $ postMeanVarB: num [1:2] 7.27e-05 7.77e-05
#>  $ postProb    : num [1:2] 0.505 0.495
#>  $ samplesVarB : num [1:6000, 1:2] 0.000135 0.000135 0.00014 0.000132 0.000138 ...
#>  $ samplesB    : num [1:6000, 1:2450] 0.01657 -0.00717 0.01233 0.00897 -0.00304 ...
#>  $ samplesVarE : num [1:6000] 0.837 0.839 0.852 0.841 0.815 ...
getCor(XXt_tst, Xyt_tst, yyt_tst, fm_BMM$b)
#> [1] 0.273657
