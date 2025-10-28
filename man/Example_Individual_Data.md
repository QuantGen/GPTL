### GPTL using individual genetype and phenotype data

The following example illustrate how GPTL software works when one has access to individual genetype and phenotype data. Here we use the [wheat](https://doi.org/10.1104/pp.105.063438) data set collected from CIMMYT's Global Wheat Program, including 599 wheat lines genotype (1279 variants) and phenotype (average grain yield).

**1. Data Preparation**

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
Xs=scale(wheat.X[CLUSTER$cluster == 1,], center=TRUE, scale=FALSE);ys=y[CLUSTER$cluster == 1]
Xt=scale(wheat.X[CLUSTER$cluster == 2,], center=TRUE, scale=FALSE);yt=y[CLUSTER$cluster == 2]
```

We estimated prior effects from the source data set using a Bayesian shrinkage estimation method (a Bayesian model with a Gaussian prior centered at zero, model ‘BRR’ in the **BGLR** R-package). Alternatively, if only sufficient statistics (**X'X** and **X'y**) for the source data set are provided, one can use *BLRCross()* function in the **BGLR** R-package.

```R
ETA=list(list(X=Xs, model="BRR"))
fm=BGLR(y=ys, ETA = ETA, response_type = "gaussian", nIter = 12000, burnIn = 2000, verbose = FALSE)
prior=fm$ETA[[1]]$b
names(prior)=colnames(Xs)
```

We further split the target data set into (i) a training set (40%), (ii) a calibration set (40%), and (iii) a testing set (20%), and compute the sufficient statistics (**X'X** and **X'y**) for the each of sets.

```R
set.seed(1234)
sets=as.integer(as.factor(cut(runif(nrow(Xt)),breaks=c(0,quantile(runif(nrow(Xt)),prob=c(.4,.8)),1.1))))
Xt_train=Xt[sets==1,];yt_train=yt[sets==1];XXt_train=crossprod(Xt_train);Xyt_train=crossprod(Xt_train, yt_train)
Xt_cali=Xt[sets==2,];yt_cali=yt[sets==2];XXt_cali=crossprod(Xt_cali);Xyt_cali=crossprod(Xt_cali, yt_cali);yyt_cali=crossprod(yt_cali)
Xt_test=Xt[sets==3,];yt_test=yt[sets==3];XXt_test=crossprod(Xt_test);Xyt_test=crossprod(Xt_test, yt_test);yyt_test=crossprod(yt_test)
```

**2. PGS Estimation Using GPTL**

- #### Loading the package

```R
library(GPTL)
```

- #### Transfer Learning using Gradient Descent with Early Stopping (*TL-GDES*)

*GD()* function takes as input the sufficient statistics derived from the target population and a vector of initial values (prior). The function returns regression coefficient values over the gradient descent cycles.

```R
fm_GDES=GD(XX=XXt_train, Xy=Xyt_train, b=prior, learningRate=1/50, nIter=100, returnPath=T)
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
#> [1] 0.4318031
```

- #### Transfer Learning using penalized regressions (*TL-PR*)

*PR()* function takes as inputs the sufficient statistics derived from the target population and a vector of initial values (prior), plus, potentially, values for $\lambda$ and $\alpha$ (if these are not provided, by default $\alpha$ is set equal to zero and the model is fitted over a grid of values of $\lambda$). The function returns estimates for a grid of values of $\lambda$, enabling users to select the optimal model based on cross-validation.

```R
fm_PR=PR(XX=XXt_train, Xy=Xyt_train, b=prior, alpha=0, nLambda=100, conv_threshold=1e-4,
         maxIter=1000, returnPath=FALSE)
str(fm_PR)
#> List of 4
#>  $ B        : num [1:1279, 1:100] -0.030243 0.012444 0.028643 0.000624 -0.002359 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:1279] "wPt.0538" "wPt.8463" "wPt.6348" "wPt.9992" ...
#>   .. ..$ : chr [1:100] "lambda_219760.0169" "lambda_200237.1427" "lambda_182448.6268" "lambda_166240.3936" ...
#>  $ lambda   : num [1:100] 219760 200237 182449 166240 151472 ...
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
#> [1] 0.4202875
```

- #### Transfer Learning using Bayesian model with an informative finite mixture prior (*TL-BMM*)

*BMM()* function takes as inputs the sufficient statistics derived from the target population, a matrix (B) whose columns contain the priors (one row per variant, one column per prior source of information), and parameters that control the algorithm. The function returns posterior means and posterior SDs for variant effects and other unknown parameters (including posterior ‘inclusion’ probabilities that link each variant effect to each of the components of the mixture).

*BMM()* only requires a single run of the algorithm (when the input **X'X** is dense) because regularization parameters and variant effects are jointly inferred from the posterior distribution. Thus, this method does not require calibrating regularization parameters. We estimate the PGS effects using both the training and calibration sets, and evaluate the final prediction accuracy in the testing set.

```R
fm_BMM=BMM(XX=XXt_train+XXt_cali, Xy=Xyt_train+Xyt_cali, my=mean(c(yt_train,yt_cali)), vy=var(c(yt_train,yt_cali)), B=cbind(prior,0), n=nrow(Xt_train)+nrow(Xt_cali),
           nIter=12000, burnIn=2000, thin=5, verbose=FALSE)
str(fm_BMM)
#> List of 7
#>  $ b           : Named num [1:1279] -0.00229 0.0161 0.0066 0.0024 0.00101 ...
#>   ..- attr(*, "names")= chr [1:1279] "wPt.0538" "wPt.8463" "wPt.6348" "wPt.9992" ...
#>  $ POST.PROB   : num [1:1279, 1:2] 0.441 0.532 0.477 0.503 0.508 ...
#>  $ postMeanVarB: num [1:2] 0.00169 0.00185
#>  $ postProb    : num [1:2] 0.498 0.502
#>  $ samplesVarB : num [1:12000, 1:2] 0.000672 0.000675 0.00075 0.000775 0.000821 ...
#>  $ samplesB    : num [1:12000, 1:1279] -0.0232 -0.0401 -0.0131 -0.0179 -0.0306 ...
#>  $ samplesVarE : num [1:12000] 0.72 0.657 0.66 0.674 0.661 ...
getCor(XXt_test, Xyt_test, yyt_test, fm_BMM$b)
#> [1] 0.410943
```

[Back to Homepage](https://github.com/QuantGen/GPTL/blob/main/README.md)
