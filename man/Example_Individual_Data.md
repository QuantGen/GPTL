### GPTL using individual genetype and phenotype data

In the following example, we use demo data generated from [link](https://github.com/QuantGen/GPTL/blob/main/man/DemoData_Preparation.md) to illustrate how GPTL software works when one has access to individual genetype and phenotype data.

**1. Data Loading**

```R
library(GPTL)
data(Ind_DemoData)
```

`GENO.Source` and `PHENO.Source` consist of genotype and phenotype data for the source population (346 samples, 1270 variants). `GENO.Target` and `PHENO.Target` consist of genotype and phenotype data for the target population (253 samples, 1270 variants), with samples splitted into 3 sets, marking in `PHENO.Target`.

**2. Single-ancestry PGS**

We use the source population data to construct a cross-ancestry PGS, and use the target population tarining and calibrating sets to construct a within-ancestry PGS. We estimated effects using a Bayesian shrinkage estimation method (a Bayesian model with a Gaussian prior centered at zero, model ‘BRR’ in the **BGLR** R-package).

```R
library(BGLR)
ETA=list(list(X=GENO.Source, model="BRR"))
fm_Cross=BGLR(y=PHENO.Source, ETA = ETA, response_type = "gaussian", nIter = 12000, burnIn = 2000, verbose = FALSE)
B_Cross=fm_Cross$ETA[[1]]$b

ETA=list(list(X=GENO.Target[PHENO.Target$sets %in% c('trn', 'cal'),], model="BRR"))
fm_Within=BGLR(y=PHENO.Target$y[PHENO.Target$sets %in% c('trn', 'cal')], ETA = ETA, response_type = "gaussian", nIter = 12000, burnIn = 2000, verbose = FALSE)
B_Within=fm_Within$ETA[[1]]$b
```

We evaluate the prediction accuracy in the testing set.

```R
Cor_Cross=cor(GENO.Target[PHENO.Target$sets=='tst',]%*%B_Cross, PHENO.Target$y[PHENO.Target$sets=='tst'])
Cor_Within=cor(GENO.Target[PHENO.Target$sets=='tst',]%*%B_Within, PHENO.Target$y[PHENO.Target$sets=='tst'])
```

**3. PGS Estimation Using GPTL**

- #### Transfer Learning using Gradient Descent with Early Stopping (*TL-GDES*)

We first compute the sufficient statistics (*X'X* and *X'y*)

*GD()* function takes as input the sufficient statistics derived from the target population and a vector of initial values (prior). The function returns regression coefficient values over the gradient descent cycles.

```R
fm_GDES=GD(XX=XXt_train, Xy=Xyt_train, b=prior, learningRate=1/50, nIter=100, returnPath=T)
dim(fm_GDES)
#> [1] 1279  100
```

We evaluate the prediction accuracy in the calibration set to select the optimal number of gradient descent cycles (nIter).

```R
Cor_GDES=getCor(XXt_cali, Xyt_cali, yyt_cali, fm_GDES)
plot(Cor_GDES, xlab='iteration', ylab='Prediction Corr.', pch=20)
opt_nIter=which.max(Cor_GDES)
```

<p align="left">
    <img src="https://github.com/QuantGen/GPTL/blob/main/man/plots/E1_GDES_wheat.png" alt="Description" width="400">
</p>

We then re-estimate the PGS effects using the training set, with the optimal shrinkage parameter, and evaluate the final prediction accuracy in the testing set.

```R
fm_GDES_final=GD(XX=XXt_train, Xy=Xyt_train, b=prior, learningRate=1/50, nIter=opt_nIter, returnPath=F)
getCor(XXt_test, Xyt_test, yyt_test, fm_GDES_final)
#> [1] 0.6364298
```

- #### Transfer Learning using penalized regressions (*TL-PR*)

*PR()* function takes as inputs the sufficient statistics derived from the target population and a vector of initial values (prior), plus, potentially, values for $\lambda$ and $\alpha$ (if these are not provided, by default $\alpha$ is set equal to zero, i.e., L2 penalty, and the model is fitted over a grid of values of $\lambda$). The function returns estimates for a grid of values of $\lambda$, enabling users to select the optimal model based on cross-validation.

```R
fm_PR=PR(XX=XXt_train, Xy=Xyt_train, b=prior, alpha=0, nLambda=100, convThreshold=1e-4,
         maxIter=1000, returnPath=FALSE)
str(fm_PR)
#> List of 4
#>  $ B        : num [1:1279, 1:100] 0.015507 0.017122 -0.007943 -0.000096 0.003377 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:1279] "wPt.0538" "wPt.8463" "wPt.6348" "wPt.9992" ...
#>   .. ..$ : chr [1:100] "lambda_256131.6607" "lambda_233377.6299" "lambda_212645.0045" "lambda_193754.2083" ...
#>  $ lambda   : num [1:100] 256132 233378 212645 193754 176542 ...
#>  $ alpha    : num 0
#>  $ conv_iter: num [1:100] 3 3 3 3 3 3 3 3 3 3 ...
```

Next, we evaluate the prediction accuracy in the calibration set to select the optimal shrinkage parameter (lambda).

```R
Cor_PR=getCor(XXt_cali, Xyt_cali, yyt_cali, fm_PR$B)
plot(x=log(fm_PR$lambda), y=Cor_PR, xlab='log(lambda)', ylab='Prediction Corr.', pch=20)
opt_lambda=fm_PR$lambda[which.max(Cor_PR)]
```

<p align="left">
    <img src="https://github.com/QuantGen/GPTL/blob/main/man/plots/E1_PR_wheat.png" alt="Description" width="400">
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
```

[Back to Homepage](https://github.com/QuantGen/GPTL/blob/main/README.md)
