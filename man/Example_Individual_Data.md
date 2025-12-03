### GPTL using individual genetype and phenotype data

In the following example, we use demo data generated from [link](https://github.com/QuantGen/GPTL/blob/main/man/DemoData_Preparation.md) to illustrate how GPTL software works when one has access to individual genetype and phenotype data.

**1. Data Loading**

```R
library(GPTL)
data(Ind_DemoData)
```

`GENO.Source` and `PHENO.Source` consist of genotype and phenotype data for the source population (346 samples, 1270 variants). `GENO.Target` and `PHENO.Target` consist of genotype and phenotype data for the target population (253 samples, 1270 variants), with samples splitted into 3 sets, marking in `PHENO.Target`.

To simplify the script below, we assign vectors for each of the sets.

```R
trn=which(PHENO.Target$sets=='trn')
cal=which(PHENO.Target$sets=='cal')
tst=which(PHENO.Target$sets=='tst')
```

**2. Single-ancestry PGS**

We use the source population data to construct a cross-ancestry PGS, and use the target population tarining and calibrating sets to construct a within-ancestry PGS. We estimated effects using a Bayesian shrinkage estimation method (a Bayesian model with a Gaussian prior centered at zero, model ‘BRR’ in the **BGLR** R-package).

```R
library(BGLR)
ETA=list(list(X=GENO.Source, model="BRR"))
fm_Cross=BGLR(y=PHENO.Source, ETA = ETA, response_type = "gaussian", nIter = 12000, burnIn = 2000, verbose = FALSE)
B_Cross=fm_Cross$ETA[[1]]$b

ETA=list(list(X=GENO.Target[c(trn, cal),], model="BRR"))
fm_Within=BGLR(y=PHENO.Target$y[c(trn, cal)], ETA = ETA, response_type = "gaussian", nIter = 12000, burnIn = 2000, verbose = FALSE)
B_Within=fm_Within$ETA[[1]]$b
```

We evaluate the prediction accuracy in the testing set.

```R
Cor_Cross=cor(GENO.Target[tst,]%*%B_Cross, PHENO.Target$y[tst])
Cor_Within=cor(GENO.Target[tst,]%*%B_Within, PHENO.Target$y[tst])
```

**3. PGS Estimation Using GPTL**

- #### Transfer Learning using Gradient Descent with Early Stopping (*TL-GDES*)

*GD()* function takes as input the sufficient statistics (**X'X** and **X'y**) derived from the target population and a vector of initial values (effects estimated from the source population—`B_Cross`). The function returns regression coefficient values over the gradient descent cycles.

```R
X=scale(GENO.Target[trn,],center=TRUE,scale=FALSE)
XX=crossprod(X)
Xy=crossprod(X,PHENO.Target$y[trn])

fm_GDES=GD(XX=XX, Xy=Xy, b=B_Cross, learningRate=1/50, nIter=100, returnPath=T)
dim(fm_GDES)
#> [1] 1270  100
B_GDES=cbind(B_Cross, fm_GDES)
```

We evaluate the prediction accuracy in the calibrating set to select the optimal number of gradient descent cycles (nIter).

```R
Cors_GDES=cor(GENO.Target[cal,]%*%B_GDES, PHENO.Target$y[cal])
plot(Cors_GDES, xlab='iteration', ylab='Prediction Corr.', pch=20, type='o')
opt_nIter=which.max(Cors_GDES)-1
```

<p align="left">
    <img src="https://github.com/QuantGen/GPTL/blob/main/man/plots/E1_GDES_wheat.png" alt="Description" width="400">
</p>

We then re-estimate the PGS effects using both the training and calibrating sets, with the optimal shrinkage parameter, and evaluate the final prediction accuracy in the testing set.

```R
X=scale(GENO.Target[c(trn,cal),],center=TRUE,scale=FALSE)
XX=crossprod(X)
Xy=crossprod(X,PHENO.Target$y[c(trn,cal)])

fm_GDES_final=GD(XX=XX, Xy=Xy, b=B_Cross, learningRate=1/500, nIter=opt_nIter, returnPath=F)
Cor_GDES=cor(GENO.Target[tst,]%*%fm_GDES_final, PHENO.Target$y[tst])
Cor_GDES
#> [1] 0.5955236
```

- #### Transfer Learning using penalized regressions (*TL-PR*)

*PR()* function takes as inputs the sufficient statistics (**X'X** and **X'y**) derived from the target population and a vector of initial values (effects estimated from the source population—`B_Cross`), plus, potentially, values for $\lambda$ and $\alpha$ (if these are not provided, by default $\alpha$ is set equal to zero, i.e., L2 penalty, and the model is fitted over a grid of values of $\lambda$). The function returns estimates for a grid of values of $\lambda$, enabling users to select the optimal model based on cross-validation.

```R
X=scale(GENO.Target[trn,],center=TRUE,scale=FALSE)
XX=crossprod(X)
Xy=crossprod(X,PHENO.Target$y[trn])

fm_PR=PR(XX=XX, Xy=Xy, b=B_Cross, alpha=0, nLambda=100, convThreshold=1e-4,
         maxIter=1000, returnPath=FALSE)
str(fm_PR)
#> List of 4
#>  $ B        : num [1:1270, 1:100] 0.01348 0.03301 -0.00341 0.00475 0.00989 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:1270] "wPt.1171" "c.312549" "c.306034" "c.346957" ...
#>   .. ..$ : chr [1:100] "lambda_155558.9764" "lambda_141739.5457" "lambda_129147.7951" "lambda_117674.6609" ...
#>  $ lambda   : num [1:100] 155559 141740 129148 117675 107221 ...
#>  $ alpha    : num 0
#>  $ conv_iter: num [1:100] 3 3 3 3 3 3 3 3 3 3 ...
B_PR=fm_PR$B
```

Next, we evaluate the prediction accuracy in the calibrating set to select the optimal shrinkage parameter (lambda).

```R
Cors_PR=getCor(XXt_cali, Xyt_cali, yyt_cali, fm_PR$B)
plot(x=log(fm_PR$lambda), y=Cor_PR, xlab='log(lambda)', ylab='Prediction Corr.', pch=20)
opt_lambda=fm_PR$lambda[which.max(Cor_PR)]

Cors_GDES=cor(GENO.Target[tst,]%*%B_GDES, PHENO.Target$y[tst])
plot(Cors_GDES, xlab='iteration', ylab='Prediction Corr.', pch=20, type='o')
opt_nIter=which.max(Cors_GDES)-1
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
