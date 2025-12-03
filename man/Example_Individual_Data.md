### GPTL using individual genetype and phenotype data

In the following example, we use a demo data set (see this [link](https://github.com/QuantGen/GPTL/blob/main/man/DemoData_Preparation.md) for more details) to illustrate how to calibrate Polygenic Scores (PGS) using GPTL when one has access to individual (genotype and phenotype) data.

A complete pipleine may include:

 1. [Loading the data](https://github.com/QuantGen/GPTL/blob/main/man/Example_Individual_Data.md#1-data-loading)
 2. [Estimating effects in the source population](https://github.com/QuantGen/GPTL/blob/main/man/Example_Individual_Data.md#2-single-ancestry-pgs-source-population). This may not be needed if one uses a pre-trained PGS as prior information.
 3. Calibrating a PGS using only training (and calibration) data from the target populations–this is not strictly needed but it is a useful benchmark to have.

#### 1. Data Loading

```R
library(GPTL)
data(Ind_DemoData)
```

The demo data set has genotype and phenotype data for the source (n=346 samples, 1270 variants) and target (n=253 samples, 1270 variants) populations. The objects loded are:

 - `GENO.Source`,
 - `GENO.Target`,
 - `PHENO.Source`,
 - `PHENO.Target`.

In addition to phenotypic data,  `PHENO.Target` includes a variable (`sets`) defining the training, calibration, and testing subsets.

```R
trn=which(PHENO.Target$sets=='trn')
cal=which(PHENO.Target$sets=='cal')
tst=which(PHENO.Target$sets=='tst')
```

#### 2. Single-ancestry PGS: Source population

We use the source population data to construct a cross-ancestry PGS, and use the target population tarining and calibrating sets to construct a within-ancestry PGS. We estimated effects using a Bayesian shrinkage estimation method (a Bayesian model with a Gaussian prior centered at zero, model ‘BRR’ in the **BGLR** R-package).

```R
library(BGLR)
fm_Cross=BGLR(y=PHENO.Source$y, ETA=list(list(X=GENO.Source, model="BRR")), nIter = 6000, burnIn = 1000, verbose = FALSE)
B_Cross=fm_Cross$ETA[[1]]$b

fm_Within=BGLR(y=PHENO.Target$y[c(trn, cal)], ETA=list(list(X=GENO.Target[c(trn, cal),], model="BRR")), nIter = 6000, burnIn = 1000, verbose = FALSE)
B_Within=fm_Within$ETA[[1]]$b
```

We evaluate the prediction accuracy in the testing set.

```R
Cor_Cross=cor(GENO.Target[tst,]%*%B_Cross, PHENO.Target$y[tst])
Cor_Cross
#> [1] 0.4656199
Cor_Within=cor(GENO.Target[tst,]%*%B_Within, PHENO.Target$y[tst])
Cor_Within
#> [1] 0.4411236
```

#### 3. PGS Estimation Using GPTL

- #### Transfer Learning using Gradient Descent with Early Stopping (*TL-GDES*)

*GD()* function takes as input the sufficient statistics (**X'X** and **X'y**) derived from the target population and a vector of initial values (effects estimated from the source population—`B_Cross`). The function returns regression coefficient values over the gradient descent cycles.

```R
X=scale(GENO.Target[trn,],center=TRUE,scale=FALSE)
XX=crossprod(X)
Xy=crossprod(X,PHENO.Target$y[trn])

fm_GDES=GD(XX=XX, Xy=Xy, b=B_Cross, learningRate=1/200, nIter=100, returnPath=T)
dim(fm_GDES)
#> [1] 1270  100
B_GDES=cbind(B_Cross, fm_GDES)
```

We evaluate the prediction accuracy in the calibrating set to select the optimal number of gradient descent cycles (nIter).

```R
Cors_GDES=cor(GENO.Target[cal,]%*%B_GDES, PHENO.Target$y[cal])
opt_nIter=which.max(Cors_GDES)-1
plot(Cors_GDES, xlab='iteration', ylab='Prediction Corr.', pch=20, type='o');abline(v=opt_nIter, lty=2)
```

<p align="left">
    <img src="https://github.com/QuantGen/GPTL/blob/main/man/plots/E1_GDES.png" alt="Description" width="400">
</p>

We then re-estimate the PGS effects using both the training and calibrating sets, with the optimal shrinkage parameter, and evaluate the final prediction accuracy in the testing set.

```R
X=scale(GENO.Target[c(trn,cal),],center=TRUE,scale=FALSE)
XX=crossprod(X)
Xy=crossprod(X,PHENO.Target$y[c(trn,cal)])

fm_GDES_final=GD(XX=XX, Xy=Xy, b=B_Cross, learningRate=1/100, nIter=opt_nIter, returnPath=F)
Cor_GDES=cor(GENO.Target[tst,]%*%fm_GDES_final, PHENO.Target$y[tst])
Cor_GDES
#> [1] 0.5717656
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
#>  $ B        : num [1:1270, 1:100] 0.0155 0.03296 -0.0029 0.00448 0.01094 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:1270] "wPt.1171" "c.312549" "c.306034" "c.346957" ...
#>   .. ..$ : chr [1:100] "lambda_156552.6478" "lambda_142644.9421" "lambda_129972.7586" "lambda_118426.3369" ...
#>  $ lambda   : num [1:100] 156553 142645 129973 118426 107906 ...
#>  $ alpha    : num 0
#>  $ conv_iter: num [1:100] 3 3 3 3 3 3 3 3 3 3 ...
B_PR=fm_PR$B
```

Next, we evaluate the prediction accuracy in the calibrating set to select the optimal shrinkage parameter (lambda).

```R
Cors_PR=cor(GENO.Target[cal,]%*%B_PR, PHENO.Target$y[cal])
opt_lambda=fm_PR$lambda[which.max(Cors_PR)]
plot(x=log(fm_PR$lambda), y=Cors_PR, xlab='log(lambda)', ylab='Prediction Corr.', pch=20, type='o');abline(v=log(opt_lambda), lty=2)
```

<p align="left">
    <img src="https://github.com/QuantGen/GPTL/blob/main/man/plots/E1_PR.png" alt="Description" width="400">
</p>

We then re-estimate the PGS effects using both the training and calibrating sets, with the optimal shrinkage parameter, and evaluate the final prediction accuracy in the testing set.

```R
X=scale(GENO.Target[c(trn,cal),],center=TRUE,scale=FALSE)
XX=crossprod(X)
Xy=crossprod(X,PHENO.Target$y[c(trn,cal)])

fm_PR_final=PR(XX=XX, Xy=Xy, b=B_Cross, alpha=0, lambda=opt_lambda, convThreshold=1e-4,
               maxIter=1000, returnPath=FALSE)
Cor_PR=cor(GENO.Target[tst,]%*%fm_PR_final$B, PHENO.Target$y[tst])
Cor_PR
#> [1] 0.629457
```

- #### Transfer Learning using Bayesian model with an informative finite mixture prior (*TL-BMM*)

*BMM()* function takes as inputs the sufficient statistics (**X'X** and **X'y**) derived from the target population, a matrix (B) whose columns contain the priors (one row per variant, one column per prior source of information), and parameters that control the algorithm. The function returns posterior means and posterior SDs for variant effects and other unknown parameters (including posterior ‘inclusion’ probabilities that link each variant effect to each of the components of the mixture).

*BMM()* only requires a single run of the algorithm (when the input **X'X** is dense) because regularization parameters and variant effects are jointly inferred from the posterior distribution. Thus, this method does not require calibrating regularization parameters. We estimate the PGS effects using the training and calibrating sets, and evaluate the final prediction accuracy in the testing set.

```R
X=scale(GENO.Target[c(trn,cal),],center=TRUE,scale=FALSE)
XX=crossprod(X)
Xy=crossprod(X,PHENO.Target$y[c(trn,cal)])

fm_BMM=BMM(XX=XX, Xy=Xy, my=mean(PHENO.Target$y[c(trn,cal)]), vy=var(PHENO.Target$y[c(trn,cal)]), B=cbind(B_Cross,0),
           n=nrow(GENO.Target[c(trn,cal),]), nIter=12000, burnIn=2000, thin=5, verbose=FALSE)
str(fm_BMM)
#> List of 7
#>  $ b           : Named num [1:1270] 5.70e-03 -5.89e-04 -2.80e-02 -6.81e-03 -7.02e-05 ...
#>   ..- attr(*, "names")= chr [1:1270] "wPt.1171" "c.312549" "c.306034" "c.346957" ...
#>  $ POST.PROB   : num [1:1270, 1:2] 0.484 0.435 0.498 0.494 0.494 ...
#>  $ postMeanVarB: num [1:2] 0.00194 0.00196
#>  $ postProb    : num [1:2] 0.5 0.5
#>  $ samplesVarB : num [1:12000, 1:2] 0.000531 0.000548 0.000556 0.000595 0.000679 ...
#>  $ samplesB    : num [1:12000, 1:1270] 0.013872 -0.011184 0.022909 -0.042009 -0.000111 ...
#>  $ samplesVarE : num [1:12000] 0.587 0.588 0.53 0.612 0.62 ...
Cor_BMM=cor(GENO.Target[tst,]%*%fm_BMM$b, PHENO.Target$y[tst])
Cor_BMM
#> [1] 0.6385827
```

#### 4. Prediction Accuracy Summary

| Method | Prediction Squared Corr. |
| --- | --- |
| cross-ancestry | 0.0564 |
| within-ancestry | 0.4012 |
| TL-GDES | 0.3269 |
| TL-PR | 0.3962 |
| TL-BMM | 0.4078 |



[Back to Homepage](https://github.com/QuantGen/GPTL/blob/main/README.md)
