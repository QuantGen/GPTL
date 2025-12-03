### GPTL using LD reference panels and GWAS results

In the following example, we use a demo data set (see this [link](https://github.com/QuantGen/GPTL/blob/main/man/DemoData_Preparation.md) for more details) to illustrate how to construct Polygenic Scores (PGS) using GPTL when using an LD reference panel (internal or external) and GWAS summary statistics.

A complete pipleine may include:

 1. [Loading the data](https://github.com/QuantGen/GPTL/blob/main/man/Example_Individual_Data.md#1-data-loading).
 2. [Deriving sufficient statistics]().
 3. [Estimating PGS using GPTL](https://github.com/QuantGen/GPTL/blob/main/man/Example_Individual_Data.md#3-pgs-estimation-using-gptl). This includes calibrating model parameters, which is not strictly needed but it is a useful benchmark.
 4. [Evaluating PGS prediction accuracy](https://github.com/QuantGen/GPTL/blob/main/man/Example_Individual_Data.md#4-prediction-accuracy-summary).

#### 1. Data Loading

```R
library(GPTL)
data(Sum_DemoData)
```

The demo data set has prior effects estimated from a source population, a sparse LD reference matrix, GWAS summary statistics, and genotype and phenotype data (for the calibrating and testing purposes) from the target population. The objects loded are:

 - `PRIOR`,
 - `LD`,
 - `GWAS`,
 - `GENO.Target`,
 - `PHENO.Target`,

In addition to phenotypic data, `PHENO.Target` includes a variable (`sets`) defining the calibrating, and testing subsets.

```R
cal=which(PHENO.Target$sets=='cal')
tst=which(PHENO.Target$sets=='tst')
```

#### 2. Deriving Sufficient Statistics

We use *getSS()* function to compute the sufficient statistics (**X'X** and **X'y**) based on the LD matrix, GWAS statistics (and prior effects, to algin the variants).

```R
SS=getSS(ld=LD,gwas=GWAS,B=PRIOR)
#>  There were 1270 variant in common between the LD reference panel, the GWAS, and the prior.

str(SS)
#>  List of 4
#>   $ XX:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>    .. ..@ i       : int [1:30477] 0 1 2 3 0 1 2 3 0 1 ...
#>    .. ..@ p       : int [1:1280] 0 4 8 12 16 20 24 28 32 38 ...
#>    .. ..@ Dim     : int [1:2] 1279 1279
#>    .. ..@ Dimnames:List of 2
#>    .. .. ..$ : chr [1:1279] "wPt.1171" "c.312549" "c.306034" "c.346957" ...
#>    .. .. ..$ : chr [1:1279] "wPt.1171" "c.312549" "c.306034" "c.346957" ...
#>    .. ..@ x       : num [1:30477] 82.53 27.85 4.75 4.88 27.85 ...
#>    .. ..@ factors : list()
#>   $ Xy: num [1:1279, 1] -2.638 -10.77 -31.066 -21.927 0.166 ...
#>    ..- attr(*, "dimnames")=List of 2
#>    .. ..$ : chr [1:1279] "wPt.1171" "c.312549" "c.306034" "c.346957" ...
#>    .. ..$ : NULL
#>   $ n : num 167
#>   $ B :'data.frame':	1279 obs. of  1 variable:
#>    ..$ B: num [1:1279] 0.01324 0.03401 -0.00575 0.00472 0.00775 ...
```



#### 3. PGS Estimation Using GPTL

- #### Transfer Learning using Gradient Descent with Early Stopping (*TL-GDES*)

*GD()* function takes as input the sufficient statistics (**X'X** and **X'y**) derived from the target population and a vector of initial values (effects estimated from the source population—`B_Cross`). The function returns regression coefficient values over the gradient descent cycles.

```R
X=scale(GENO.Target[trn,],center=TRUE,scale=FALSE)
XX=crossprod(X)
Xy=crossprod(X,PHENO.Target$y[trn])

fm_GDES=GD(XX=XX, Xy=Xy, b=B_Cross, learningRate=1/100, nIter=100, returnPath=T)
dim(fm_GDES)
#> [1] 1270  100
B_GDES=cbind(B_Cross, fm_GDES)
```

We evaluate the prediction accuracy in the calibrating set to select the optimal number of gradient descent cycles (nIter).

```R
Cors_GDES=cor(GENO.Target[cal,]%*%B_GDES, PHENO.Target$y[cal])
opt_nIter=which.max(Cors_GDES)
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
#> [1] 0.5258604
```










**2. PGS Estimation Using GPTL**

- #### Transfer Learning using Gradient Descent with Early Stopping (*TL-GDES*)

*GD()* function takes as input the above derived sufficient statistics and prior effects. The function returns regression coefficient values over the gradient descent cycles.

```R
fm_GDES=GD(XX=SS$XX, Xy=SS$Xy, b=SS$B, learningRate=1/100, nIter=100, returnPath=T)
dim(fm_GDES)
#> [1] 1279  100
```

We evaluate the prediction accuracy in the calibration set to select the optimal shrinkage parameter (nIter).

```R
Cor_GDES=cor(wheat_VLD.X %*% fm_GDES, wheat_VLD.y)
plot(Cor_GDES, xlab='iteration', ylab='Prediction Corr.', pch=20)
opt_nIter=which.max(Cor_GDES)
```

<p align="left">
    <img src="https://github.com/QuantGen/GPTL/blob/main/man/plots/E2_GDES.png" alt="Description" width="400">
</p>

We then re-estimate the PGS effects using the training set, with the optimal shrinkage parameter, and evaluate the final prediction accuracy in the testing set.

```R
fm_GDES_final=GD(XX=SS$XX, Xy=SS$Xy, b=SS$B, learningRate=1/100, nIter=opt_nIter, returnPath=F)
cor(wheat_TST.X %*% fm_GDES_final, wheat_TST.y)
#> [1] 0.606333
```

- #### Transfer Learning using penalized regressions (*TL-PR*)

*PR()* function takes as inputs the above derived sufficient statistics and prior effects, plus, potentially, values for $\lambda$ and $\alpha$ (if these are not provided, by default $\alpha$ is set equal to zero and the model is fitted over a grid of values of $\lambda$). The function returns estimates for a grid of values of $\lambda$, enabling users to select the optimal model based on cross-validation.

```R
fm_PR=PR(XX=SS$XX, Xy=SS$Xy, b=SS$B, alpha=0, nLambda=100, convThreshold=1e-4,
         maxIter=1000, returnPath=FALSE)
str(fm_PR)
#> List of 4
#>  $ B        : num [1:1279, 1:100] 0.01322 0.03397 -0.00585 0.00464 0.00775 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:1279] "wPt.1171" "c.312549" "c.306034" "c.346957" ...
#>   .. ..$ : chr [1:100] "lambda_617943.6592" "lambda_563047.2476" "lambda_513027.682" "lambda_467451.7167" ...
#>  $ lambda   : num [1:100] 617944 563047 513028 467452 425925 ...
#>  $ alpha    : num 0
#>  $ conv_iter: num [1:100] 3 3 3 3 3 3 3 3 3 3 ...
```

We evaluate the prediction accuracy in the calibration set to select the optimal shrinkage parameter (lambda).

```R
Cor_PR=cor(wheat_VLD.X %*% fm_PR$B, wheat_VLD.y)
plot(x=log(fm_PR$lambda), y=Cor_PR, xlab='log(lambda)', ylab='Prediction Corr.', pch=20)
opt_lambda=fm_PR$lambda[which.max(Cor_PR)]
```

<p align="left">
    <img src="https://github.com/QuantGen/GPTL/blob/main/man/plots/E2_PR.png" alt="Description" width="400">
</p>

We then re-estimate the PGS effects using the training set, with the optimal shrinkage parameter, and evaluate the final prediction accuracy in the testing set.

```R
fm_PR_final=PR(XX=SS$XX, Xy=SS$Xy, b=SS$B, alpha=0, lambda=opt_lambda, convThreshold=1e-4,
            maxIter=1000, returnPath=FALSE)
cor(wheat_TST.X %*% fm_PR_final$B, wheat_TST.y)
#> [1] 0.5645487
```

- #### Transfer Learning using Bayesian model with an informative finite mixture prior (*TL-BMM*)

*BMM()* function takes as inputs the above derived sufficient statistics prior, a matrix (B) whose columns contain the priors (one row per variant, one column per prior source of information), and parameters that control the algorithm. The function returns posterior means and posterior SDs for variant effects and other unknown parameters (including posterior ‘inclusion’ probabilities that link each variant effect to each of the components of the mixture).

When using sparse **X'X** as inputs for *BMM()*, variance parameters cannot be properly updated due to ignoring the LD between LD blocks. Therefore, in this circumstance, we suggest fixing the error and effects variances (by setting function parameters *fixVarE=FALSE, fixVarB=rep(FALSE,nComp)*). The error and effects variances will be computed based on *R2* parameters and we can profile *R2* over a grid of values (e.g., {0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5}). The optimal *R2* value can be selected in a calibration set.

```R
R2s=c(0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)
B_BMM=matrix(nrow=nrow(SS$XX), ncol=length(R2s))
for (i in 1:length(R2s)) {
    fm_BMM=BMM(XX=SS$XX, Xy=SS$Xy, B=cbind(SS$B,0), n=SS$n, R2=R2s[i],
           nIter=12000, burnIn=2000, thin=5, fixVarE=TRUE, fixVarB=rep(TRUE,2), verbose=FALSE)
    B_BMM[,i]=fm_BMM$b
}
```

We evaluate the prediction accuracy in the calibration set to select the optimal *R2* parameter.

```R
Cor_BMM=cor(wheat_VLD.X %*% B_BMM, wheat_VLD.y)
plot(x=R2s, y=Cor_BMM, xlab='R2', ylab='Prediction Corr.', pch=20)
opt_R2=R2s[which.max(Cor_BMM)]
```

<p align="left">
    <img src="https://github.com/QuantGen/GPTL/blob/main/man/plots/E2_BMM.png" alt="Description" width="400">
</p>

We then re-estimate the PGS effects using the training set, with the optimal *R2* parameter, and evaluate the final prediction accuracy in the testing set.

```R
fm_BMM_final=BMM(XX=SS$XX, Xy=SS$Xy, B=cbind(SS$B,0), n=SS$n, R2=opt_R2,
           nIter=12000, burnIn=2000, thin=5, fixVarE=TRUE, fixVarB=rep(TRUE,2), verbose=FALSE)
cor(wheat_TST.X %*% fm_BMM_final$b, wheat_TST.y)
#> [1] 0.5363532
```

[Back to Homepage](https://github.com/QuantGen/GPTL/blob/main/README.md)


