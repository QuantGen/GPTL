### GPTL using LD reference panels and GWAS results

In the following example, we use a demo data set (see this [link](https://github.com/QuantGen/GPTL/blob/main/man/DemoData_Preparation.md) for more details) to illustrate how to construct Polygenic Scores (PGS) using GPTL when using an LD reference panel (internal or external) and GWAS summary statistics.

A complete pipleine may include:

 1. [Loading the data](https://github.com/QuantGen/GPTL/blob/main/man/Example_LD_GWAS.md#1-data-loading).
 2. [Deriving sufficient statistics](https://github.com/QuantGen/GPTL/blob/main/man/Example_LD_GWAS.md#2-deriving-sufficient-statistics).
 3. [Estimating PGS using GPTL](). This includes calibrating model parameters, which is not strictly needed but it is a useful benchmark.
 4. [Evaluating PGS prediction accuracy]().

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
#> List of 4
#>  $ XX:Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>   .. ..@ i       : int [1:15960] 0 1 2 3 0 1 2 3 0 1 ...
#>   .. ..@ p       : int [1:1271] 0 4 8 12 16 20 24 28 32 38 ...
#>   .. ..@ Dim     : int [1:2] 1270 1270
#>   .. ..@ Dimnames:List of 2
#>   .. .. ..$ : chr [1:1270] "wPt.1171" "c.312549" "c.306034" "c.346957" ...
#>   .. .. ..$ : chr [1:1270] "wPt.1171" "c.312549" "c.306034" "c.346957" ...
#>   .. ..@ x       : num [1:15960] 124.95 45.5 5.02 1.68 45.5 ...
#>   .. ..@ factors : list()
#>  $ Xy: num [1:1270, 1] 37.97 -4.02 18.18 23.76 19.76 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:1270] "wPt.1171" "c.312549" "c.306034" "c.346957" ...
#>   .. ..$ : NULL
#>  $ n : num 253
#>  $ B :'data.frame':	1270 obs. of  1 variable:
#>   ..$ B: num [1:1270] 0.01417 0.00519 -0.0164 -0.00537 0.02633 ...
```

#### 3. PGS Estimation Using GPTL

- #### Transfer Learning using Gradient Descent with Early Stopping (*TL-GDES*)

*GD()* function takes as input the above derived sufficient statistics and prior effects. The function returns regression coefficient values over the gradient descent cycles.

```R
fm_GDES=GD(XX=SS$XX, Xy=SS$Xy, b=SS$B, learningRate=1/100, nIter=100, returnPath=T)
dim(fm_GDES)
#> [1] 1270  100
B_GDES=as.matrix(cbind(SS$B, fm_GDES))
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

We then evaluate the final prediction accuracy in the testing set, with the PGS effects estimated using the optimal iteration parameter.

```R
Cor_GDES=cor(GENO.Target[tst,]%*%B_GDES[,opt_nIter], PHENO.Target$y[tst])
Cor_GDES
#> [1] 0.6372282
```

- #### Transfer Learning using penalized regressions (*TL-PR*)

*PR()* function takes as inputs the above derived sufficient statistics and prior effects, plus, potentially, values for $\lambda$ and $\alpha$ (if these are not provided, by default $\alpha$ is set equal to zero and the model is fitted over a grid of values of $\lambda$). The function returns estimates for a grid of values of $\lambda$, enabling users to select the optimal model based on cross-validation.

```R
fm_PR=PR(XX=SS$XX, Xy=SS$Xy, b=SS$B, alpha=0, nLambda=100, convThreshold=1e-4,
         maxIter=1000, returnPath=FALSE)
str(fm_PR)
#> List of 4
#>  $ B        : num [1:1270, 1:100] 0.01425 0.00518 -0.01635 -0.00532 0.02637 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:1270] "wPt.1171" "c.312549" "c.306034" "c.346957" ...
#>   .. ..$ : chr [1:100] "lambda_933935.2028" "lambda_850966.9734" "lambda_775369.4128" "lambda_706487.7312" ...
#>  $ lambda   : num [1:100] 933935 850967 775369 706488 643725 ...
#>  $ alpha    : num 0
#>  $ conv_iter: num [1:100] 3 3 3 3 3 3 3 3 3 3 ...

B_PR=fm_PR$B
```

We evaluate the prediction accuracy in the calibrating set to select the optimal shrinkage parameter (lambda).

```R
Cors_PR=cor(GENO.Target[cal,]%*%B_PR, PHENO.Target$y[cal])
opt_lambda=fm_PR$lambda[which.max(Cors_PR)]
plot(x=log(fm_PR$lambda), y=Cors_PR, xlab='log(lambda)', ylab='Prediction Corr.', pch=20, type='o');abline(v=log(opt_lambda), lty=2)
```

<p align="left">
    <img src="https://github.com/QuantGen/GPTL/blob/main/man/plots/E2_PR.png" alt="Description" width="400">
</p>

We then evaluate the final prediction accuracy in the testing set, with the PGS effects estimated using the optimal shrinkage parameter.

```R
Cor_PR=cor(GENO.Target[tst,]%*%B_PR[,which.max(Cors_PR)], PHENO.Target$y[tst])
Cor_PR
#> [1] 0.6427939
```

- #### Transfer Learning using Bayesian model with an informative finite mixture prior (*TL-BMM*)

*BMM()* function takes as inputs the above derived sufficient statistics prior, a matrix (B) whose columns contain the priors (one row per variant, one column per prior source of information), and parameters that control the algorithm. The function returns posterior means and posterior SDs for variant effects and other unknown parameters (including posterior ‘inclusion’ probabilities that link each variant effect to each of the components of the mixture).

When using sparse **X'X** as inputs for *BMM()*, variance parameters cannot be properly updated due to ignoring the LD between LD blocks. Therefore, in this circumstance, we suggest fixing the error and effects variances (by setting function parameters *fixVarE=FALSE, fixVarB=rep(FALSE,nComp)*). The error and effects variances will be computed based on *R2* parameters and we can profile *R2* over a grid of values (e.g., {0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5}). The optimal *R2* value can be selected in a calibration set.

```R
R2s=c(0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)
B_BMM=matrix(nrow=nrow(SS$XX), ncol=length(R2s))
for (i in 1:length(R2s)) {
    fm_BMM=BMM(XX=SS$XX, Xy=SS$Xy, B=cbind(SS$B,0), n=SS$n, R2=R2s[i], my=mean(PHENO.Target$y), vy=var(PHENO.Target$y),
           nIter=12000, burnIn=2000, thin=5, fixVarE=TRUE, fixVarB=rep(TRUE,2), verbose=FALSE)
    B_BMM[,i]=fm_BMM$b
}
```

We evaluate the prediction accuracy in the calibration set to select the optimal *R2* parameter.

```R
Cors_BMM=cor(GENO.Target[cal,]%*%B_BMM, PHENO.Target$y[cal])
opt_R2=R2s[which.max(Cors_BMM)]
plot(x=R2s, y=Cors_BMM, xlab='R2', ylab='Prediction Corr.', pch=20, type='o');abline(v=opt_R2, lty=2)
```

<p align="left">
    <img src="https://github.com/QuantGen/GPTL/blob/main/man/plots/E2_BMM.png" alt="Description" width="400">
</p>

We then evaluate the final prediction accuracy in the testing set, with the PGS effects estimated using the optimal *R2* parameter.

```R
Cor_BMM=cor(GENO.Target[tst,]%*%B_BMM[,which.max(Cors_BMM)], PHENO.Target$y[tst])
Cor_BMM
#> [1] 0.6277782
```

[Back to Homepage](https://github.com/QuantGen/GPTL/blob/main/README.md)


