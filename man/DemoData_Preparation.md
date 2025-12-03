### Simulating genotype and phenotype data as examples for GPTL

The following script generates demo datasets to illustrate how GPTL software works when one has access to individual genetype and phenotype data or using an LD reference panel (internal or external) and GWAS summary statistics. The generated data sets are stored as **.RData** formats:

- `Ind_DemoData.RData`, including individual genotype and phenotype data for a source and a target population, with the target population data being subset into training, calibrating, and testing sets.

- `Sum_DemoData.RData`, including prior effects estimated from a source population, and LD reference panel, GWAS summary statistics, and individual calibrating/testing data from a target population.

Here we use the [wheat](https://doi.org/10.1104/pp.105.063438) data set collected from CIMMYT's Global Wheat Program, including 599 wheat lines genotype (1279 variants) and phenotype (average grain yield).

This data set has two clear clusters, we use this to illustrate how to transfer learning from one source population to improve prediction accuracy in another target population.

```R
library(BGLR)
data(wheat)
y=wheat.Y[,1]
X=scale(wheat.X, center=TRUE, scale=FALSE)

CLUSTER=kmeans(X,centers=2,nstart=100)
table(CLUSTER$cluster)
#>   1   2 
#> 346 253 
```

We use samples in cluster 1 as the source data set (where information is transferred) and samples in cluster 2 as the target data set (where the PGS will be used). 

```R
GENO.Source=scale(wheat.X[CLUSTER$cluster == 1,], center=TRUE, scale=FALSE)
GENO.Target=scale(wheat.X[CLUSTER$cluster == 2,], center=TRUE, scale=FALSE)
PHENO.Source=wheat.Y[CLUSTER$cluster == 1,1]
PHENO.Target=wheat.Y[CLUSTER$cluster == 2,1]
```

We further split the target data set into (i) a training set (40%), (ii) a calibrating set (30%), and (iii) a testing set (30%).

```R
set.seed(10)
sets <- cut(runif(nrow(GENO.Target)), breaks = c(0, 0.4, 0.7, 1), labels = c("trn","cal","tst"))
PHENO.Target=data.frame(y=PHENO.Target, sets=sets)
```

`GENO.Source` and `PHENO.Source` consist of genotype and phenotype data for the source population (346 samples, 1279 variants). `GENO.Target` and `PHENO.Target` consist of genotype and phenotype data for the target population (253 samples, 1279 variants), with samples splitted into 3 sets, marking in `PHENO.Target`.

The above-generated data sets are saved in `Ind_DemoData.RData`, and can be loaded by `data(Ind_DemoData)` with the `GPTL` package.

```R
save(GENO.Source, PHENO.Source, GENO.Target, PHENO.Target, file='Ind_SimData.RData')
```

We continue to prepare demo data 2. We estimate prior effects from the source population using a Bayesian shrinkage estimation method (a Bayesian model with a Gaussian prior centered at zero, model ‘BRR’ in the **BGLR** R-package). Alternatively, if only sufficient statistics (**X'X** and **X'y**) for the source data set are provided, one can use *BLRCross()* function in the **BGLR** R-package.

```R
ETA=list(list(X=GENO.Source, model="BRR"))
fm=BGLR(y=PHENO.Source, ETA = ETA, response_type = "gaussian", nIter = 12000, burnIn = 2000, verbose = FALSE)
PRIOR=fm$ETA[[1]]$b
```

We define LD block boundaries and generate the sparse LD reference matrix for the target population.

```R
source()

```














The following script simulates genotype and phenotype data that are used as examples. The simulated data sets are stored as **.RData** formats:

- **Ind_SimData.RData**, including individual genotype and phenotype data for a source and a target population.

- **Sum_SimData.RData**, including prior effects estimated from a source population, and LD reference panel, GWAS results, and individual calibrating/testing data for a target population.

Setting simulation parameters.

```R
library(GPTL)

nChr=5 # number of chromosomes
p=500  # number of loci per chromosome
nSource=10000 # sample size for the source populations
nTarget=4000 # sample size for the target pouplation
nQTN=50 # number of causal variants
h2=0.3
 
set.seed(1950)
```

Simulating gentoypes for the source population.

```R
X1=matrix(nrow=nSource,ncol=0)
for(i in 1:nChr){
    DPrime=rbeta(shape1=20,shape2=10,n=p) # DPrime parameter to control LD
    # shape1 and shape2 are used control the distribution of allele frequencies;
    X1=cbind(X1,sampleGenomes(n=nSource,nLoci=p,Dprime=DPrime,shape1=2,shape2=6))
}
```

Simulating gentoypes for the target population. Relative to source, in target we simulate similar LD with a different pattern.

```R
X2=matrix(nrow=nTarget,ncol=0)
for(i in 1:nChr){
    DPrime=rbeta(shape1=20,shape2=10,n=p)
    X2=cbind(X2,sampleGenomes(n=nTarget,nLoci=p,Dprime=DPrime,shape1=2,shape2=6))
}
```

Simulating genetic values and phenotypes.

```R
bSource=rgamma(nQTN,shape=20,scale=1)
bTarget=bSource/2+rgamma(nQTN,shape=20,scale=1)/2
 
QTN=as.integer(seq(from=50,to=p*nChr,length=nQTN))
gSource=X1[,QTN]%*%bSource
gTarget=X2[,QTN]%*%bTarget
 
SDg=sd(gSource)
gSource=scale(gSource,center=TRUE,scale=FALSE)*(sqrt(h2)/SDg)
gTarget=scale(gTarget,center=TRUE,scale=FALSE)*(sqrt(h2)/SDg)
 
ys=gSource+rnorm(n=nSource,sd=sqrt(1-h2))
yt=gTarget+rnorm(n=nTarget,sd=sqrt(1-h2))
 
# Removing causal variants
Xs=X1[,-QTN]
Xt=X2[,-QTN]
colnames(Xs)=paste0('SNP_', 1:ncol(Xs))
colnames(Xt)=paste0('SNP_', 1:ncol(Xt))
```

Splitting the target data into training, calibrating, and testing sets.

```R
sets=c(rep('trn', 3000), rep('cal', 500), rep('tst', 500))
Xt_trn=Xt[sets=='trn',];Xt_cal=Xt[sets=='cal',];Xt_tst=Xt[sets=='tst',]
yt_trn=yt[sets=='trn'];yt_cal=yt[sets=='cal'];yt_tst=yt[sets=='tst']
```

Saving above-generated data sets.

```R
save(Xs, ys, Xt_trn, Xt_cal, Xt_tst, yt_trn, yt_cal, yt_tst, file='Ind_SimData.RData')
```







We use a toy data set wheatSumStats of LD matrix wheat_LD, GWAS results wheat_GWAS, prior effects wheat_PRIOR, and individual validation/testing data wheat_VLD.X, wheat_VLD.y, wheat_TST.X, wheat_TST.y. These statistics were generated based on the wheat data set collected from CIMMYT's Global Wheat Program, including 599 wheat lines genotype (1279 variants) and phenotype (average grain yield).


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

We further split the target data set into (i) a training set (60%), (ii) a calibration set (20%), and (iii) a testing set (20%), and compute the sufficient statistics (**X'X** and **X'y**) for the each of sets.

```R
set.seed(1234)
sets=as.integer(as.factor(cut(runif(nrow(Xt)),breaks=c(0,quantile(runif(nrow(Xt)),prob=c(.6,.8)),1.1))))
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

We evaluate the prediction accuracy in the calibration set to select the optimal number of gradient descent cycles (nIter).

```R
Cor_GDES=getCor(XXt_cali, Xyt_cali, yyt_cali, fm_GDES)
plot(Cor_GDES, xlab='iteration', ylab='Prediction Corr.', pch=20)
opt_nIter=which.max(Cor_GDES)
```

<p align="left">
    <img src="https://github.com/QuantGen/GPTL/blob/main/man/plots/E1_GDES.png" alt="Description" width="400">
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
```

[Back to Homepage](https://github.com/QuantGen/GPTL/blob/main/README.md)
