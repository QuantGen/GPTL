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

#### Loading the package

```R
library(GPTL)
```

#### Transfer Learning using Gradient Descent with Early Stopping (*TL-GDES*)

GD() function takes as input the sufficient statistics derived from the target population and a vector of initial values (prior). The function returns regression coefficient values over the gradient descent cycles.

```R
fm_GDES=GD(XX=XXt_train, Xy=Xyt_train, b=prior, learning_rate=1/50, nIter=100, returnPath=T)
dim(fm_GDES)
#> [1] 1279  100
```

We evaluate the prediction accuracy in the calibration set to select the optimal shrinkage parameter (nIter).

```R
Cor_GDES=getCor(XXt_cali, Xyt_cali, yyt_cali, fm_GDES)
plot(Cor_GDES, xlab='iteration', ylab='Prediction Corr.', pch=20)
opt_nIter=which.max(Cor_GDES)
```

[
    <img
        src="MY_SRC_HERE" 
        width=70%
        title="My Image"
        alt="My Image"
    />
](https://github.com/QuantGen/GPTL/blob/main/man/plots/GDES_plot.png)

We then re-estimate the PGS effects using both the training and calibration sets, with the optimal shrinkage parameter, and evaluate the final prediction accuracy in the testing set.

```R
fm_GDES_final=GD(XX=XXt_train+XXt_cali, Xy=Xyt_train+Xyt_cali, b=prior, learning_rate=1/50, nIter=opt_nIter, returnPath=F)
getCor(XXt_test, Xyt_test, yyt_test, fm_GDES_final)
#> [1] 0.4334093
```




#### Transfer Learning using penalized regressions (*TL-PR*)

PR.SS() function takes as inputs the sufficient statistics derived from the target population and a vector of initial values (prior), plus, potentially, values for $\lambda$ and $\alpha$ (if these are not provided, by default $\alpha$ is set equal to zero and the model is fitted over a grid of values of $\lambda$). The function returns estimates for a grid of values of $\lambda$ and $\alpha$, enabling users to select the optimal model based on cross-validation.

```R
fm_PR=PR(XX=XXt_train, Xy=Xyt_train, b=prior, alpha=0, nLambda=100, conv_threshold=1e-4,
         maxIter=1000, returnPath=FALSE)
str(fm_PR)
#> List of 4
#>  $ B        : num [1:1279, 1:100] 0.002611 0.003227 -0.00547 -0.00045 -0.000568 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:1279] "wPt.0538" "wPt.8463" "wPt.6348" "wPt.9992" ...
#>   .. ..$ : chr [1:100] "lambda_4614.6071" "lambda_4407.4852" "lambda_4209.5668" "lambda_4020.443" ...
#>  $ lambda   : num [1:100] 4615 4407 4210 4020 3840 ...
#>  $ alpha    : num 0
#>  $ conv_iter: num [1:100] 10 10 10 10 11 11 12 12 12 13 ...
```

```R
Cor_PR=getCor(XXt_cali, Xyt_cali, yyt_cali, fm_PR$B)

```


#### Transfer Learning using Bayesian model with an informative finite mixture prior (*TL-BMM*)

BMM.SS() function takes as inputs the sufficient statistics derived from the target population, a matrix (B) whose columns contain the priors (one row per variant, one column per prior source of information), and parameters that control the algorithm. The function returns posterior means and posterior SDs for variant effects and other unknown parameters (including posterior ‘inclusion’ probabilities that link each variant effect to each of the components of the mixture).

```R
fm_BMM=BMM(XX=XXt_train, Xy=Xyt_train, my=mean(yt_train), vy=var(yt_train), B=cbind(prior,0), n=nrow(Xt_train),
           nIter=12000, burnIn=2000, thin=5, verbose=FALSE)
str(fm_BMM)
#> List of 7
#>  $ b           : Named num [1:1279] -0.00131 0.01974 0.00557 0.00333 -0.00142 ...
#>   ..- attr(*, "names")= chr [1:1279] "wPt.0538" "wPt.8463" "wPt.6348" "wPt.9992" ...
#>  $ POST.PROB   : num [1:1279, 1:2] 0.429 0.537 0.482 0.482 0.49 ...
#>  $ postMeanVarB: num [1:2] 0.00159 0.00194
#>  $ postProb    : num [1:2] 0.499 0.501
#>  $ samplesVarB : num [1:12000, 1:2] 0.000703 0.000793 0.000812 0.0009 0.000918 ...
#>  $ samplesB    : num [1:12000, 1:1279] 0.00358 -0.07363 -0.01596 -0.02647 -0.06138 ...
```

```R
Cor_BMM=getCor(XXt_cali, Xyt_cali, yyt_cali, fm_BMM$b)

```















###### Not Using


Alternatively, one can provide a LD reference panel (i.e., a matrix of correlations between the variants) and GWAS results (including variants allele frequencies, estimated effects, and SEs) if the sufficient statistics (**X'X** and **X'y**) or individual genotype and phenotype (**X** and **y**) for the target population are not available.

```R
str(ld)
#> num [1:1279, 1:1279] 1 0.12 0.102 0.256 0.241 ...
#> - attr(*, "dimnames")=List of 2
#>  ..$ : chr [1:1279] "wPt.0538" "wPt.8463" "wPt.6348" "wPt.9992" ...
#>  ..$ : chr [1:1279] "wPt.0538" "wPt.8463" "wPt.6348" "wPt.9992" ...

head(gwas, 3)
#>                beta       se         n          allele_freq
#> wPt.0538       0.00235    0.01269    15000      0.55476
#> wPt.8463       -0.01228   0.01206    15000      0.44582
#> wPt.6348       0.00989    0.01162    15000      0.63647
```

In this study we generated LD reference panels for African American and Hispanic ancestry populations using the All of Us (CDRv7, Controlled Tier) data, available at [Link](https://doi.org/10.5281/zenodo.16923734).

getSS() function takes as input the LD reference panel and GWAS results. The function returns the sufficient statistics (**X'X** and **X'y**).

```R
SS=getSS(ld, gwas)
str(SS)
#> List of 3
#>  $ XX:Formal class 'dgeMatrix' [package "Matrix"] with 4 slots
#>   .. ..@ Dim     : int [1:2] 1279 1279
#>   .. ..@ Dimnames:List of 2
#>   .. .. ..$ : chr [1:1279] "wPt.0538" "wPt.8463" "wPt.6348" "wPt.9992" ...
#>   .. .. ..$ : chr [1:1279] "wPt.0538" "wPt.8463" "wPt.6348" "wPt.9992" ...
#>   .. ..@ x       : num [1:1635841] 446.16 3.73 25.69 -69.47 -3.84 ...
#>   .. ..@ factors : list()
#>  $ Xy: Named num [1:1279] -14.1 -20.8 -54.7 38.9 -56.8 ...
#>   ..- attr(*, "names")= chr [1:1279] "wPt.0538" "wPt.8463" "wPt.6348" "wPt.9992" ...
#>  $ n : num 15000
```




#>  $ samplesVarE : num [1:12000] 0.799 0.667 0.664 0.672 0.721 ...
```
