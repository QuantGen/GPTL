### Simulating genotype and phenotype data as examples for GPTL

The following script generates demo datasets to illustrate how GPTL software works when one has access to individual genetype and phenotype data or using an LD reference panel (internal or external) and GWAS summary statistics. The generated data sets are stored as **.RData** formats:

- `Ind_DemoData.RData`, including individual genotype and phenotype data for a source and a target population, with the target population data being subset into training, calibrating, and testing sets.

- `Sum_DemoData.RData`, including prior effects estimated from a source population, and LD reference panel, GWAS summary statistics, and individual calibrating/testing data from a target population.

Here we use the [wheat](https://doi.org/10.1104/pp.105.063438) data set collected from CIMMYT's Global Wheat Program, including 599 wheat lines genotype (1279 variants) and phenotype (average grain yield).

```R
library(BGLR)
data(wheat)
```

We first reorder the variants according to hierarchical clustering.

```R
library(corrplot)
Z=scale(wheat.X,center=TRUE,scale=TRUE)/sqrt(nrow(wheat.X)-1)
COR=crossprod(Z)/(nrow(Z)-1)
colOrder= corrMatOrder(abs(COR), order = 'hclust')
X=wheat.X[,colOrder]
```

This data set has two clear clusters, we use this to illustrate how to transfer learning from one source population to improve prediction accuracy in another target population.

```R
CLUSTER=kmeans(scale(X, center=TRUE, scale=FALSE),centers=2,nstart=100)
table(CLUSTER$cluster)
#>   1   2 
#> 346 253 
```

We use samples in cluster 1 as the source data set (where information is transferred) and samples in cluster 2 as the target data set (where the PGS will be used). We also remove the monomorphic variants.

```R
GENO.Source=scale(wheat.X[CLUSTER$cluster == 1,], center=TRUE, scale=FALSE)
GENO.Target=scale(wheat.X[CLUSTER$cluster == 2,], center=TRUE, scale=FALSE)

monomorphic=which(matrixStats::colVars(GENO.Source)==0 | matrixStats::colVars(GENO.Target)==0)
GENO.Source=GENO.Source[,-monomorphic]
GENO.Target=GENO.Target[,-monomorphic]

PHENO.Source=wheat.Y[CLUSTER$cluster == 1,1]
PHENO.Target=wheat.Y[CLUSTER$cluster == 2,1]
```

We further split the target data set into (i) a training set (40%), (ii) a calibrating set (30%), and (iii) a testing set (30%).

```R
set.seed(1234)
sets <- cut(runif(nrow(GENO.Target)), breaks = c(0, 0.4, 0.7, 1), labels = c("trn","cal","tst"))
PHENO.Target=data.frame(y=PHENO.Target, sets=sets)
table(PHENO.Target$sets)
#> trn cal tst 
#> 102  80  71
```

`GENO.Source` and `PHENO.Source` consist of genotype and phenotype data for the source population (346 samples, 1270 variants). `GENO.Target` and `PHENO.Target` consist of genotype and phenotype data for the target population (253 samples, 1270 variants), with samples splitted into 3 sets, marking in `PHENO.Target`.

The above-generated data sets are saved in `Ind_DemoData.RData`, and can be loaded by `data(Ind_DemoData)` with the `GPTL` package. Example pipeline using this demo data can be found [here](https://github.com/QuantGen/GPTL/blob/main/man/Example_Individual_Data.md)

```R
save(GENO.Source, PHENO.Source, GENO.Target, PHENO.Target, file='Ind_DemoData.RData')
```

We continue to prepare demo data 2. We estimate prior effects from the source population using a Bayesian shrinkage estimation method (a Bayesian model with a Gaussian prior centered at zero, model ‘BRR’ in the **BGLR** R-package). Alternatively, if only sufficient statistics (**X'X** and **X'y**) for the source data set are provided, one can use *BLRCross()* function in the **BGLR** R-package.

```R
ETA=list(list(X=GENO.Source, model="BRR"))
fm=BGLR(y=PHENO.Source, ETA = ETA, response_type = "gaussian", nIter = 12000, burnIn = 2000, verbose = FALSE)
PRIOR=fm$ETA[[1]]$b
```

We define LD block boundaries and generate the sparse LD reference matrix for the target population.

```R
library(Matrix)
source('https://raw.githubusercontent.com/QuantGen/GPTL/refs/heads/main/misc/getBlock.R')

Z=scale(GENO.Target,center=TRUE,scale=TRUE)
COR=crossprod(Z)/(nrow(Z)-1)
R2=COR^2
blocks<-findBlocks(R2,threshold=0.1,maxGap=3)
blocks<-mergeBlocks(R2,blocks,minSize = 3)

LD=sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(ncol(GENO.Target), ncol(GENO.Target)))
for (blk in 1: max(blocks)) {
    index=which(blocks == blk)
    LD[index,index]=R2[index,index]
}
colnames(LD)=colnames(R2)
rownames(LD)=rownames(R2)
```

We run GWAS for the target population training set.

```R
library(BGData)
SMR=GWAS(y~1,data=BGData(geno=GENO.Target[PHENO.Target$sets=='trn',],pheno=PHENO.Target[PHENO.Target$sets=='trn',]),method='rayOLS')
GWAS=as.data.frame(SMR[,c(1,5,6)])
colnames(GWAS)=c('beta', 'n', 'allele_freq')
```

We retain the calibrating and testing sets for demo data 2.

```R
GENO.Target=GENO.Target[PHENO.Target$sets!='trn',]
PHENO.Target=PHENO.Target[PHENO.Target$sets!='trn',]
```

`PRIOR`, `LD`, and `GWAS` consist of prior estimates from the source population, LD reference panel and GWAS summary statistics from the target population. `GENO.Target` and `PHENO.Target` consist of genotype and phenotype data for the target population (151 samples, 1270 variants), with samples splitted into calibrating and testing sets, marking in `PHENO.Target`.

The above-generated data sets are saved in `Sum_DemoData.RData`, and can be loaded by `data(Sum_DemoData)` with the `GPTL` package. Example pipeline using this demo data can be found [here](https://github.com/QuantGen/GPTL/blob/main/man/Example_LD_GWAS.md)

```R
save(PRIOR, LD, GWAS, GENO.Target, PHENO.Target, file='Sum_DemoData.RData')
```

[Back to Homepage](https://github.com/QuantGen/GPTL/blob/main/README.md)
