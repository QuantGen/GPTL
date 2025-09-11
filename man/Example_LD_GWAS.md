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

