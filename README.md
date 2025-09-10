GPTL: Genomic Prediction Using Transfer Learning
================================================

GPTL is an R package that implements several methods for Genomic Prediction using Transfer Learning, including
 -  *TL-GDES*: A Gradient Descent Algorithm with Early Stopping.
 -  *TL-PR*: A Penalized regression with srhinkage towards a prior mean (e.g., estimates derived from another data set).
 -  *TL-BMM*: A Bayesian model with a finite mixture prior that allows TL from multiple prior sources of information.

Installation
------------

To install from Github:

```R
# install.packages("remotes")
remotes::install_github("QuantGen/GPTL")
```

This should take less than 30 seconds for a typical computer when including installation of all non-base dependencies.

This package should be compatible with Windows, Mac, and Linux operating systems and has been tested on Windows 7 & 10, macOS Sequoia & Sonoma, and Linux CentOS 7.

System Requirements
-------------------

Depends: R (>= 3.5.0)

Examples
--------

 - [1. GPTL using individual genotype and phenotype data](https://github.com/QuantGen/GPTL/blob/main/man/Example_Individual_Data.Rd)
  
 - [2. GPTL using LD reference panel and GWAS results](https://github.com/QuantGen/GPTL/blob/main/man/Example_LD_GWAS.Rd)
