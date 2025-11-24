# GPTL: Genomic Prediction Using Transfer Learning

GPTL is an R package that implements several methods for Genomic Prediction using Transfer Learning, including
-  *TL-GDES*: A Gradient Descent Algorithm with Early Stopping.
-  *TL-PR*: A Penalized regression with shrinkage towards a prior mean (e.g., estimates derived from another data set).
-  *TL-BMM*: A Bayesian model with a finite mixture prior that allows TL from multiple prior sources of information.

## Installation

- Install the package from GitHub

```R
# install.packages("remotes")
remotes::install_github("QuantGen/GPTL")
```

## Using GPTL   
   
GPTL offers three polygenic score methods using Transfer Learning. The functions `GD()`, `PR()`, and `BMM()` implement GDES, PR, and BMM, respectively.

These functions take as input the prior estimates from the source population(s) and sufficient statistics for the target population. The sufficient statistics (𝑿′𝑿, 𝑿′𝒚) can be computed from individual
genotype-phenotype data or from GWAS results and an LD reference panel.

- Download the LD reference panels and extract files:

    LD reference panels constructed using the All of Us (CDRv7, Controlled Tier) samples:
    
     [AA reference](https://zenodo.org/records/17686189/files/AA_AOU.tar.gz) (~3.5G);
     `wget https://zenodo.org/records/17686189/files/AA_AOU.tar.gz`; `tar -zxvf AA_AOU.tar.gz`
     
     [Hispanic reference](https://zenodo.org/records/17686189/files/HIS_AOU.tar.gz) (~2.7G);
     `wget https://zenodo.org/records/17686189/files/HIS_AOU.tar.gz`; `tar -zxvf HIS_AOU.tar.gz`
  
    LD reference panels constructed using the UK Biobank samples:
    
     [AFR reference](https://zenodo.org/records/17686189/files/AFR_UKB.tar.gz) (~5.1G);
     `wget https://zenodo.org/records/17686189/files/AFR_UKB.tar.gz`; `tar -zxvf AFR_UKB.tar.gz`
     
     [AMR reference](https://zenodo.org/records/17686189/files/AMR_UKB.tar.gz) (~4.3G);
     `wget https://zenodo.org/records/17686189/files/AMR_UKB.tar.gz`; `tar -zxvf AMR_UKB.tar.gz`



## Examples

- [1. GPTL using individual genotype and phenotype data](https://github.com/QuantGen/GPTL/blob/main/man/Example_Individual_Data.md)
- [2. GPTL using LD reference panel and GWAS results](https://github.com/QuantGen/GPTL/blob/main/man/Example_LD_GWAS.md)

## System Requirements

- Depends: R (>= 3.5.0)
- This package is compatible with Windows, Mac, and Linux operating systems and has been tested on Windows 7 & 10, macOS Sequoia & Sonoma, and Linux CentOS 7.

## Support

Please direct any problems or questions to Hao Wu (hwuwu95@gmail.com) and Gustavo de los Campos (gustavoc@msu.edu).
