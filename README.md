# GPTL: Genomic Prediction Using Transfer Learning

GPTL is an R package that implements several methods for Genomic Prediction using Transfer Learning, including
-  *TL-GDES*: A Gradient Descent Algorithm with Early Stopping.
-  *TL-PR*: A Penalized regression with shrinkage towards a prior mean (e.g., estimates derived from another data set).
-  *TL-BMM*: A Bayesian model with a finite mixture prior that allows TL from multiple prior sources of information.

## Getting Started

- Install the package from GitHub

```R
# install.packages("remotes")
remotes::install_github("QuantGen/GPTL")
```

- Download the LD reference panels and extract files:

    LD reference panels constructed using the All of Us Verison 7 samples:
    
     [AA reference](https://zenodo.org/records/16923735/files/LD_MAP_AA.csv) (~3.5G);
     `wget https://zenodo.org/records/16923735/files/LD_MAP_AA.csv`
     
     [Hispanic reference](https://) (~3.84G);
     `wget `
  
    LD reference panels constructed using the UK Biobank samples:
    
     [AMR reference](https://) (~3.5G);
     `wget `
     
     [AFR reference](https://) (~3.84G);
     `wget `

     [EUR reference](https://) (~3.84G);
     `wget `

## Using PROSPER   
   
efovbsjlbfvlws

## Examples

- [1. GPTL using individual genotype and phenotype data](https://github.com/QuantGen/GPTL/blob/main/man/Example_Individual_Data.md)
- [2. GPTL using LD reference panel and GWAS results](https://github.com/QuantGen/GPTL/blob/main/man/Example_LD_GWAS.md)

## System Requirements

- Depends: R (>= 3.5.0)
- This package is compatible with Windows, Mac, and Linux operating systems and has been tested on Windows 7 & 10, macOS Sequoia & Sonoma, and Linux CentOS 7.

## Support

Please direct any problems or questions to Hao Wu (hwuwu95@gmail.com) and Gustavo de los Campos (gustavoc@msu.edu).
