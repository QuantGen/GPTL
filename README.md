# GPTL: Genomic Prediction Through Transfer Learning

GPTL is an R package that implements several methods for Genomic Prediction using Transfer Learning, including
-  *TL-GDES*: A Gradient Descent Algorithm with Early Stopping.
-  *TL-PR*: A Penalized regression with shrinkage towards a prior mean (e.g., estimates derived from another data set).
-  *TL-BMM*: A Bayesian model with a finite mixture prior that allows TL from multiple prior sources of information.

## Citation

Please cite [Wu et. al](https://doi.org/10.1038/s41467-026-68696-7) for the GPTL package.

## Installation

The latest version of GPTL can be installed from GitHub using:

```r
remotes::install_github("QuantGen/GPTL")
```

## Previous Versions and Reproducibility

To ensure the reproducibility of published research, previous versions of GPTL are preserved through GitHub releases and archived on Zenodo.

### Version 1.0.1 (Published Paper Version)

This version corresponds to the version of GPTL used in our published paper and is preserved to ensure reproducibility of the reported analyses and results.

- [GitHub Release (v1.0.1)](https://github.com/QuantGen/GPTL/releases/tag/v1.0.1)
- [Zenodo Archive](https://doi.org/10.5281/zenodo.24591039)

To install this specific version, use:

```r
remotes::install_github("QuantGen/GPTL@v1.0.1")
```

## Development Version

Ongoing development, including new features, improvements, bug fixes, and documentation updates, takes place in the [development branch](https://github.com/QuantGen/GPTL/tree/development).

To install the development version, use:

```r
remotes::install_github("QuantGen/GPTL@development")
```

Please note that the development version is under active development and may contain experimental features or changes that are not backward compatible with previous releases.

## Using GPTL   
   
GPTL offers three polygenic score methods using Transfer Learning (TL). The functions `GD()`, `PR()`, and `BMM()` implement Gradient Descent with Early Stopping, Penalized Regression, and Bayesian Mixture model, respectively.

These functions take as input SNP effects estimates from a source population (used as prior values to the TL algorithm) and sufficient (or summary) statistics from the target population. The sufficient statistics (𝑿′𝑿, 𝑿′𝒚) can be computed from individual genotype-phenotype data or reconstructed from GWAS results and an LD reference panel.

We provide below two links to human LD reference panels derived using All of Us and UK-Biobank, and [examples](#Examples) illustrating how to use each of the functions included in GPTL. 

## Human LD reference panels

#### LD reference panels constructed using the All of Us data (CDRv7, Controlled Tier)
    
  - **[African American](https://zenodo.org/records/17686189/files/AA_AOU.tar.gz)** (~3.5G)
    
   ```bash
     wget https://zenodo.org/records/17686189/files/AA_AOU.tar.gz
     tar -zxvf AA_AOU.tar.gz
   ```

  - **[Hispanic](https://zenodo.org/records/17686189/files/HIS_AOU.tar.gz)** (~2.7G)

   ```bash
     wget https://zenodo.org/records/17686189/files/HIS_AOU.tar.gz
     tar -zxvf HIS_AOU.tar.gz
   ```
 
#### LD reference panels constructed using the UK Biobank data
    
   - **[African](https://zenodo.org/records/17686189/files/AFR_UKB.tar.gz)** (~5.1G)

   ```bash
     wget https://zenodo.org/records/17686189/files/AFR_UKB.tar.gz
     tar -zxvf AFR_UKB.tar.gz
   ```
     
   - **[American](https://zenodo.org/records/17686189/files/AMR_UKB.tar.gz)** (~4.3G)

   ```bash
     wget https://zenodo.org/records/17686189/files/AMR_UKB.tar.gz
     tar -zxvf AMR_UKB.tar.gz
   ```

## Examples

- **[1. GPTL using individual genotype-phenotype data](https://github.com/QuantGen/GPTL/blob/main/man/Example_Individual_Data.md)**
- **[2. GPTL using LD reference panel and GWAS results](https://github.com/QuantGen/GPTL/blob/main/man/Example_LD_GWAS.md)**

## System Requirements

- Depends: R (>= 3.5.0)
- This package is compatible with Windows, Mac, and Linux operating systems and has been tested on Windows 7 & 10, macOS Sequoia & Sonoma, and Linux CentOS 7.

## Support

Please direct any problems or questions to Hao Wu (hwuwu95@gmail.com) or Gustavo de los Campos (gustavoc@msu.edu).
