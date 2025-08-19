\name{getSS}
\alias{getSS}
\title{Compute sufficient statistics based on LD and GWAS}
\description{
    Compute sufficient statistics (XX and Xy) and prior (if provided), taking LD reference panel, GWAS results, and prior estimates (if provided) as inputs.
}
\usage{
getSS(ld, gwas, B=NULL, verbose=TRUE)
}
\arguments{
    \item{ld}{A (p x p) matrix, with row and column names being variant IDs. LD reference matrix.}
    \item{gwas}{A (p x 4) matrix or data.frame, with row names being variant IDs. GWAS results must consist of four columns: beta (variant effects), se (variant standard errors), n (sample sizes for GWAS), allele_freq (variant allele frequency).}
    \item{B}{(optional) A numeric vector of length p, with names being variant IDs. Alternatively, B can be a (p x k) matrix or data.frame, with row names being variant IDs, and each column being a source of prior information.}
    \item{verbose}{(logical) if TRUE the algorithm history is printed, default TRUE.}
}
\value{
    If returnPath==TRUE, a (p x nIter) matrix includes the estimated effects for each predictor (rows) for each iteration.

    If returnPath==FALSE, a numeric vector of length p includes the final estimated effects for each predictor.
}
