\name{getCor}
\alias{getCor}
\title{Compute prediction correlation based on sufficient statistics and estimated effects}
\description{
    Compute prediction correlation, taking sufficient statistics (XX, Xy, and yy) and estimated PGS effects as inputs.
}
\usage{
getCor(XX, Xy, yy, B, verbose=FALSE)
}
\arguments{
    \item{ld}{A (p x p) matrix, with row and column names being variant IDs. LD reference matrix.}
    \item{gwas}{A (p x 4) matrix or data.frame, with row names being variant IDs. GWAS results must consist of four columns: beta (variant effects), se (variant standard errors), n (sample sizes for GWAS), allele_freq (variant allele frequency).}
    \item{B}{(optional) A numeric vector of length p, with names being variant IDs. Alternatively, \code{B} can be a (p x k) matrix or data.frame, with row names being variant IDs, and each column being a source of prior information.}
    \item{verbose}{(logical) if TRUE the algorithm history is printed, default TRUE.}
}
\value{
    Sufficient statistics (XX and Xy) and prior (if provided) will be computed for all the common variants in \code{ld}, \code{gwas}, and \code{B} (if provided).

    If \code{B} is provided, a list with four entries:
    \itemize{
        \item \code{XX}: A matrix of XX.
        \item \code{Xy}: A numeric vector of Xy.
        \item \code{n}: Sample size.
        \item \code{B}: A matrix of B.
    }

    If \code{B} is not provided, a list with three entries:
    \itemize{
        \item \code{XX}: A matrix of XX.
        \item \code{Xy}: A numeric vector of Xy.
        \item \code{n}: Sample size.
    }
}
