\name{GD}
\alias{GD}
\title{Transfer Learning using Gradient Descent with Early Stopping}
\description{
    Performs Gradient Descent with Early Stopping, taking LD reference panel and GWAS results as inputs.
}
\usage{
GD(ld, gwas, b, nIter = 10, learning_rate = 1/50, returnPath = FALSE)
}
\arguments{
    \item{ld}{A (p x p) matrix, with row and column names being variant IDs. LD reference panel.}
    \item{gwas}{A data frame. Must provide GWAS results that consist of columns id (variant IDs), beta (variant effects), se (variant standard errors), n (sample sizes for GWAS), allele_freq (variant allele frequency).}
    \item{b}{A numeric column vector of length p, with row names being variant IDs. Prior estimates.}
    \item{nIter}{(integer) Number of iterations that the algorithm will run.}
    \item{learning_rate}{(integer) Step size of the algorithm.}
    \item{returnPath}{(logical) if TRUE the iteration history is printed, default TRUE.}
}
\value{
    If returnPath==TRUE, a (p x nIter) matrix includes the estimated effects for each predictor (rows) for each iteration.

    If returnPath==FALSE, a numeric vector of length p includes the final estimated effects for each predictor.
}
