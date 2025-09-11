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
    \item{XX}{A (p x p) matrix, with row and column names being variant IDs. XX=crossprod(X), with X an genotype matrix of dimension n times p.}
    \item{Xy}{A numeric vector of length p, with names being variant IDs. Xy=crossprod(X,y).}
    \item{yy}{A numeric value. yy=crossprod(y).}
    \item{B}{(matrix) A matrix includes the estimated effects for each estimates (one row per SNP, one column per estimates), with row names being variant IDs.}
    \item{verbose}{(logical) if TRUE the history is printed, default FALSE.}
}
\value{
    Prediction correlationwill be computed for all the common variants in \code{XX}, \code{Xy}, and \code{B}.

    A list with four entries:
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
