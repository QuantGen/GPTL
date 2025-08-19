\name{GD}
\alias{GD}
\title{Transfer Learning using Gradient Descent with Early Stopping}
\description{
    Performs Gradient Descent with Early Stopping, taking sufficient statistics as inputs.
}
\usage{
GD.SS(XX, Xy, b, nIter = 10, learning_rate = 1/50, returnPath = FALSE ,verbose=TRUE)
}
\arguments{
    \item{XX}{A (p x p) matrix, with row and column names being variant IDs. XX=crossprod(X), with X an genotype matrix of dimension n times p.}
    \item{Xy}{A numeric vector of length p, with names being variant IDs. Xy=crossprod(X,y).}
    \item{b}{A numeric vector of length p, with names being variant IDs. Alternatively, b can be a (p x 1) matrix or data.frame, with row names being variant IDs. Prior estimates.}
    \item{nIter}{(integer) Number of iterations that the algorithm will run.}
    \item{learning_rate}{(integer) Step size of the algorithm.}
    \item{returnPath}{(logical) if TRUE the iteration history is printed, default FASLE.}
    \item{verbose}{(logical) if TRUE the algorithm history is printed, default TRUE.}
}
\value{
    If returnPath==TRUE, a (p x nIter) matrix includes the estimated effects for each predictor (rows) for each iteration.

    If returnPath==FALSE, a numeric vector of length p includes the final estimated effects for each predictor.
}
