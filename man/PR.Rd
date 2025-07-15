\name{PR.SS}
\alias{PR.SS}
\title{Transfer Learning using penalized regressions}
\description{
    Performs penalized regressions that shrinkage estimates towards prior value, taking sufficient statistics as inputs.
}
\usage{
PR(XX, Xy, b, lambda = NULL, nLambda = 30, alpha = 0, conv_threshold= 1e-04, maxIter = 500, returnPath = TRUE)
}
\arguments{
    \item{XX}{A (p x p) matrix, with row and column names being variant IDs. XX=crossprod(X), with X an genotype matrix of dimension n times p.}
    \item{Xy}{A numeric column vector of length p, with row names being variant IDs. Xy=crossprod(X,y).}
    \item{b}{A numeric column vector of length p, with row names being variant IDs. Prior estimates.}
    \item{lambda}{numeric, or a numeric vector of shrinkage parameters. If NULL a grid of values will be generated based on XX and Xy.}
    \item{nLambda}{(integer) Number of shrinkage parameters to be generated. If lambda is given, nLambda will be overriden.}
    \item{alpha}{(numeric) The Elastic Net mixing parameter, with \eqn{0 \le \alpha \le 1}.}
    \item{conv_threshold}{(numeric) A tolerance parameter to control when to stop the Coordinate Descent algorithm.}
    \item{maxIter}{(integer) The maximum number of iterations for the Coordinate Descent algorithm (only used when the algorithm is not stopped by the tolerance parameter).}
    \item{returnPath}{(logical) if TRUE the iteration history is printed, default TRUE.}
}
\value{
    If returnPath==TRUE, a list with four entries:
    \itemize{
        \item \code{B}: (p x nLambda) includes the estimated effects for each predictor (rows) for each shrinkage parameter.
        \item \code{lambda}: A numeric vector of shrinkage parameters.
        \item \code{alpha}: (numeric) Elastic Net mixing parameter.
        \item \code{conv_iter}: A numeric vector of numbers of iterations where the algorithm converged for each shrinkage parameter.
    }
    If returnPath==FALSE, a list with four entries:
    \itemize{
        \item \code{B}: (p x maxIter x nLambda) includes the estimated effects for each predictor for each shrinkage parameter in each iteration.
        \item \code{lambda}: A numeric vector of shrinkage parameters.
        \item \code{alpha}: (numeric) Elastic Net mixing parameter.
        \item \code{conv_iter}: A numeric vector of numbers of iterations where the algorithm converged for each shrinkage parameter.
    }
}
