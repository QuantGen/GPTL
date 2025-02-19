\name{PR}
\alias{PR}
\title{Transfer Learning using penalized regressions}
\description{
    Performs penalized regressions that shrinkage estimates towards prior value.
}
\usage{
PR(XX, Xy, p=ncol(XX), b=rep(0,p),lambda=NULL,nLambda=30,alpha=0,
    conv_threshold=1e-4,maxIter=500,returnPath=TRUE)
}
\arguments{
    \item{XX}{A matrix, XX=crossprod(X), with X an incidence matrix of dimension n times p.}
    \item{Xy}{A numeric vector of length p, Xy=crossprod(X,y).}
    \item{p}{(integer) Number of variants. If NULL you must provide XX.}
    \item{b}{A numeric vector of length p. Prior estimates.}
    \item{lambda}{numeric, or a numeric vector of shrinkage parameters. If NULL a grid of values will be generated based on XX and Xy.}
    \item{nLambda}{(integer) Number of shrinkage parameters to be generated. If lambda is given, nLambda will be overriden.}
    \item{alpha}{(numeric) Sparase parameter.}
    \item{conv_threshold}{(numeric) A tolerance parameter to control when to stop the Coordinate Descent algorithm.}
    \item{maxIter}{(integer) The maximum number of iterations for the Coordinate Descent algorithm (only used when the algorithm is not stopped by the tolerance parameter).}
    \item{returnPath}{(logical) if TRUE the iteration history is printed, default TRUE.}
}
\value{
    A list with four entries:
    \itemize{
        \item \code{B}: (p x nLambda) includes the estimated effects for each predictor (rows) for each shrinkage parameter.
        \item \code{lambda}: A numeric vector of shrinkage parameters.
        \item \code{alpha}: (numeric) Sparase parameter.
        \item \code{conv_iter}: A numeric vector of numbers of iterations where the algorithm converged for each shrinkage parameter.
    }
}
