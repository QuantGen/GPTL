\name{GD}
\alias{GD}
\title{Transfer Learning using Gradient Descent with Early Stopping}
\description{
    Performs Gradient Descent with Early Stopping.
}
\usage{
GD(XX,Xy,p=ncol(XX),b=rep(0,p),nIter=10,learning_rate=1/50,returnPath=FALSE)
}
\arguments{
    \item{XX}{A matrix, XX=crossprod(X), with X an incidence matrix of dimension n times p.}
    \item{Xy}{A numeric vector of length p, Xy=crossprod(X,y).}
    \item{p}{(integer) Number of variants. If NULL you must provide XX.}
    \item{b}{A numeric vector of length p. Prior estimates.}
    \item{nIter}{(integer) Number of iteration that the algorithm will run.}
    \item{learning_rate}{(integer) Step size of the algorithm.}
    \item{returnPath}{(logical) if TRUE the iteration history is printed, default TRUE.}
}
\value{
    A (p x nIter) matrix includes the estimated effects for each predictor (rows) for each iteration.
}
