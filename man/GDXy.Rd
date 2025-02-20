\name{GDXy}
\alias{GDXy}
\title{Transfer Learning using Gradient Descent with Early Stopping (X and y)}
\description{
    Performs Gradient Descent with Early Stopping that takes X and y, instead of XX and Xy.
}
\usage{
GDXy(X,y,centerX=TRUE,scaleX=FALSE,sortGradient=FALSE,...)
}
\arguments{
    \item{X}{A matrix, genotype, with X an incidence matrix of dimension n times p.}
    \item{y}{A numeric vector of length p, phenotype.}
    \item{p}{(integer) Number of variants. If NULL you must provide XX.}
    \item{b}{A numeric vector of length p. Prior estimates.}
    \item{nIter}{(integer) Number of iteration that the algorithm will run.}
    \item{learning_rate}{(integer) Step size of the algorithm.}
    \item{returnPath}{(logical) if TRUE the iteration history is printed, default TRUE.}
}
\value{
    A (p x nIter) matrix includes the estimated effects for each predictor (rows) for each iteration.
}
