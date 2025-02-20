\name{GDXy}
\alias{GDXy}
\title{Transfer Learning using Gradient Descent with Early Stopping (X and y)}
\description{
    Performs Gradient Descent with Early Stopping that takes X and y, instead of XX and Xy.
}
\usage{
GDXy(X,y,centerX=TRUE,scaleX=FALSE,...)
}
\arguments{
    \item{X}{A matrix, genotype, with X an incidence matrix of dimension n times p.}
    \item{y}{A numeric vector of length p, phenotype.}
    \item{centerX}{(logical) Whether centers the columns of X, default TRUE.}
    \item{scaleX}{(logical) Whether scales the columns of X, default FALSE.}
    \item{...}{Other arguments passed on to \code{\link[=GD]{GD()}}'s params argument.}
}
\value{
    A (p x nIter) matrix includes the estimated effects for each predictor (rows) for each iteration.
}
