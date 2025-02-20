\name{BMM}
\alias{BMM}
\title{Transfer Learning using Bayesian model with an informative finite mixture prior}
\description{
    Performs Bayesian regressions with an informative finite mixture prior that shrinkage estimates towards prior values.
}
\usage{
BMM(C,rhs,my,vy,n,B0=matrix(nrow=ncol(C),ncol=1,0),nIter=150,burnIn=50,thin=5,R2=0.25,
    nComp=matrix(ncol(B0)),K=1/nComp, df0.E=5,S0.E=vy*(1-R2)*df0.E,df0.b=rep(10,nComp), 
    priorProb=rep(1/nComp,nComp),priorCounts=rep(2*nComp,nComp),verbose=TRUE)
}
\arguments{
    \item{C}{A matrix, XX=crossprod(X), with X an incidence matrix of dimension n times p.}
    \item{rhs}{A numeric vector of length p, Xy=crossprod(X,y).}
    \item{my}{(numeric) Mean of y.}
    \item{vy}{(numeric) Variance of y.}
    \item{B0}{(matrix) A matrix includes the estimated effects for each prior information (one row per SNP, one column per prior source of information).}
    \item{nIter}{(integer) The number of iterations.}
    \item{burnIn}{(integer) The number of burn-in.}
    \item{thin}{(integer) The number of thinning.}
    \item{R2}{(numeric, 0<R2<1) }
    \item{nComp}{}
    \item{K}{}
    \item{df0.E}{}
    \item{S0.E}{}
    \item{df0.b}{}
    \item{priorProb}{}
    \item{priorCounts}{}
    \item{verbose}{(logical) if TRUE the iteration history is printed, default TRUE.}
}
\value{
    A list with seven entries:
    \itemize{
        \item \code{b}: A numeric vector of the estimated effects for each predictor.
        \item \code{POST.PROB}: A matrix of posterior probalities for SNP effects.
        \item \code{postMeanVarB}: A numeric vector of posterior variances for SNP effects.
        \item \code{postProb}: A numeric vector of mean posterior probalities for SNP effects.
        \item \code{samplesVarB}: A matrix of sampling variance of SNP effects.
        \item \code{samplesB}: A matrix of sampling SNP effects.
        \item \code{samplesVarE}: A matrix of variance of sampling errors.
    }
}
