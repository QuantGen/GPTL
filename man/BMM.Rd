\name{BMM}
\alias{BMM}
\title{Transfer Learning using Bayesian model with an informative finite mixture prior}
\description{
    Performs Bayesian regressions with an informative finite mixture prior that shrink estimates towards prior values, taking sufficient statistics as inputs.
}
\usage{
BMM(XX, Xy, B, my, vy, n, nIter=150, burnIn=50, thin=5, R2=0.25,
       nComp=matrix(ncol(B)), K=1/nComp, df0.E=5, S0.E=vy*(1-R2)*df0.E, df0.b=rep(10,nComp), 
       priorProb=rep(1/nComp,nComp), priorCounts=rep(2*nComp,nComp), fixVarE=FALSE, fixVarB=rep(FALSE,nComp), verbose=TRUE)
}
\arguments{
    \item{XX}{A (p x p) matrix, with row and column names being variant IDs. XX=crossprod(X), with X an genotype matrix of dimension n times p.}
    \item{Xy}{A numeric vector of length p, with names being variant IDs. Xy=crossprod(X,y).}
    \item{B}{(matrix) A matrix includes the estimated effects for each prior information (one row per SNP, one column per prior source of information), with row names being variant IDs.}
    \item{my}{(numeric) Mean of y.}
    \item{vy}{(numeric) Variance of y.}
    \item{n}{(integer) Sample size.}
    \item{nIter}{(integer) The number of iterations.}
    \item{burnIn}{(integer) The number of burn-in.}
    \item{thin}{(integer) The number of thinning.}
    \item{R2}{(numeric, 0<R2<1) The proportion of variance that one expects, a priori, to be explained by the regression.}
    \item{nComp}{(integer) The number of prior information in B0.}
    \item{K}{(numeric) Inverse of the number of prior information in B0.}
    \item{df0.E}{(numeric) The degree of freedom for the scaled inverse-chi squared prior assigned to the residual variance. The default value for the df parameter is 5.}
    \item{S0.E}{(numeric) The scale parameter for the scaled inverse-chi squared prior assigned to the residual variance. If the scale is not specified a value is calculated so that the prior mode of the residual variance equals var(y)*(1-R2)*df0.E.}
    \item{df0.b}{A numeric vector of length nComp, the parameter for the Dirichlet prior assigned to the prior information.}
    \item{priorProb}{A numeric vector of length nComp, the mixture prior (\eqn{\pi}).}
    \item{priorCounts}{(integer) Counts of the mixture prior (\eqn{\pi}).}
    \item{fixVarE}{(logical) if TRUE the error variance is fixed and computed depending on R2, default FALSE.}
    \item{fixVarB}{(logical) if TRUE the effects variance for each prior component is fixed and computed depending on R2, default FALSE.}
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
