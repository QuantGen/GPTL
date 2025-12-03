\name{Ind_DemoData}
\alias{Ind_DemoData}
\title{Demo data - individual genotype and phenotype}
\description{

Information from a collection of 599 historical CIMMYT wheat lines. The wheat data set is from 
CIMMYT's Global Wheat Program. Historically, this program has conducted numerous international 
trials across a wide variety of wheat-producing environments.

We split the data set into two subgroups (source group - 346 samples, target group - 253 samples), 
serving the illustration of the three TL-PGS methods.
  
We split the target group into (i) a training set (n=53), (ii) a calibration set (n=50), and (iii) a 
testing set (n=150). We conducted GWAS and derived LD matrix using the training set.

}

\usage{
  data(Ind_DemoData)
}

\format{

The data.frame GENO.Source contains genetypes for the source population.

The data.frame PHENO.Source contains phenotypes for the source population.

The data.frame GENO.Target contains training, calibrating and testing genetypes for the target population.

The data.frame PHENO.Target contains training, calibrating and testing phenotypes for the target population.

}

\source{
  International Maize and Wheat Improvement Center (CIMMYT), Mexico.
}

\references{
McLaren, C. G., L. Ramos, C. Lopez, and W. Eusebio. 2000. ``Applications of the geneaology manegment system.'' 
In \emph{International Crop Information System. Technical  Development Manual, version VI}, edited by McLaren, C. G., J.W. White 
and P.N. Fox. pp. 5.8-5.13. CIMMyT, Mexico: CIMMyT and IRRI. 

McLaren, C. G., R. Bruskiewich, A.M. Portugal, and A.B. Cosico. 2005. The International Rice Information System. 
A platform for meta-analysis of rice crop data. \emph{Plant Physiology} \bold{139}: 637-642.
}
\keyword{datasets}
