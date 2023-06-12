\name{setpar}
\alias{setpar}
\title{Split Distribution Parameter}
\description{Split the List of Parameters.}
\usage{
setpar(distname,allpar)}
\arguments{
  \item{distname}{character string specifying the pdf (or cdf) of the parent distribution}
  \item{allpar}{named list of the all distribution parameters}.
}
\value{named list of the parent distribution's parameters, named list of the new parameters.}
\author{Mohamed Hussein}
\examples{
  setpar("pbeta",allpar=list(a=5,b=5,shape1=3,shape2=2))
}
