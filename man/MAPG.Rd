\name{MAPG}
\alias{dMAPG}
\alias{pMAPG}
\alias{qMAPG}
\alias{rMAPG}
\alias{mgofMAPG}
\title{Modified alpha-power G Distribution}
\description{Find the maximum goodness-of-fit estimates of the distribution's parameters with cdf
\deqn{F(x) = {a^G(x) b^{G^2(x)}-1}/(ab-1), a,b \ge 1, ab \ne 1,}
with any parent distribution with valid cdf G.
Also computes the pdf, cdf, approximate quantile and generates random deviates.}
\usage{
dMAPG(x,par,distr, log = FALSE)
pMAPG(q,par,distr, lower.tail = TRUE, log.p = FALSE )
qMAPG(p,par,distr, lower = 0, upper = 1e5, lower.tail = TRUE , log.p = FALSE)
rMAPG(n,par,distr)
mgofMAPG(data, distr, start , gofs = "CvM")}
\arguments{
  \item{x}{scaler or vector of values at which the pdf will be computed}
  \item{q}{scaler or vector of values at which the cdf will be computed}
  \item{p}{scalar or vector of probabilities.}
  \item{n}{number of random deviates to be generated}
  \item{par}{named list of the distribution parameters.}
  \item{distr}{character string specifying the parent distribution.}
  \item{log}{logical; if TRUE, log(pdf) is returned}
  \item{lower.tail}{logical; if TRUE, probability are P(X <= x), otherwise, P(X > x).}
  \item{log.p}{logical; if TRUE, probabilities p are given as log(p).}
   \item{lower}{the lower end point of the interval to be searched.}
  \item{upper}{the upper end point of the interval to be searched.}
  \item{start}{named numeric vector of the initial values of the distribution parameters
   \eqn{(a, \theta)}, \eqn{\theta} is the parameter(s) of \code{G}.}
  \item{data}{a vector of data values for which the distribution is to be fitted .}
  \item{gofs}{goodness-of-fit statistic used in etimation procedure. It can be one of \code{"CvM"}, \code{"KS"}, \code{"AD"}, \code{"ADR"}, \code{"ADL"}, \code{"AD2R"}, \code{"AD2l"} or \code{"AD2"}.  The default is \code{"CvM"}.}
}
\value{An object of the same length as \code{x} and \code{q}, giving the pdf or cdf values computed at \code{x} and \code{q} or an object of the same length as \code{p}, giving the quantile values computed at \code{p} or an object of the same length as \code{n}, giving the random deviates generated. An object giving the maximum goodness-of-fit estimates, Akaike information criterion (AC), consistent Akaike information criterion (CAC), Bayesian information criterion (BC), Hannan-Quinn information criterion (HC), the value of the goodness-of-fit statistic used and the value log-likelihood function.}
\references{Hussein M, Elsayed H, Cordeiro GM. A new family of continuous distributions: properties and estimation. Symmetry 2022; 14(2):276. https://doi.org/10.3390/sym14020276.}
\author{Mohamed Hussein}
\examples{
  dMAPG(3,list(a = 1.4 , b = 2 , rate = 0.2),"exp")
  dMAPG(0.75,list(a = 2.5 , b = 3),"norm")
  dMAPG(7,list(a = 5 , b = 2 , rate = 0.5),"exp",log = TRUE)
  curve(dMAPG(x,list(a = 1.2 , b = 4 , rate = 2),"exp"),xlim = c(0,5))
  dMAPG(x = seq(0,5),list(a = 3 , b = 2.5 , shape = 2,scale = 3),"weibull")
  pMAPG(3,list(a = 2.4 , b = 3 , rate = 0.2),"exp")
  pMAPG(0.75,list(a = 5, b = 2),"norm")
  curve(pMAPG(x,list(a = 5 , b = 2 , mean = 2 , sd = 2),"norm"),xlim = c(-2,5))
  qMAPG(0.25,list(a = 3.4 , b = 4 , rate = 0.2),"exp")
  qMAPG(0.75,list(a = 3.4 , b = 4),"norm",upper = 1e7)
  qMAPG(0.10,list(a = 3.5 , b = 2 , rate = 0.5),"exp",lower = 0,upper = 1000,lower.tail = FALSE)
  qMAPG(p = c(0.25,0.5),list(a = 2 , b = 3 , rate = 2),"exp")
  qMAPG(p = 0.5,list(a = 2 , b = 3 , shape = 2 , scale = 3),"weibull")
  x <- rMAPG(100,list(a = 3 , b = 2 , shape = 2 , scale = 3),"weibull")
  data <- rMAPG(150,list(a = 5, b = 4, shape = 2,scale = 0.3),"weibull")
  mgofMAPG(x,distr="weibull", c(a=1.1,b=1.1,shape=0.1,scale=0.1))
  mgofMAPG(data,distr="weibull", c(a=1.1,b=1.1,shape=0.1,scale=0.1),gofs = "AD2R")
}
\keyword{Alpha Power distribution, family of distributions}
