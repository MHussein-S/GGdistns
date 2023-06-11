\name{PG}
\alias{dPG}
\alias{pPG}
\alias{qPG}
\alias{rPG}
\alias{mgofPG}
\title{Poisson G Distribution}
\description{Find the maximum goodness-of-fit estimates of the distribution's parameters with cdf
\deqn{F(x) = {1-exp[-aG^b(x)]}/{1-exp{-a}} , a,b > 0,}
with any parent distribution with valid cdf G.
Also computes the pdf, cdf, quantile and generates random deviates.}
\usage{
dPG(x,par,distr, log = FALSE)
pPG(q,par,distr, lower.tail = TRUE, log.p = FALSE )
qPG(p,par,distr, lower.tail = TRUE, log.p = FALSE)
rPG(n,par,distr)
mgofPG(data, distr, start , gofs = "CvM")}
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
  \item{start}{named numeric vector of the initial values of the distribution parameters
   \eqn{(a, \theta)}, \eqn{\theta} is the parameter(s) of \code{G}.}
  \item{data}{a vector of data values for which the distribution is to be fitted .}
  \item{gofs}{goodness-of-fit statistic used in etimation procedure. It can be one of \code{"CvM"}, \code{"KS"}, \code{"AD"}, \code{"ADR"}, \code{"ADL"}, \code{"AD2R"}, \code{"AD2l"} or \code{"AD2"}.  The default is \code{"CvM"}.}
}
\value{An object of the same length as \code{x} and \code{q}, giving the pdf or cdf values computed at \code{x} and \code{q} or an object of the same length as \code{p}, giving the quantile values computed at \code{p} or an object of the same length as \code{n}, giving the random deviates generated. An object giving the maximum goodness-of-fit estimates, Akaike information criterion (AC), consistent Akaike information criterion (CAC), Bayesian information criterion (BC), Hannan-Quinn information criterion (HC), the value of the goodness-of-fit statistic used and the value log-likelihood function.}
\references{Ristic MM, Nadarajah S. A new lifetime distribution. J. of Statistical Computation and Simul. 2013;84(1):135–150. http://dx.doi.org/10.1080/00949655.2012.697163.}
\author{Mohamed Hussein}
\examples{
  dPG(3,list(a = 4, b = 5 , rate = 0.2),"exp")
  dPG(0.75,list(a = 0.5 , b = 2),"norm")
  dPG(2,list(a = 5 , b = 0.5 ,rate = 0.5),"exp",log = TRUE)
  curve(dPG(x,list(a = 2 , b = 4, rate = 2),"exp"),xlim = c(0,2))
  dPG(x = seq(0,5),list(a = 3, b = 7, shape = 2,scale = 3),"weibull")
  pPG(3,list(a = 4 , b = 0.1 , rate = 0.2),"exp")
  pPG(0.75 , list(a = 0.5 , b = 0.3 , mean = 2 , sd = 2),"norm")
  qPG(0.25,list(a = 4, b = 0.1 , rate = 0.2),"exp")
  qPG(0.75,list(a = 0.5, b = 0.5),"norm")
  qPG(log(0.5),list(a = 5 , b = 2 , rate = 0.5),"exp",log.p = TRUE)
  qPG(0.10,list(a = 1.5 , b = 4.1, rate = 0.5),"exp",lower.tail = FALSE)
  qPG(p = c(0.25,0.5),list(a = 0.2 , b = 5 , rate = 2),"exp")
  qPG(p = 0.5,list(a = 3 , b = 1 , shape = 2 , scale = 3),"weibull")
  x <- rPG(100,list(a = 3 , b = 4 , rate = 3),"exp")
  mgofPG(x,distr="exp", c(a = 0.1 , b = 0.1 , rate = 0.1),"ADL")
  mgofPG(x,distr="exp", c(a = 0.1 , b = 0.1 , rate = 0.1),gofs = "AD2")
}
\keyword{Poisson G distributions, family of distributions}
