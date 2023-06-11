\name{MOG}
\alias{dMOG}
\alias{pMOG}
\alias{qMOG}
\alias{rMOG}
\alias{mgofMOG}
\title{Marshall Olkin G Distribution}
\description{Find the maximum goodness-of-fit estimates of the distribution's parameters with cdf
\deqn{F(x) = G(x)/[a+(1-a) G(x)], a > 0,}
with any parent distribution with valid cdf G.
Also computes the pdf, cdf, quantile and generates random deviates.}
\usage{
dMOG(x,par,distr, log = FALSE)
pMOG(q,par,distr, lower.tail = TRUE, log.p = FALSE )
qMOG(p,par,distr, lower.tail = TRUE, log.p = FALSE)
rMOG(n,par,distr)
mgofMOG(data, distr, start , gofs = "CvM")}
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
\references{Marshall AW, Olkin I. A new method for adding a parameter to a family of distributions with application to the exponential and Weibull families. Biometrika 1997;84:641–652. https://doi.org/10.1093/biomet/84.3.641.}
\author{Mohamed Hussein}
\examples{
  dMOG(3,list(a = 4,rate = 0.2),"exp")
  dMOG(0.75,list(a = 0.5),"norm")
  dMOG(7,list(a = 5,rate = 0.5),"exp",log = TRUE)
  curve(dMOG(x,list(a = 2,rate = 2),"exp"),xlim = c(0,5))
  dMOG(x = seq(0,5),list(a = 3,shape = 2,scale = 3),"weibull")
  pMOG(3,list(a = 4,rate = 0.2),"exp")
  pMOG(0.75,list(a = 0.5),"norm")
  curve(pMOG(x,list(a = 0.5,mean = 2 , sd = 2),"norm"),xlim = c(-2,5))
  qMOG(0.25,list(a = 4,rate = 0.2),"exp")
  qMOG(0.75,list(a = 0.5),"norm")
  qMOG(log(0.5),list(a = 5,rate = 0.5),"exp",log.p = TRUE)
  qMOG(0.10,list(a = 1.5,rate = 0.5),"exp",lower.tail = FALSE)
  qMOG(p = c(0.25,0.5),list(a = 0.2,rate = 2),"exp")
  qMOG(p = 0.5,list(a = 3,shape = 2,scale = 3),"weibull")
  x <- rMOG(100,list(a = 3,shape = 2,scale = 3),"weibull")
  data <- rMOG(150,list(a = 2,shape = 2,scale = 3),"weibull")
  mgofMOG(x,distr="weibull", c(a=0.1,shape=0.1,scale=0.1))
  mgofMOG(data,distr="weibull", c(a=0.1,shape=0.1,scale=0.1),gofs = "AD")
}
\keyword{Marshall Olkin G distribution, family of distributions}
