\name{MPoG}
\alias{dMPoG}
\alias{pMPoG}
\alias{qMPoG}
\alias{rMPoG}
\alias{mgofMPoG}
\title{Modified power G Distribution}
\description{Find the maximum goodness-of-fit estimates of the distribution's parameters with cdf
\deqn{F(x) = a^{G(x)-1} G(x), a \ge exp(-1) }
with any parent distribution with valid cdf G.
Also computes the pdf, cdf, approximate quantile and generates random deviates.}
\usage{
dMPoG(x,par,distr, log = FALSE)
pMPoG(q,par,distr, lower.tail = TRUE, log.p = FALSE )
qMPoG(p,par,distr, lower = 0, upper = 1e5, lower.tail = TRUE , log.p = FALSE)
rMPoG(n,par,distr)
mgofMPoG(data, distr, start , gofs = "CvM")}
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
\value{An object of the same length as \code{x} and \code{q}, giving the pdf or cdf values computed at \code{x} and \code{q} or an object of the same length as \code{p}, giving an approximate quantile values computed at \code{p} or an object of the same length as \code{n}, giving the random deviates generated. An object giving the maximum goodness-of-fit estimates, Akaike information criterion (AC), consistent Akaike information criterion (CAC), Bayesian information criterion (BC), Hannan-Quinn information criterion (HC), the value of the goodness-of-fit statistic used and the value log-likelihood function.}
\references{Hussein M, Cordeiro GM. A modified power family of distributions: properties, simulations and applications. Mathematics 2022; 10(7):1035. https://doi.org/10.3390/math10071035.}
\author{Mohamed Hussein}
\examples{
  dMPoG(3,list(a = 1.5 , rate = 0.2),"exp")
  dMPoG(0.75,list(a = 2.5 ),"norm")
  dMPoG(7,list(a = 0.5 , rate = 0.5),"exp",log = TRUE)
  curve(dMPoG(x,list(a = 1.2 , rate = 2),"exp"),xlim = c(0,5))
  dMPoG(x = seq(0,5),list(a = 3 , shape = 2,scale = 3),"weibull")
  pMPoG(3,list(a = 2.4 , rate = 0.2),"exp")
  pMPoG(0.75,list(a = 5),"norm")
  curve(pMPoG(x,list(a = 5 , mean = 2 , sd = 2),"norm"),xlim = c(-2,5))
  qMPoG(0.25,list(a = 1.5 , rate = 0.2),"exp")
  qMPoG(0.75,list(a = 5),"norm")
  qMPoG(log(0.5),list(a = 3, rate = 0.5),"exp",log.p = TRUE)
  qMPoG(0.10,list(a = 15,rate = 0.5),"exp",lower = 0,upper = 1000,lower.tail = FALSE)
  qMPoG(p=c(0.25,0.5),list(a = 1.2 , rate=2),"exp")
  qMPoG(p = 0.5,list(a = 3,shape = 2,scale = 3),"weibull")
  x <- rMPoG(100,list(a = 3 , shape = 2 , scale = 3),"weibull")
  data <- rMPoG(150,list(a = 5 , shape = 2,scale = 0.3),"weibull")
  mgofMPoG(x,distr="weibull", c(a = 0.4 , shape = 0.1 , scale = 0.1))
  mgofMPoG(data,distr="weibull", c(a = 1 , shape = 0.1 , scale = 0.1),gofs = "AD2R")
}
\keyword{Power distributions, family of distributions}
