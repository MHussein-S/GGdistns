\name{NMEPAG}
\alias{dNMEPAG}
\alias{pNMEPAG}
\alias{qNMEPAG}
\alias{rNMEPAG}
\alias{mgofNMEPAG}
\title{New Modified Exponent Power Alpha G Distribution}
\description{Find the maximum goodness-of-fit estimates of the distribution's parameters with cdf
\deqn{F(x) = G(x) exp^[(1-a)^2-(1-a G(x))^2], a > 0 }
with any parent distribution with valid cdf G.
Also computes the pdf, cdf, approximate quantile and generates random deviates.}
\usage{
dNMEPAG(x,par,distr, log = FALSE)
pNMEPAG(q,par,distr, lower.tail = TRUE, log.p = FALSE )
qNMEPAG(p,par,distr, lower = 0, upper = 1e5, lower.tail = TRUE , log.p = FALSE)
rNMEPAG(n,par,distr)
mgofNMEPAG(data, distr, start , gofs = "CvM")}
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
\references{Shah Z, Khan DM, Khan Z, Shafiq M, Choi J. A new modified exponent power alpha family of distributions with applications in Reliability Engineering. Processes 2022, 10(11):2250. https://doi.org/10.3390/pr10112250.}
\author{Mohamed Hussein}
\examples{
  dNMEPAG(3,list(a = 1.5 , rate = 0.2),"exp")
  dNMEPAG(0.75,list(a = 2.5 ),"norm")
  dNMEPAG(7,list(a = 0.5 , rate = 0.5),"exp",log = TRUE)
  curve(dNMEPAG(x,list(a = 1.2 , rate = 2),"exp"),xlim = c(0,5))
  dNMEPAG(x = seq(0,5),list(a = 3 , shape = 2,scale = 3),"weibull")
  pNMEPAG(3,list(a = 2.4 , rate = 0.2),"exp")
  pNMEPAG(0.75,list(a = 5),"norm")
  curve(pNMEPAG(x,list(a = 5 , mean = 2 , sd = 2),"norm"),xlim = c(0,10))
  qNMEPAG(0.25,list(a = 1.5 , rate = 0.2),"exp")
  qNMEPAG(0.75,list(a = 5),"norm")
  qNMEPAG(log(0.5),list(a = 3, rate = 0.5),"exp",log.p = TRUE)
  qNMEPAG(0.10,list(a = 15,rate = 0.5),"exp",lower = 0,upper = 1000,lower.tail = FALSE)
  qNMEPAG(p=c(0.25,0.5),list(a = 1.2 , rate=2),"exp")
  qNMEPAG(p = 0.5,list(a = 3,shape = 2,scale = 3),"weibull")
  x <- rNMEPAG(100,list(a = 3 , shape = 2 , scale = 3),"weibull")
  data <- rNMEPAG(150,list(a = 0.3 , shape = 2,scale = 0.3),"weibull")
  mgofNMEPAG(x,distr="weibull", c(a = 0.1 , shape = 0.1 , scale = 0.1),"KS")
  mgofNMEPAG(data,distr="weibull", c(a = 0.1 , shape = 0.1 , scale = 0.1),gofs = "AD")
}
\keyword{power distribution, family of distributions}
