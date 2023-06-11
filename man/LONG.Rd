\name{LONG}
\alias{dLONG}
\alias{pLONG}
\alias{qLONG}
\alias{rLONG}
\alias{mgofLONG}
\title{Log-odd normal G Distribution}
\description{Find the maximum goodness-of-fit estimates of the distribution's parameters with cdf
\deqn{F(x) = 1/2 {1+erf [ln G(x)/(1-G(x))/(a\sqrt 2)} , a > 0,}
with any parent distribution with valid cdf G.
Also computes the pdf, cdf, quantile and generates random deviates.}
\usage{
dLONG(x,par,distr, log = FALSE)
pLONG(q,par,distr, lower.tail = TRUE, log.p = FALSE )
qLONG(p,par,distr, lower.tail = TRUE, log.p = FALSE)
rLONG(n,par,distr)
mgofLONG(data, distr, start , gofs = "CvM")}
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
\references{Zubair M, Pogany TK., Cordeiro GM, Tahir MH. The log-odd normal generalized family of distributions with application.
Anais da Academia Brasileira de Cienc 2019;91(02). https://doi.org/10.1590/0001-3765201920180207.}
\author{Mohamed Hussein}
\examples{
  dLONG(3,list(a = 4,rate = 0.2),"exp")
  dLONG(0.75,list(a = 0.5),"norm")
  dLONG(7,list(a = 5,rate = 0.5),"exp",log = TRUE)
  curve(dLONG(x,list(a = 2,rate = 2),"exp"),xlim = c(0,5))
  dLONG(x = seq(0,5),list(a = 3,shape = 2,scale = 3),"weibull")
  pLONG(3,list(a = 4,rate = 0.2),"exp")
  pLONG(0.75,list(a = 0.5),"norm")
  curve(pLONG(x,list(a = 0.5,mean = 2 , sd = 2),"norm"),xlim = c(-2,5))
  qLONG(0.25,list(a = 4,rate = 0.2),"exp")
  qLONG(0.75,list(a = 0.5),"norm")
  qLONG(log(0.5),list(a = 5,rate = 0.5),"exp",log.p = TRUE)
  qLONG(0.10,list(a = 1.5,rate = 0.5),"exp",lower.tail = FALSE)
  qLONG(p = c(0.25,0.5),list(a = 0.2,rate = 2),"exp")
  qLONG(p = 0.5,list(a = 3,shape = 2,scale = 3),"weibull")
  x <- rLONG(100,list(a = 3,shape = 2,scale = 3),"weibull")
  data <- rLONG(150,list(a = 2,shape = 2,scale = 3),"weibull")
  mgofLONG(x,distr="weibull", c(a=0.1,shape=0.1,scale=0.1))
  mgofLONG(data,distr="weibull", c(a=0.4,shape=0.1,scale=0.1),gofs = "KS")
}
\keyword{log-odd distribution, family of distributions}
