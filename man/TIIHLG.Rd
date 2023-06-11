\name{TIIHLG}
\alias{dTIIHLG}
\alias{pTIIHLG}
\alias{qTIIHLG}
\alias{rTIIHLG}
\alias{mgofTIIHLG}
\title{Type II half-logistic G Distribution}
\description{Find the maximum goodness-of-fit estimates of the distribution's parameters with cdf
\deqn{F(x) = [2G^a(x)]/[1+G^a(x)] , a > 0,}
with any parent distribution with valid cdf G.
Also computes the pdf, cdf, quantile and generates random deviates.}
\usage{
dTIIHLG(x,par,distr, log = FALSE)
pTIIHLG(q,par,distr, lower.tail = TRUE, log.p = FALSE )
qTIIHLG(p,par,distr, lower.tail = TRUE, log.p = FALSE)
rTIIHLG(n,par,distr)
mgofTIIHLG(data, distr, start , gofs = "CvM")}
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
\references{Hassan AS, Elgarhy M, Shakil M. Type II half logistic family of distributions with applications. Pak. J. Stat. Oper. Res. 2017; 13:245–264.
https://doi.org/10.18187/pjsor.v13i2.1560.}
\author{Mohamed Hussein}
\examples{
  dTIIHLG(3,list(a = 4,rate = 0.2),"exp")
  dTIIHLG(0.75,list(a = 0.5),"norm")
  dTIIHLG(7,list(a = 5,rate = 0.5),"exp",log = TRUE)
  curve(dTIIHLG(x,list(a = 2,rate = 2),"exp"),xlim = c(0,5))
  dTIIHLG(x = seq(0,5),list(a = 3,shape = 2,scale = 3),"weibull")
  pTIIHLG(3,list(a = 4,rate = 0.2),"exp")
  pTIIHLG(0.75,list(a = 0.5),"norm")
  curve(pTIIHLG(x,list(a = 0.5,mean = 2 , sd = 2),"norm"),xlim = c(-2,5))
  qTIIHLG(0.25,list(a = 4,rate = 0.2),"exp")
  qTIIHLG(0.75,list(a = 0.5),"norm")
  qTIIHLG(log(0.5),list(a = 5,rate = 0.5),"exp",log.p = TRUE)
  qTIIHLG(0.10,list(a = 1.5,rate = 0.5),"exp",lower.tail = FALSE)
  qTIIHLG(p = c(0.25,0.5),list(a = 0.2,rate = 2),"exp")
  qTIIHLG(p = 0.5,list(a = 3,shape = 2,scale = 3),"weibull")
  x <- rTIIHLG(100,list(a = 3,shape = 2,scale = 3),"weibull")
  data <- rTIIHLG(150,list(a = 2,shape = 2,scale = 3),"weibull")
  mgofTIIHLG(x,distr="weibull", c(a=0.1,shape=0.1,scale=0.1))
  mgofTIIHLG(data,distr="weibull", c(a=0.1,shape=0.1,scale=0.1),gofs = "KS")
}
\keyword{half logistic distribution, family of distributions}
