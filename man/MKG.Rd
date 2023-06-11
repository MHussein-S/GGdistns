\name{MKG}
\alias{dMKG}
\alias{pMKG}
\alias{qMKG}
\alias{rMKG}
\alias{mgofMKG}
\title{Modefied Kies G Distribution}
\description{Find the maximum goodness-of-fit estimates of the distribution's parameters with cdf
\deqn{F(x) = 1-exp{-[{G(x)}/{1-G(x)}]^a}, a > 0,}
with any parent distribution with valid cdf G.
Also computes the pdf, cdf, quantile and generates random deviates.}
\usage{
dMKG(x,par,distr, log = FALSE)
pMKG(q,par,distr, lower.tail = TRUE, log.p = FALSE )
qMKG(p,par,distr, lower.tail = TRUE, log.p = FALSE)
rMKG(n,par,distr)
mgofMKG(data, distr, start , gofs = "CvM")}
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
\references{Al-Babtain AA, Shakhatreh MK, Nassar M, Afify AZ. A new modified Kies family: properties, estimation under complete and type-II censored samples, and engineering applications. Mathematics 2020; 8(8):1345. https://doi.org/10.3390/math8081345.}
\author{Mohamed Hussein}
\examples{
  dMKG(3,list(a = 4,rate = 0.2),"exp")
  dMKG(0.75,list(a = 0.5),"norm")
  dMKG(7,list(a = 5,rate = 0.5),"exp",log = TRUE)
  curve(dMKG(x,list(a = 2,rate = 2),"exp"),xlim = c(0,1))
  dMKG(x = seq(0,5),list(a = 3,shape = 2,scale = 3),"weibull")
  pMKG(3,list(a = 4,rate = 0.2),"exp")
  pMKG(0.75,list(a = 0.5),"norm")
  curve(pMKG(x,list(a = 0.5,mean = 2 , sd = 2),"norm"),xlim = c(-2,10))
  qMKG(0.25,list(a = 4,rate = 0.2),"exp")
  qMKG(0.75,list(a = 0.5),"norm")
  qMKG(log(0.5),list(a = 5,rate = 0.5),"exp",log.p = TRUE)
  qMKG(0.10,list(a = 1.5,rate = 0.5),"exp",lower.tail = FALSE)
  qMKG(p = c(0.25,0.5),list(a = 0.2,rate = 2),"exp")
  qMKG(p = 0.5,list(a = 3,shape = 2,scale = 3),"weibull")
  x <- rMKG(100,list(a = 0.4,shape = 0.6,scale = 3),"weibull")
  data <- rMKG(150,list(a = 0.8,shape = 2,scale = 3),"weibull")
  mgofMKG(x,distr="weibull", c(a=0.1,shape=0.1,scale=0.1),"ADR")
  mgofMKG(data,distr="weibull", c(a=0.1,shape=0.1,scale=0.1),gofs = "AD")
}
\keyword{Kies distribution, family of distributions}
