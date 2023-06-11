\name{EXPG}
\alias{dEXPG}
\alias{pEXPG}
\alias{qEXPG}
\alias{rEXPG}
\alias{mgofEXPG}
\title{Exponentiated G Distribution}
\description{Find the maximum goodness-of-fit estimates of the distribution's parameters with cdf
\deqn{F(x) = G(x)^a, a > 0,}
with any parent distribution with valid cdf G.
Also computes the pdf, cdf, quantile and generates random deviates.}
\usage{
dEXPG(x,par,distr, log = FALSE)
pEXPG(q,par,distr, lower.tail = TRUE, log.p = FALSE )
qEXPG(p,par,distr, lower.tail = TRUE, log.p = FALSE)
rEXPG(n,par,distr)
mgofEXPG(data, distr, start , gofs = "CvM")}
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
\references{Gupta RC, Gupta PL, Gupta RD. Modeling failure time data by Lehman alternatives. Commu. in Stat.-Theory and Methods 1998; 27:887-904. http://doi.org/10.1080/03610929808832134.}
\author{Mohamed Hussein}
\examples{
  dEXPG(3,list(a = 4,rate = 0.2),"exp")
  dEXPG(0.75,list(a = 0.5),"norm")
  dEXPG(7,list(a = 5,rate = 0.5),"exp",log = TRUE)
  curve(dEXPG(x,list(a = 2,rate = 2),"exp"),xlim = c(0,5))
  dEXPG(x = seq(0,5),list(a = 3,shape = 2,scale = 3),"weibull")
  pEXPG(3,list(a = 4,rate = 0.2),"exp")
  pEXPG(0.75,list(a = 0.5),"norm")
  curve(pEXPG(x,list(a = 0.5,mean = 2 , sd = 2),"norm"),xlim = c(-2,5))
  qEXPG(0.25,list(a = 4,rate = 0.2),"exp")
  qEXPG(0.75,list(a = 0.5),"norm")
  qEXPG(log(0.5),list(a = 5,rate = 0.5),"exp",log.p = TRUE)
  qEXPG(0.10,list(a = 1.5,rate = 0.5),"exp",lower.tail = FALSE)
  qEXPG(p = c(0.25,0.5),list(a = 0.2,rate = 2),"exp")
  qEXPG(p = 0.5,list(a = 3,shape = 2,scale = 3),"weibull")
  x <- rEXPG(100,list(a = 3,shape = 2,scale = 3),"weibull")
  data <- rEXPG(150,list(a = 2,shape = 2,scale = 3),"weibull")
  mgofEXPG(x,distr="weibull", c(a=0.1,shape=0.1,scale=0.1))
  mgofEXPG(data,distr="weibull", c(a=0.1,shape=0.1,scale=0.1),gofs = "AD")
}
\keyword{exponentiated distribution, family of distributions}
