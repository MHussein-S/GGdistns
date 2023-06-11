\name{MG}
\alias{dMG}
\alias{pMG}
\alias{qMG}
\alias{rMG}
\alias{mgofMG}
\title{Muth G Distribution}
\description{Find the maximum goodness-of-fit estimates of the distribution's parameters with cdf
\deqn{F(x) = G^{-a}(x)exp{-1/a [G^{-a}(x)-1]},a \le 1, a \ne 0,}
with any parent distribution with valid cdf G.
Also computes the pdf, cdf, approximate quantile and generates random deviates.}
\usage{
dMG(x,par,distr, log = FALSE)
pMG(q,par,distr, lower.tail = TRUE, log.p = FALSE )
qMG(p,par,distr, lower = 0, upper = 1e5, lower.tail = TRUE , log.p = FALSE)
rMG(n,par,distr)
mgofMG(data, distr, start , gofs = "CvM")}
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
\references{Almarashi AM , Elgarhy M. A new muth generated family of distributions with applications. J. of Nonlinear Sci. and Appli. 2018; 11(10):1171-1184. http://doi.org/10.22436/jnsa.011.10.06.}
\author{Mohamed Hussein}
\examples{
  dMG(3,list(a = 0.4,rate = 0.2),"exp")
  dMG(0.75,list(a = -0.5),"norm")
  dMG(7,list(a = 0.5,rate = 0.5),"exp",log = TRUE)
  curve(dMG(x,list(a = -0.2,rate = 2),"exp"),xlim = c(0,5))
  dMG(x = seq(0,5),list(a = 0.3,shape = 2,scale = 3),"weibull")
  pMG(3,list(a = 0.4,rate = 0.2),"exp")
  pMG(0.75,list(a = -0.5),"norm")
  curve(pMG(x,list(a = 0.5,mean = 2 , sd = 2),"norm"),xlim = c(-2,5))
  qMG(0.25,list(a = -0.4,rate = 0.2),"exp")
  qMG(0.75,list(a = -0.5),"norm",upper = 1e7)
  qMG(log(0.5),list(a = -0.25,rate = 0.5),"exp",log.p = TRUE)
  qMG(0.10,list(a = -0.35,rate = 0.5),"exp",lower = 0,upper = 1000,lower.tail = FALSE)
  qMG(p = c(0.25,0.5),list(a = 0.2,rate = 2),"exp")
  qMG(p = 0.5,list(a = -0.3,shape = 2,scale = 3),"weibull")
  x <- rMG(100,list(a = 0.3,shape = 2,scale = 3),"weibull")
  data <- rMG(150,list(a = -0.2,shape = 2,scale = 3),"weibull")
  mgofMG(x,distr="weibull", c(a=0.1,shape=0.1,scale=0.1))
  mgofMG(data,distr="weibull", c(a=0.1,shape=0.1,scale=0.1),gofs = "AD")
}
\keyword{Muth distribution, family of distributions}
