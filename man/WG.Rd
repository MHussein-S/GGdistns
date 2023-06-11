\name{WG}
\alias{dWG}
\alias{pWG}
\alias{qWG}
\alias{rWG}
\alias{mgofWG}
\title{Weibull G Distribution}
\description{Find the maximum goodness-of-fit estimates of the distribution's parameters with cdf
\deqn{F(x) = 1-exp {-[-log[1-G(x)]/b]^a}, a,b > 0,}
with any parent distribution with valid cdf G.
Also computes the pdf, cdf, quantile and generates random deviates.}
\usage{
dWG(x,par,distr, log = FALSE)
pWG(q,par,distr, lower.tail = TRUE, log.p = FALSE )
qWG(p,par,distr, lower.tail = TRUE, log.p = FALSE)
rWG(n,par,distr)
mgofWG(data, distr, start , gofs = "CvM")}
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
\references{Alzaatreh A, Lee C, Famoye F. A new method for generating families of continuous distributions. Metron 2013; 71:63–79. http://doi.org/10.1007/s40300-013-0007-y.}
\author{Mohamed Hussein}
\examples{
  dWG(3,list(a = 4, b = 5 , rate = 0.2),"exp")
  dWG(0.75,list(a = 0.5 , b = 2),"norm")
  dWG(2,list(a = 5 , b = 0.5 ,rate = 0.5),"exp",log = TRUE)
  curve(dWG(x,list(a = 2 , b = 4, rate = 2),"exp"),xlim = c(0,5))
  dWG(x = seq(0,5),list(a = 3, b = 7, shape = 2,scale = 3),"weibull")
  pWG(3,list(a = 4 , b = 0.1 , rate = 0.2),"exp")
  pWG(0.75 , list(a = 0.5 , b = 0.3 , mean = 2 , sd = 2),"norm")
  qWG(0.25,list(a = 4, b = 0.1 , rate = 0.2),"exp")
  qWG(0.75,list(a = 0.5, b = 0.5),"norm")
  qWG(log(0.5),list(a = 5 , b = 2 , rate = 0.5),"exp",log.p = TRUE)
  qWG(0.10,list(a = 1.5 , b = 4.1, rate = 0.5),"exp",lower.tail = FALSE)
  qWG(p = c(0.25,0.5),list(a = 0.2 , b = 5 , rate = 2),"exp")
  qWG(p = 0.5,list(a = 3 , b = 1 , shape = 2 , scale = 3),"weibull")
  x <- rWG(100,list(a = 3 , b = 4 , rate = 3),"exp")
  mgofWG(x,distr="exp", c(a = 0.1 , b = 0.1 , rate = 0.1))
  mgofWG(x,distr="exp", c(a = 0.1 , b = 0.1 , rate = 0.1),gofs = "AD2L")
}
\keyword{Weibull generated distributions, family of distributions}
