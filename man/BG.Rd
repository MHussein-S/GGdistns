\name{BG}
\alias{dBG}
\alias{pBG}
\alias{qBG}
\alias{rBG}
\alias{mgofBG}
\title{Beta G Distribution}
\description{Find the maximum goodness-of-fit estimates of the distribution's parameters with cdf
\deqn{F(x) = I_G(x)(a,b), a,b > 0,}
with any parent distribution with valid cdf G.
Also computes the pdf, cdf, quantile and generates random deviates.}
\usage{
dBG(x,par,distr, log = FALSE)
pBG(q,par,distr, lower.tail = TRUE, log.p = FALSE )
qBG(p,par,distr, lower.tail = TRUE, log.p = FALSE)
rBG(n,par,distr)
mgofBG(data, distr, start , gofs = "CvM")}
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
\references{Eugene N, Lee C, Famoye F. Beta-normal distribution and its applications. Commu. in Stat.-Theory and Methods 2002; 31.  http://doi.org/10.1081/STA-120003130.}
\author{Mohamed Hussein}
\examples{
  dBG(3,list(a = 4, b = 5 , rate = 0.2),"exp")
  dBG(0.75,list(a = 0.5 , b = 2),"norm")
  dBG(7,list(a = 5 , b = 0.5 ,rate = 0.5),"exp",log = TRUE)
  curve(dBG(x,list(a = 2 , b = 0.4, rate = 2),"exp"),xlim = c(0,5))
  dBG(x = seq(0,5),list(a = 3, b = 7, shape = 2,scale = 3),"weibull")
  pBG(3,list(a = 4 , b = 0.1 , rate = 0.2),"exp")
  pBG(0.75 , list(a = 0.5 , b = 0.3 , mean = 2 , sd = 2),"norm")
  qBG(0.25,list(a = 4, b = 0.1 , rate = 0.2),"exp")
  qBG(0.75,list(a = 0.5, b = 0.5),"norm")
  qBG(log(0.5),list(a = 5 , b = 2 , rate = 0.5),"exp",log.p = TRUE)
  qBG(0.10,list(a = 1.5 , b = 4.1, rate = 0.5),"exp",lower.tail = FALSE)
  qBG(p = c(0.25,0.5),list(a = 0.2 , b = 5 , rate = 2),"exp")
  qBG(p = 0.5,list(a = 3 , b = 1 , shape = 2 , scale = 3),"weibull")
  x <- rBG(100,list(a = 0.3 , b = 4 , shape = 2 , scale = 3),"weibull")
  mgofBG(x,distr="weibull", c(a = 0.1 , b = 0.1 , shape = 0.1 , scale = 0.1))
  mgofBG(x,distr="weibull", c(a = 0.1 , b = 0.1 , shape = 0.1 , scale = 0.1),gofs = "ADR")
}
\keyword{beta generated distribution, family of distributions}
