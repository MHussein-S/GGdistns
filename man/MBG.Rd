\name{MBG}
\alias{dMBG}
\alias{pMBG}
\alias{qMBG}
\alias{rMBG}
\alias{mgofMBG}
\title{Modified beta G Distribution}
\description{Find the maximum goodness-of-fit estimates of the distribution's parameters with cdf
\deqn{F(x) = I_{c G(x)/[1+(c-1)G(x)]}(a,b) , a,b,c > 0,}
with any parent distribution with valid cdf G.
Also computes the pdf, cdf, quantile and generates random deviates.}
\usage{
dMBG(x,par,distr, log = FALSE)
pMBG(q,par,distr, lower.tail = TRUE, log.p = FALSE )
qMBG(p,par,distr, lower.tail = TRUE, log.p = FALSE)
rMBG(n,par,distr)
mgofMBG(data, distr, start , gofs = "CvM")}
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
\references{Nadarajah S, Teimouri M, Shih SH. Modified beta distributions. Sankhya B 2014;76(1):19–48. https://doi.org/10.1007/s13571-013-0077-0.}
\author{Mohamed Hussein}
\examples{
  dMBG(3,list(a = 4, b = 5 , c = 0.2 , rate = 0.2),"exp")
  dMBG(0.75,list(a = 0.5 , b = 2 , c = 3),"norm")
  dMBG(2,list(a = 5 , b = 0.5 , c = 0.9 ,rate = 0.5),"exp",log = TRUE)
  curve(dMBG(x,list(a = 0.2 , b = 4 , c = 5, rate = 2),"exp"),xlim = c(0,1))
  dMBG(x = seq(0,5),list(a = 0.3, b = 7 , c = 0.8, shape = 2,scale = 3),"weibull")
  pMBG(3,list(a = 4 , b = 0.1 , c = 2.5, rate = 0.2),"exp")
  pMBG(0.75 , list(a = 0.5 , b = 0.3 , c = 2.5, mean = 2 , sd = 2),"norm")
  qMBG(0.25,list(a = 4, b = 0.1 , c = 3, rate = 0.2),"exp")
  qMBG(0.75,list(a = 0.5, b = 0.5 , c = 1.8),"norm")
  qMBG(log(0.5),list(a = 5 , b = 2 , c = 4, rate = 0.5),"exp",log.p = TRUE)
  qMBG(0.10,list(a = 1.5 , b = 4.1 , c = 0.2, rate = 0.5),"exp",lower.tail = FALSE)
  qMBG(p = c(0.25,0.5),list(a = 0.2 , b = 5, c = 2 , rate = 2),"exp")
  qMBG(p = 0.5,list(a = 3 , b = 1 , c = 6, shape = 2 , scale = 3),"weibull")
  x <- rMBG(100,list(a = 3 , b = 4 , c = 5, rate = 3),"exp")
  mgofMBG(x,distr="exp", c(a = 0.1 , b = 0.1 , c = 0.1 , rate = 0.1))
  mgofMBG(x,distr="exp", c(a = 0.1 , b = 0.1 , c = 0.1 , rate = 0.1),gofs = "KS")
}
\keyword{beta generated distributions, family of distributions}
