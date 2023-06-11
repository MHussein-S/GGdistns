\name{GBG}
\alias{dGBG}
\alias{pGBG}
\alias{qGBG}
\alias{rGBG}
\alias{mgofGBG}
\title{Generalized beta G Distribution}
\description{Find the maximum goodness-of-fit estimates of the distribution's parameters with cdf
\deqn{F(x) = I_{G^c(x)}(a,b) , a,b,c > 0,}
with any parent distribution with valid cdf G.
Also computes the pdf, cdf, quantile and generates random deviates.}
\usage{
dGBG(x,par,distr, log = FALSE)
pGBG(q,par,distr, lower.tail = TRUE, log.p = FALSE )
qGBG(p,par,distr, lower.tail = TRUE, log.p = FALSE)
rGBG(n,par,distr)
mgofGBG(data, distr, start , gofs = "CvM")}
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
\references{Alexander C, Cordeiro GM, Ortega EMM, Sarabia JM. Generalized beta-generated distributions. Comput. Stat.\&Data Analysis 2012; 56:1880–1897. https://doi.org/10.1016/j.csda.2011.11.015.}
\author{Mohamed Hussein}
\examples{
  dGBG(3,list(a = 4, b = 5 , c = 2, rate = 0.2),"exp")
  dGBG(0.75,list(a = 0.5 , b = 2 , c = 0.5),"norm")
  dGBG(2,list(a = 5 , b = 0.5 , c = 1, rate = 0.5),"exp",log = TRUE)
  curve(dGBG(x,list(a = 0.2 , b = 4 , c = 3 , rate = 2),"exp"),xlim = c(0,1))
  dGBG(x = seq(0,5),list(a = 0.3, b = 7 , c = 3, shape = 2,scale = 3),"weibull")
  pGBG(3,list(a = 4 , b = 0.1 , c = 1 , rate = 0.2),"exp")
  pGBG(0.75 , list(a = 0.5 , b = 0.3 , c = 5 , mean = 2 , sd = 2),"norm")
  qGBG(0.25,list(a = 4, b = 0.1 , c = 0.7, rate = 0.2),"exp")
  qGBG(0.75,list(a = 0.5, b = 0.5 , c = 0.8),"norm")
  qGBG(log(0.5),list(a = 5 , b = 2 , c = 6, rate = 0.5),"exp",log.p = TRUE)
  qGBG(0.10,list(a = 1.5 , b = 4.1 , c = 2 , rate = 0.5),"exp",lower.tail = FALSE)
  qGBG(p = c(0.25,0.5),list(a = 2 , b = 5 , c = 1, rate = 2),"exp")
  qGBG(p = 0.5,list(a = 3 , b = 1 , c = 0.7 , shape = 2 , scale = 3),"weibull")
  x <- rGBG(100,list(a = 3 , b = 4 , c = 3 , rate = 4),"exp")
  mgofGBG(x,distr="exp", c(a = 0.1 , b = 0.1 , c = 0.1 , rate = 0.1))
  mgofGBG(x,distr="exp", c(a = 0.1 , b = 0.1 , c = 0.1 , rate = 0.1),gofs = "KS")
}
\keyword{beta generated distributions, family of distributions}
