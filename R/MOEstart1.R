#------- new version of MOEstart ---
#---------------1-MOEstart------------
MOEstart1<- function(x, distr)
{
pdistname <- paste("p", distr, sep = "")
if (!exists(pdistname, mode = "function")) 
  stop(paste("The ", ddistname, " function must be defined"))
n<-length(x)
alpha0<-2
#objfn<-function(obs,alpha)
#{
#n<-length(obs)
#minusloglik <- -n*log(alpha)+sum(obs)+2*sum(log(1-(1-alpha)*exp(-obs)))
#return(minusloglik)
#}
if (distr == "norm") 
{
sd0 <- sqrt((n - 1)/n) * sd(x)
mx <- mean(x)
par1<-list(mean=mx, sd=sd0)
}
else if (distr == "exp") 
{
rate0=1/mean(x)
par1<-list(rate=rate0)
}
else if (distr == "lnorm") 
{
    if (any(x <= 0)) 
   stop("values must be positive to fit a lognormal distribution")
    lx <- log(x)
    sd0 <- sqrt((n - 1)/n) * sd(lx)
    ml <- mean(lx)
    par1<- list(meanlog=ml, sdlog=sd0)
  }
else if (distr == "gamma") 
{
    if (any(x < 0)) 
      stop("values must be positive to fit an gamma  distribution")
    n <- length(x)
    m <- mean(x)
    v <- (n - 1)/n*var(x)
    par1<- list(shape=m^2/v, rate=m/v)
  }
else if (distr == "beta") 
{
    if (any(x < 0) | any(x > 1)) 
      stop("values must be in [0-1] to fit a beta distribution")
    n <- length(x)
    m <- mean(x)
    v <- (n - 1)/n*var(x)
    aux <- m*(1-m)/v - 1
    par1<- list(shape1=m*aux, shape2=(1-m)*aux)
  }
else if (distr == "weibull") 
{
    if (any(x < 0)) 
      stop("values must be positive to fit an Weibull  distribution")
    m <- mean(log(x))
    v <- var(log(x))
    shape <- 1.2/sqrt(v)
    scale <- exp(m + 0.572/shape)
    par1<- list(shape = shape, scale = scale)
  }
else if (distr == "logis") 
{
    n <- length(x)
    m <- mean(x)
    v <- (n - 1)/n*var(x)
    par1<- list(location=m, scale=sqrt(3*v)/pi)
}
else if (distr == "cauchy") 
{
    par1<- list(location=median(x), scale=IQR(x)/2)
}
else if (distr == "unif")
{
par1<- list(min=0, max=1)
}
else if (distr == "invgamma")
{
    if (any(x < 0)) 
      stop("values must be positive to fit an inverse gamma  distribution")
    #http://en.wikipedia.org/wiki/Inverse-gamma_distribution
    m1 <- mean(x)
    m2 <- mean(x^2)
    shape <- (2*m2-m1^2)/(m2-m1^2)
    scale <- m1*m2/(m2-m1^2)
    par1<- list(shape=shape, scale=scale)
  }
else if (distr == "llogis")
 {
    if (any(x < 0)) 
      stop("values must be positive to fit a log-logistic  distribution")
    p25 <- as.numeric(quantile(x, 0.25))
    p75 <- as.numeric(quantile(x, 0.75))
    shape <- 2*log(3)/(log(p75)-log(p25))
    scale <- exp(log(p75)+log(p25))/2
    par1<- list(shape=shape, scale=scale)
}
else if (distr == "invweibull")
{
    if (any(x < 0)) 
    stop("values must be positive to fit an inverse Weibull distribution")
    g <- log(log(4))/(log(log(4/3)))
    p25 <- as.numeric(quantile(x, 0.25))
    p75 <- as.numeric(quantile(x, 0.75))
    shape <- exp((g*log(p75)-log(p25))/(g-1))
    scale <-log(log(4))/(log(shape)-log(p25))
    par1<- list(shape=shape, scale=max(scale, 1e-9))
}
else if (distr == "pareto1")
 {
    if (any(x < 0)) 
    stop("values must be positive to fit a Pareto distribution")
    #http://www.math.umt.edu/gideon/pareto.pdf
    x1 <- min(x)
    m1 <- mean(x)
    n <- length(x)
    shape <- (n*m1-x1)/(n*(m1-x1))
    min <- x1*(n*shape - 1)/(n*shape)
    par1<- list(shape=shape, min=min)
  }
else if (distr == "pareto")
  {
    if (any(x < 0)) 
    stop("values must be positive to fit a Pareto distribution")
    m1 <- mean(x)
    m2 <- mean(x^2)
    scale <- (m1*m2)/(m2-2*m1^2)
    shape <- 2*(m2-m1^2)/(m2-2*m1^2)
    par1<- list(shape=shape, scale=scale)
  }
else if (distr == "lgamma")
  {
    if (any(x < 0)) 
    stop("values must be positive to fit a log-gamma distribution")
    #p228 of Klugmann and Hogg (1984)
    m1 <- mean(log(x))
    m2 <- mean(log(x)^2)
    alpha <- m1^2/(m2-m1^2)
    lambda <- m1/(m2-m1^2)
    par1<- list(shapelog=alpha, ratelog=lambda)
  }
else if (distr == "trgamma") 
{
    if (any(x < 0)) 
    stop("values must be positive to fit an trans-gamma  distribution")
    #same as gamma with shape2=tau=1
    n <- length(x)
    m <- mean(x)
    v <- (n - 1)/n*var(x)
    par1<- list(shape1=m^2/v, shape2=1, rate=m/v)
  }
else if (distr == "invtrgamma") 
{
    if (any(x < 0)) 
    stop("values must be positive to fit an inverse trans-gamma  distribution")
    #same as gamma with shape2=tau=1
    n <- length(1/x)
    m <- mean(1/x)
    v <- (n - 1)/n*var(1/x)
    par1<- list(shape1=m^2/v, shape2=1, rate=m/v)
  }
else
stop(paste0("Unknown starting values for distribution ", distr, "."))
CDF<-do.call(pdistname, c(list(median(data)), as.list(par1)))
alpha<-CDF/(1-CDF)
start <- c(tilt=alpha,par1)
return(start)
} 
