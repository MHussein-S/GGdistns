#Maximum goodness-of-fit estimates of Marshall-Olkin family
MOEmgof <- function (data, distr, start,gof = "CvM",...)
{
  if (!(is.numeric(data) & length(data)>1))
    stop("data must be a numeric vector of length greater than 1")
  if(!is.list(start))
    start<- as.list(start)
  if (is.null(names(start)))
    stop("'start' must be a named list or vector")
  if (!is.character(distr))
    stop("distr must be a character string naming a distribution")
  #------
  tiltpar<-match("tilt",names(start))
  if (is.na(tiltpar))
    stop("Starting value of \"Tilt\" is not defined")
  alpha<-start$tilt
  par<-start[-tiltpar]
  #--------
  pdistname <- paste("p",distr,sep="")
  ddistname <- paste("d",distr,sep="")
  if (!exists(ddistname, mode="function"))
    stop(paste("The ", distr, " distribution is not defined"))
  parnm <- names(par)
  args <- names(formals(ddistname))
  m <- match(parnm,args)
  if (any(is.na(m)))
    stop("'start' specifies names which are not arguments to ",distr)
  #--------
  data<-sort(data)

  gof <- match.arg(gof, c("CvM", "KS", "AD", "ADR", "ADL", "AD2R", "AD2L", "AD2"))

  if (gof == "CvM")
    fnobj <- function(allpar,pdistnam,data)
    {
      if(!is.list(allpar))
        allpar<- as.list(allpar)
      tiltpar<-match("tilt",names(allpar))
      alpha<-allpar$tilt
      distpar<-allpar[-tiltpar]
      alpha_bar<-1-alpha
      CDF<-do.call(pdistnam,c(list(data),as.list(distpar)))
      theop<-CDF/(1-alpha_bar*(1-CDF))
      n<-length(theop)
      1/(12*n) + sum( ( theop - (2 * 1:n - 1)/(2 * n) )^2 )
      }
  else
    if (gof == "KS")
      fnobj <- function(allpar,pdistnam,data)
      {
        if(!is.list(allpar))
          allpar<- as.list(allpar)
        tiltpar<-match("tilt",names(allpar))
        alpha<-allpar$tilt
        distpar<-allpar[-tiltpar]
        alpha_bar<-1-alpha
        CDF<-do.call(pdistnam,c(list(data),as.list(distpar)))
        theop<-CDF/(1-alpha_bar*(1-CDF))
        n<-length(theop)
        obspu <- seq(1,n)/n
        obspl <- seq(0,n-1)/n
        max(pmax(abs(theop-obspu),abs(theop-obspl)))
        }
  else
    if (gof == "AD")
      fnobj <- function(allpar,pdistnam,data)
      {
        if(!is.list(allpar))
          allpar<- as.list(allpar)
        tiltpar<-match("tilt",names(allpar))
        alpha<-allpar$tilt
        distpar<-allpar[-tiltpar]
        alpha_bar<-1-alpha
        CDF<-do.call(pdistnam,c(list(data),as.list(distpar)))
        theop<-CDF/(1-alpha_bar*(1-CDF))
        n <- length(theop)
        - n - mean( (2 * 1:n - 1) * (log(theop) + log(1 - rev(theop))))
        }
  else
    if (gof == "ADR")
      fnobj <- function(allpar,pdistnam,data)
      {
        if(!is.list(allpar))
          allpar<- as.list(allpar)
        tiltpar<-match("tilt",names(allpar))
        alpha<-allpar$tilt
        distpar<-allpar[-tiltpar]
        alpha_bar<-1-alpha
        CDF<-do.call(pdistnam,c(list(data),as.list(distpar)))
        theop<-CDF/(1-alpha_bar*(1-CDF))
        n <- length(theop)
        n/2 - 2 * sum(theop) - mean ( (2 * 1:n - 1) * log(1 - rev(theop)))
        }
  else
    if (gof == "ADL")
      fnobj <- function(allpar,pdistnam,data)
      {
        if(!is.list(allpar))
          allpar<- as.list(allpar)
        tiltpar<-match("tilt",names(allpar))
        alpha<-allpar$tilt
        distpar<-allpar[-tiltpar]
        alpha_bar<-1-alpha
        CDF<-do.call(pdistnam,c(list(data),as.list(distpar)))
        theop<-CDF/(1-alpha_bar*(1-CDF))
        n <- length(theop)
        -3*n/2 + 2 * sum(theop) - mean ( (2 * 1:n - 1) * log(theop))
        }
  else
    if (gof == "AD2R")
      fnobj <- function(allpar,pdistnam,data)
      {
        if(!is.list(allpar))
          allpar<- as.list(allpar)
        tiltpar<-match("tilt",names(allpar))
        alpha<-allpar$tilt
        distpar<-allpar[-tiltpar]
        alpha_bar<-1-alpha
        CDF<-do.call(pdistnam,c(list(data),as.list(distpar)))
        theop<-CDF/(1-alpha_bar*(1-CDF))
        n <- length(theop)
        2 * sum(log(1 - theop)) + mean ( (2 * 1:n - 1) / (1 - rev(theop)))
        }
  else
    if (gof == "AD2L")
      fnobj <- function(allpar,pdistnam,data)
        {
        if(!is.list(allpar))
          allpar<- as.list(allpar)
        tiltpar<-match("tilt",names(allpar))
        alpha<-allpar$tilt
        distpar<-allpar[-tiltpar]
        alpha_bar<-1-alpha
        CDF<-do.call(pdistnam,c(list(data),as.list(distpar)))
        theop<-CDF/(1-alpha_bar*(1-CDF))
        n <- length(theop)
        2 * sum(log(theop)) + mean ( (2 * 1:n - 1) / theop )
        }
  else
    if (gof == "AD2")
      fnobj <- function(allpar,pdistnam,data)
        {
        if(!is.list(allpar))
          allpar<- as.list(allpar)
        tiltpar<-match("tilt",names(allpar))
        alpha<-allpar$tilt
        distpar<-allpar[-tiltpar]
        alpha_bar<-1-alpha
        CDF<-do.call(pdistnam,c(list(data),as.list(distpar)))
        theop<-CDF/(1-alpha_bar*(1-CDF))
        n <- length(theop)
        2 * sum(log(theop) + log(1 - theop) ) + mean ( ((2 * 1:n - 1) / theop) + ((2 * 1:n - 1) / (1 - rev(theop))) )
        }
  opt <- optim(par=start, fn=fnobj,pdistnam=pdistname,data=data,hessian=TRUE,...)
  #---------------
  if(!is.list(opt$par))
    opt_par<- as.list(opt$par)
  tiltpar<-match("tilt",names(opt_par))
  alpha<-opt_par$tilt
  par<-opt_par[-tiltpar]
  alpha_bar<-1-alpha

  CDF<-do.call(pdistname,c(list(data),as.list(par)))
  PDF<-do.call(ddistname,c(list(data),as.list(par)))
  theod<-(alpha*PDF)/(1-alpha_bar*(1-CDF))^2
  loglik <- sum(log(theod))


      if(is.null(names(opt$par)))
          names(opt$par) <- names(start)
  res <- list(estimate = opt$par, convergence = opt$convergence, value = opt$value,hessian = opt$hessian,counts=opt$counts, optim.message=opt$message,loglik=loglik, gof=gof)

  return(res)
}
