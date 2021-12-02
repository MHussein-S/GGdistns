#Maximum goodness-of-fit estimates of Exponentiated Exponential Poisson G distribution "using constroptim"
# start is must be a named vector ex: EEPGmgof1(data, distr="weibull", gof = "AD", start=c(a=0.5,lambda=0.5,scale=0.5,shape=0.5))
EEPGmgof1 <- function (data, distr, start,gof = "CvM",...)
{
  if (!(is.numeric(data) & length(data)>1))
    stop("data must be a numeric vector of length greater than 1")
  if(is.list(start)|is.null(names(start)))
    stop("'start' must be a named atomic vector")
  if (!is.character(distr))
    stop("distr must be a character string naming a distribution")
  #------
  apar<-match("a",names(start))
  if (is.na(apar))
    stop("Starting value of \"a\" is not defined")
  a<-start["a"]
  par<-start[-apar]

  lpar<-match("lambda",names(par))
  if (is.na(lpar))
    stop("Starting value of \"lambda\" is not defined")
  lambda<-start["lambda"]
  par<-par[-lpar]
  #--------
  ddistname <- paste("d",distr,sep="")
  pdistname <- paste("p",distr,sep="")
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
      apar<-match("a",names(allpar))
      a<-allpar$a
      distpar<-allpar[-apar]
      lpar<-match("lambda",names(distpar))
      lambda<-allpar$lambda
      distpar<-distpar[-lpar]
      G<-do.call(pdistnam,c(list(data),as.list(distpar)))
      theop<-(1-exp(-lambda))^-1*(1-exp(-lambda*G^a))
      n<-length(theop)
      1/(12*n) + sum( ( theop - (2 * 1:n - 1)/(2 * n) )^2 )
      }
  else
    if (gof == "KS")
      fnobj <- function(allpar,pdistnam,data)
      {
        if(!is.list(allpar))
          allpar<- as.list(allpar)
        apar<-match("a",names(allpar))
        a<-allpar$a
        distpar<-allpar[-apar]
        lpar<-match("lambda",names(distpar))
        lambda<-allpar$lambda
        distpar<-distpar[-lpar]
        G<-do.call(pdistnam,c(list(data),as.list(distpar)))
        theop<-(1-exp(-lambda))^-1*(1-exp(-lambda*G^a))
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
        apar<-match("a",names(allpar))
        a<-allpar$a
        distpar<-allpar[-apar]
        lpar<-match("lambda",names(distpar))
        lambda<-allpar$lambda
        distpar<-distpar[-lpar]
        G<-do.call(pdistnam,c(list(data),as.list(distpar)))
        theop<-(1-exp(-lambda))^-1*(1-exp(-lambda*G^a))
        n<-length(theop)
        - n - mean( (2 * 1:n - 1) * (log(theop) + log(1 - rev(theop))))
        }
  else
    if (gof == "ADR")
      fnobj <- function(allpar,pdistnam,data)
      {
        if(!is.list(allpar))
          allpar<- as.list(allpar)
        apar<-match("a",names(allpar))
        a<-allpar$a
        distpar<-allpar[-apar]
        lpar<-match("lambda",names(distpar))
        lambda<-allpar$lambda
        distpar<-distpar[-lpar]
        G<-do.call(pdistnam,c(list(data),as.list(distpar)))
        theop<-(1-exp(-lambda))^-1*(1-exp(-lambda*G^a))
        n<-length(theop)
        n/2 - 2 * sum(theop) - mean ( (2 * 1:n - 1) * log(1 - rev(theop)))
        }
  else
    if (gof == "ADL")
      fnobj <- function(allpar,pdistnam,data)
      {
        if(!is.list(allpar))
          allpar<- as.list(allpar)
        apar<-match("a",names(allpar))
        a<-allpar$a
        distpar<-allpar[-apar]
        lpar<-match("lambda",names(distpar))
        lambda<-allpar$lambda
        distpar<-distpar[-lpar]
        G<-do.call(pdistnam,c(list(data),as.list(distpar)))
        theop<-(1-exp(-lambda))^-1*(1-exp(-lambda*G^a))
        n<-length(theop)
        -3*n/2 + 2 * sum(theop) - mean ( (2 * 1:n - 1) * log(theop))
        }
  else
    if (gof == "AD2R")
      fnobj <- function(allpar,pdistnam,data)
      {
        if(!is.list(allpar))
          allpar<- as.list(allpar)
        apar<-match("a",names(allpar))
        a<-allpar$a
        distpar<-allpar[-apar]
        lpar<-match("lambda",names(distpar))
        lambda<-allpar$lambda
        distpar<-distpar[-lpar]
        G<-do.call(pdistnam,c(list(data),as.list(distpar)))
        theop<-(1-exp(-lambda))^-1*(1-exp(-lambda*G^a))
        n<-length(theop)
        2 * sum(log(1 - theop)) + mean ( (2 * 1:n - 1) / (1 - rev(theop)))
        }
  else
    if (gof == "AD2L")
      fnobj <- function(allpar,pdistnam,data)
      {
        if(!is.list(allpar))
          allpar<- as.list(allpar)
        apar<-match("a",names(allpar))
        a<-allpar$a
        distpar<-allpar[-apar]
        lpar<-match("lambda",names(distpar))
        lambda<-allpar$lambda
        distpar<-distpar[-lpar]
        G<-do.call(pdistnam,c(list(data),as.list(distpar)))
        theop<-(1-exp(-lambda))^-1*(1-exp(-lambda*G^a))
        n<-length(theop)
        2 * sum(log(theop)) + mean ( (2 * 1:n - 1) / theop )
        }
  else
    if (gof == "AD2")
      fnobj <- function(allpar,pdistnam,data)
      {
        if(!is.list(allpar))
          allpar<- as.list(allpar)
        apar<-match("a",names(allpar))
        a<-allpar$a
        distpar<-allpar[-apar]
        lpar<-match("lambda",names(distpar))
        lambda<-allpar$lambda
        distpar<-distpar[-lpar]
        G<-do.call(pdistnam,c(list(data),as.list(distpar)))
        theop<-(1-exp(-lambda))^-1*(1-exp(-lambda*G^a))
        n<-length(theop)
        2 * sum(log(theop) + log(1 - theop) ) + mean ( ((2 * 1:n - 1) / theop) + ((2 * 1:n - 1) / (1 - rev(theop))) )
      }
  npar <- length(start)
  Mat <- diag(npar)
  colnames(Mat) <- names(start)
  rownames(Mat) <- paste0("constr", 1:npar)
  initconstr <- Mat %*% start - seq(0,0,npar)
  if (any(initconstr < 0))
    stop("Starting values must be in the feasible region.")
  opt <- constrOptim(theta=start,f=fnobj,pdistnam=pdistname, data = data ,grad=NULL, ui=Mat,ci = seq(0,0,npar))
  #---------------
  if(!is.list(opt$par))
    opt_par<- as.list(opt$par)

  apar<-match("a",names(opt_par))
  a<-opt_par$a
  par<-opt_par[-apar]
  lpar<-match("lambda",names(par))
  lambda<-par$lambda
  par<-par[-lpar]

  G<-do.call(pdistname,c(list(data),as.list(par)))
  g<-do.call(ddistname,c(list(data),as.list(par)))
  theod<-a*lambda*(1-exp(-lambda))^-1*g*G^(a-1)*exp(-lambda*G^a)
  loglik <- sum(log(theod))


      if(is.null(names(opt$par)))
          names(opt$par) <- names(start)
  res <- list(estimate = opt$par, convergence = opt$convergence, value = opt$value,hessian = opt$hessian,counts=opt$counts, optim.message=opt$message,loglik=loglik, gof=gof)

  return(res)
}
