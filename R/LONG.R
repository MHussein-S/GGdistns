########## LONG ########################
dLONG<-function (x,par,distr, log = FALSE)
{
  if(!is.list(par)|is.null(names(par)))
    stop("'par' must be a named list")
  ddistname <- paste("d", distr, sep = "")
  pdistname <- paste("p", distr, sep = "")
  if (!exists(ddistname, mode="function"))
    stop(paste("The ", distr, " distribution is not defined"))
  apar<-match("a",names(par))
  if(is.na(apar))
    stop("the extra parameter 'a' is not specified")
  args <- names(formals(ddistname))
  a<-par$a
  if (! a > 0)
    stop("LONG distribution is only defined for a > 0")
  distparn<-setdiff(names(par),c("a"))
  distpar<-par[distparn]
  m <- match(distparn,args)
  if (any(is.na(m)))
    stop("you specified invalid names of parameters for ",distr)
  g<-do.call(ddistname, c(list(x), as.list(distpar)))
  G<-do.call(pdistname, c(list(x), as.list(distpar)))
  d<-g/(a*sqrt(pi)*G*(1-G))*exp(-1/(2*a^2)*(log(G/(1-G)))^2)
  if (log)
    d <- log(d)
  return(d)
}

pLONG<-function (q,par,distr, lower.tail = TRUE, log.p = FALSE )
{
  if(!is.list(par)|is.null(names(par)))
    stop("'par' must be a named list")
  pdistname <- paste("p", distr, sep = "")
  if (!exists(pdistname, mode="function"))
    stop(paste("The ", distr, " distribution is not defined"))
  apar<-match("a",names(par))
  if(is.na(apar))
    stop("the extra parameter 'a' is not  specified")
  args <- names(formals(pdistname))
  a<-par$a
  if (! a > 0)
    stop("LONG distribution is only defined for a > 0")
  distparn<-setdiff(names(par),c("a"))
  distpar<-par[distparn]
  m <- match(distparn,args)
  if (any(is.na(m)))
    stop("you specified invalid names of parameters for ",distr)
  G<-do.call(pdistname, c(list(q), as.list(distpar)))
  p<-1/2*(1+pracma::erf(1/(a*sqrt(2))*log(G/(1-G))))
  if (!lower.tail)
    p <- 1 - p
  if (log.p)
    p <- log(p)
  return(p)
}
qLONG<-function(p,par,distr,lower.tail=TRUE,log.p=FALSE)
{
  if(!is.list(par)|is.null(names(par)))
    stop("'par' must be a named list")
  qdistname <- paste("q", distr, sep = "")
  if (!exists(qdistname, mode="function"))
    stop(paste("The ", distr, " distribution is not defined"))
  apar<-match("a",names(par))
  if(is.na(apar))
    stop("the extra parameter 'a' is not specified")
  args <- names(formals(qdistname))
  a<-par$a
  if (! a > 0)
    stop("LONG distribution is only defined for a > 0")
  distparn<-setdiff(names(par),c("a"))
  distpar<-par[distparn]
  m <- match(distparn,args)
  if (any(is.na(m)))
    stop("you specified invalid names of parameters for ",distr)
  if (log.p==TRUE)
    p<-exp(p)
  if (lower.tail==FALSE)
    p<-1-p
  p<-exp(a*sqrt(2)*pracma::erfinv(2*p-1))/(1+exp(a*sqrt(2)*pracma::erfinv(2*p-1)))
  q<-do.call(qdistname,c(list(p), as.list(distpar)))
  return(q)
}
rLONG<-function(n,par,distr)
{
  if(!is.list(par)|is.null(names(par)))
    stop("'par' must be a named list")
  qdistname <- paste("q", distr, sep = "")
  if (!exists(qdistname, mode="function"))
    stop(paste("The ", distr, " distribution is not defined"))
  apar<-match("a",names(par))
  if(is.na(apar))
    stop("the extra parameters 'a' is not specified")
  args <- names(formals(qdistname))
  a<-par$a
  if (! a > 0)
    stop("LONG distribution is only defined for a > 0")
  distparn<-setdiff(names(par),c("a","b"))
  distpar<-par[distparn]
  m <- match(distparn,args)
  if (any(is.na(m)))
    stop("you specified invalid names of parameters for ",distr)
  y<-runif(n)
  y<-exp(a*sqrt(2)*pracma::erfinv(2*y-1))/(1+exp(a*sqrt(2)*pracma::erfinv(2*y-1)))
  do.call(qdistname,c(list(y), as.list(distpar)))
}

mgofLONG<- function (data, distr, start,gofs = "CvM")
{
  if(!is.numeric(start)|is.null(names(start))|any(names(start)==""))
    stop("'start' must be a named numeric vector'")
  if (!is.character(distr))
    stop("distr must be a character string naming the baseline distribution")
  ddistname <- paste("d",distr,sep="")
  pdistname <- paste("p",distr,sep="")
  data<-sort(data)
  gofs <- match.arg(gofs, c("CvM", "KS", "AD", "ADR", "ADL", "AD2R", "AD2L", "AD2"))
  if (gofs == "CvM")
    fnobj <- function(allpar,pdistname,data)
    {
      parset<-setpar(pdistname,allpar=as.list(allpar))
      a<-parset$extrpar$a
      distpar<-parset$distpar
      G<-do.call(pdistname,c(list(data),as.list(distpar)))
      theop<-1/2*(1+pracma::erf(1/(a*sqrt(2))*log(G/(1-G))))
      n<-length(data)
      1/(12*n) + sum( ( theop - (2 * 1:n - 1)/(2 * n) )^2 )
    }
  else
    if (gofs == "KS")
      fnobj <- function(allpar,pdistname,data)
      {
        parset<-setpar(pdistname,allpar=as.list(allpar))
        a<-parset$extrpar$a
        distpar<-parset$distpar
        G<-do.call(pdistname,c(list(data),as.list(distpar)))
        theop<-1/2*(1+pracma::erf(1/(a*sqrt(2))*log(G/(1-G))))
        n<-length(data)
        obspu <- seq(1,n)/n
        obspl <- seq(0,n-1)/n
        max(pmax(abs(theop-obspu),abs(theop-obspl)))
      }
  else
    if (gofs == "AD")
      fnobj <- function(allpar,pdistname,data)
      {
        parset<-setpar(pdistname,allpar=as.list(allpar))
        a<-parset$extrpar$a
        distpar<-parset$distpar
        G<-do.call(pdistname,c(list(data),as.list(distpar)))
        theop<-1/2*(1+pracma::erf(1/(a*sqrt(2))*log(G/(1-G))))
        n<-length(data)
        - n - mean( (2 * 1:n - 1) * (log(theop) + log(1 - rev(theop))))
      }
  else
    if (gofs == "ADR")
      fnobj <- function(allpar,pdistname,data)
      {
        parset<-setpar(pdistname,allpar=as.list(allpar))
        a<-parset$extrpar$a
        distpar<-parset$distpar
        G<-do.call(pdistname,c(list(data),as.list(distpar)))
        theop<-1/2*(1+pracma::erf(1/(a*sqrt(2))*log(G/(1-G))))
        n<-length(data)
        n/2 - 2 * sum(theop) - mean ( (2 * 1:n - 1) * log(1 - rev(theop)))
      }
  else
    if (gofs == "ADL")
      fnobj <- function(allpar,pdistname,data)
      {
        parset<-setpar(pdistname,allpar=as.list(allpar))
        a<-parset$extrpar$a
        distpar<-parset$distpar
        G<-do.call(pdistname,c(list(data),as.list(distpar)))
        theop<-1/2*(1+pracma::erf(1/(a*sqrt(2))*log(G/(1-G))))
        n<-length(data)
        -3*n/2 + 2 * sum(theop) - mean ( (2 * 1:n - 1) * log(theop))
      }
  else
    if (gofs == "AD2R")
      fnobj <- function(allpar,pdistname,data)
      {
        parset<-setpar(pdistname,allpar=as.list(allpar))
        a<-parset$extrpar$a
        distpar<-parset$distpar
        G<-do.call(pdistname,c(list(data),as.list(distpar)))
        theop<-1/2*(1+pracma::erf(1/(a*sqrt(2))*log(G/(1-G))))
        n<-length(data)
        2 * sum(log(1 - theop)) + mean ( (2 * 1:n - 1) / (1 - rev(theop)))
      }
  else
    if (gofs == "AD2L")
      fnobj <- function(allpar,pdistname,data)
      {
        parset<-setpar(pdistname,allpar=as.list(allpar))
        a<-parset$extrpar$a
        distpar<-parset$distpar
        G<-do.call(pdistname,c(list(data),as.list(distpar)))
        theop<-1/2*(1+pracma::erf(1/(a*sqrt(2))*log(G/(1-G))))
        n<-length(data)
        2 * sum(log(theop)) + mean ( (2 * 1:n - 1) / theop )
      }
  else
    if (gofs == "AD2")
      fnobj <- function(allpar,pdistname,data)
      {
        parset<-setpar(pdistname,allpar=as.list(allpar))
        a<-parset$extrpar$a
        distpar<-parset$distpar
        G<-do.call(pdistname,c(list(data),as.list(distpar)))
        theop<-1/2*(1+pracma::erf(1/(a*sqrt(2))*log(G/(1-G))))
        n<-length(data)
        2 * sum(log(theop) + log(1 - theop) ) + mean ( ((2 * 1:n - 1) / theop) + ((2 * 1:n - 1) / (1 - rev(theop))) )
      }
  npar <- length(start)
  Mat <- diag(npar)
  opt <- constrOptim(theta=start,f=fnobj,pdistname=pdistname, data = data ,grad=NULL, ui=Mat,ci = 0)
  #---------------
  parset<-setpar(ddistname,allpar=as.list(opt$par))
  a<-parset$extrpar$a
  distpar<-parset$distpar
  G<-do.call(pdistname,c(list(data),as.list(distpar)))
  g<-do.call(ddistname,c(list(data),as.list(distpar)))
  theod<-g/(a*sqrt(pi)*G*(1-G))*exp(-1/(2*a^2)*(log(G/(1-G)))^2)
  loglik <- sum(log(theod))
  n<-length(data)
  AIC<--2*loglik+2*npar
  AICc <- -2 * loglik + 2 * npar + 2 * (npar * (npar + 1))/(n - npar - 1)
  BIC<- -2*loglik+ npar*log(n)
  HQIC<- -2*loglik+2*npar*log(log(n))
  res=cbind(opt$par)
  colnames(res)=c("Estimate")
  res1=cbind(AIC,AICc,BIC, HQIC, loglik)
  colnames(res1)=c("AC","CAC","BC","HQC","log(Likelihood)")
  rownames(res1)=c("")
  res2=cbind(if(opt$convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
  colnames(res2)=c("")
  rownames(res2)=c("")
  res3=cbind(gofs,opt$value)
  colnames(res3)=c("name","value")
  rownames(res3)=c("")
  list("Estimates"=res,"Measures"=res1,"Convergence Status"=res2, "Goodness-of-fit Statistic"=res3)
}



