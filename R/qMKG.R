#Quantile function of Modified Kies G for any base distribution
qMKG<-function (p, distr,par,lower.tail=TRUE,log.p=FALSE)
{
  if(!is.list(par))
    par<- as.list(par)
  if (!is.character(distr))
    stop("distr must be a character string naming a distribution")
  if(is.null(names(par)))
  	stop("Parameters must be a named list or vector.")
  apar<-match("a",names(par))
  if (is.na(apar))
    stop("parameter \"a\" is not defined")
  apar<-match("a",names(par))
  a<-par$a
  par<-par[-apar]
#--------
  qdistname <- paste("q",distr,sep="")
  if (!exists(qdistname, mode="function"))
    stop(paste("The ", distr, " distribution is not defined"))
  parnm <- names(par)
  args <- names(formals(qdistname))
  m <- match(parnm,args)
  if (any(is.na(m)))
    stop("'par' specifies names which are not arguments to ",distr)
#--------
  if (!lower.tail)
    p <- 1 - p
  if (log.p)
    p <- exp(p)
  p1<-(-log(1-p))^(1/a)
  p1<-p1/(1+p1)
  q<-do.call(qdistname,c(list(p1), as.list(par)))
}
