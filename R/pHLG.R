#cdf of Half logistic G for any base distribution
pHLG<-function (q, distr,par,lower.tail=TRUE,log.p=FALSE)
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
  pdistname <- paste("p",distr,sep="")
  if (!exists(pdistname, mode="function"))
    stop(paste("The ", distr, " distribution is not defined"))
  parnm <- names(par)
  args <- names(formals(pdistname))
  m <- match(parnm,args)
  if (any(is.na(m)))
    stop("'par' specifies names which are not arguments to ",distr)
#--------
  G<-do.call(pdistname,c(list(q), as.list(par)))
  p<-2*G^a/(1+G^a)
  if (!lower.tail)
    p <- 1 - p
  if (log.p)
    p <- log(p)
  return(p)
}
