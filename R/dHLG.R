#pdf of Half logistic G for any base distribution
dHLG<-function (x, distr,par,log=FALSE)
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
  ddistname <- paste("d",distr,sep="")
  if (!exists(pdistname, mode="function"))
    stop(paste("The ", distr, " distribution is not defined"))
  parnm <- names(par)
  args <- names(formals(pdistname))
  m <- match(parnm,args)
  if (any(is.na(m)))
    stop("'par' specifies names which are not arguments to ",distr)
#--------
  G<-do.call(pdistname,c(list(x), as.list(par)))
  g<-do.call(ddistname,c(list(x), as.list(par)))
  d<-2*a*g*G^(a-1)
  d<-d/(1-G^a)^2
  if (log)
    d <- log(d)
  return(d)
}
