#new version (working)
plotMOE2<-function (data, distr, para, histo = TRUE, breaks = "default", demp = FALSE,HRF=FALSE,PHRF=FALSE, ...) 
{
def.par <- par(no.readonly = TRUE)
if (missing(data) || !is.vector(data, mode = "numeric")) 
stop("data must be a numeric vector")

if ((missing(distr) & !missing(para)) || (!missing(distr) &  missing(para)))   
stop("distr and para -must defined")
   
if (!histo & !demp) 
stop("one the arguments histo and demp must be put to TRUE")
xlim <- c(min(data), max(data))
s <- sort(data)
n <- length(data)
if (missing(distr)) 
{
par(mfrow = c(1, 2))    
obsp<-ppoints(s)
	if (histo) 
	{
		if (demp) 
		{
		lines(density(data)$x, density(data)$y, lty = 2, col 	= 					"black")
			if (breaks == "default") 
			h <- hist(data, freq = FALSE, xlab = "Data", main = 						"Empirical density", ...)
			else 
			h <- hist(data, freq = FALSE, xlab = "Data", main = 						"Empirical density", breaks = breaks,   ...)
 		}
            	else 
		{
                  		if (breaks == "default") 
                    		h <- hist(data, freq = FALSE, xlab = "Data", main = 						"Histogram", ...)
                  		else 
			h <- hist(data, freq = FALSE, xlab = "Data", main = 						"Histogram", 	breaks	= breaks, ...)
                	}
	}
	else 
	{
            h <- hist(data, freq = FALSE, xlab = "Data", main = "Histogram", plot = 			FALSE, ...)
            plot(density(data)$x, density(data)$y, lty = 1, col = "black", type = "l", 			xlab = 	"Data", main = paste("Empirical density"), ylab = "Density", ...)
            }
plot(s, obsp, main = paste("Cumulative distribution"), xlab = "Data", xlim = c(h$breaks[1], h$breaks[length(h$breaks)]), ylab = "CDF", ...)
}
else #distr not missing
{
	if (!is.list(para)||is.null(names(para))) 
	stop("'para' must be a named list")
alpha<-para$tilt
alpha_bar=1-alpha
tiltpar<-match("tilt",names(para))
para<-para[-tiltpar]
ddistname <- paste("d", distr, sep = "")
pdistname <- paste("p", distr, sep = "")
qdistname <- paste("q", distr, sep = "")
	if (!exists(ddistname, mode = "function")) 
	paste("The " , distr, " disribution", "is not defined")
parnm <- names(para)
args <- names(formals(ddistname))
m <- match(parnm,args)
if (any(is.na(m))) 
	  stop("'para' specifies names which are not arguments to ",distr)
par(mfrow = c(2, 2))
if (any(HRF,PHRF))
par(mfrow = c(3, 2))
	if (breaks == "default") 
	h <- hist(data, plot = FALSE)
	else 
	h <- hist(data, breaks = breaks, plot = FALSE, ...)
obsp <- ppoints(s)
xhist <- seq(min(h$breaks), max(h$breaks), length = 1000)
ycdf <- do.call(pdistname, c(list(xhist), as.list(para)))
ypdf <- do.call(ddistname, c(list(xhist), as.list(para)))
MOCDF<-ycdf/(1-alpha_bar*(1-ycdf))
MOPDF<-(alpha*ypdf)/(1-alpha_bar*(1-ycdf))^2
yhist <- MOPDF
	if (length(yhist) != length(xhist)) 
	stop("problem when computing densities.")
ymax <- ifelse(is.finite(max(yhist)), max(max(h$density), max(yhist)), max(h$density))
	if (histo) 
	{
	hist(data, freq = FALSE, xlab = "Data", ylim = c(0, ymax), breaks = h			$breaks, main = paste("Empirical and theoretical dens."),    ...)
		if (demp) 
		lines(density(data)$x, density(data)$y, lty = 2, col = "black")
           	}
            else 
	plot(density(data)$x, density(data)$y, lty = 2, col = "black", type = "l", 			xlab 	= "Data", main = paste("Empirical and theoretical dens."), ylab 			= "Density", 	xlim = c(min(h$breaks), max(h$breaks)), ...)
            if (demp) 
 	legend("topright", bty = "n", lty = c(2,1), col = c("black", "red"), legend 			= 	c("empirical","theoretical"), bg = "white", cex = 0.7)
lines(xhist, yhist, lty = 1, col = "red")
obsp1<-alpha*obsp/(1-alpha_bar*obsp)
MOQF<-do.call(qdistname, c(list(obsp1), as.list(para)))
	if (length(MOQF) != length(obsp)) 
	stop("problem when computing quantities.")
plot(MOQF, s, main = " Q-Q plot", xlab = "Theoretical quantiles",  ylab ="Empirical quantiles", ...)
abline(0, 1,col = "red")
xmin <- h$breaks[1]
xmax <- h$breaks[length(h$breaks)]
	if (length(s) != length(obsp)) 
	stop("problem when computing probabilities.")
obsp<-ppoints(s)
plot(s, obsp, main = paste("Empirical and theoretical CDFs"),  xlab = "Data", 	ylab = "CDF", xlim = c(xmin, xmax), ...)
sfin <- seq(xmin, xmax, by = (xmax - xmin)/100)
ycdf_fin<- do.call(pdistname, c(list(sfin), as.list(para)))
MOCDF_fin<-ycdf_fin/(1-alpha_bar*(1-ycdf_fin))
lines(sfin, MOCDF_fin, lty = 1, col = "red")
ycdf_s<- do.call(pdistname, c(list(s), as.list(para)))
MOCDF_s<-ycdf_s/(1-alpha_bar*(1-ycdf_s))
	if (length(MOCDF_s) != length(obsp)) 
	stop("problem when computing probabilities.")
plot(MOCDF_s, obsp, main = "P-P plot", xlab = "Theoretical probabilities", ylab = "Empirical probabilities", ...)
abline(0, 1, col="red")
	if (HRF)
	{
	xdata <- seq(min(data), max(data), length = 1000)
	ycdf <- do.call(pdistname, c(list(xdata), as.list(para)))
	ypdf <- do.call(ddistname, c(list(xdata), as.list(para)))
	MOHRF<-ypdf/((1-ycdf)*(1-alpha_bar*(1-ycdf)))
	MOHRF<-format(MOHRF,7)
	plot(xdata, MOHRF,typ='l',col='red', main = "Hazard Rate Function", xlab = 	"Data", ylab = "hazard rate", ...)
	}
	if (PHRF)
	{
	xdata <- seq(min(data), max(data), length = 1000)
	MOPHRF<-alpha*ypdf/(ycdf*(1-alpha_bar*(1-ycdf)))
	MOPHRF<-format(MOPHRF,7)
	plot(xdata, MOPHRF,typ='l',col='red', main = "Proportional Hazard Rate 	Function", xlab = "Data", ylab = "p. hazard rate", ...)
	}
}
par(def.par)
invisible()
}
