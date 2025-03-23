plotMDSUsingSEs <- function(y,top=500,labels=NULL,pch=NULL,cex=1,dim.plot=c(1,2),xlab=NULL,ylab=NULL,plot=TRUE,var.explained=TRUE,...)
#	Multi-dimensional scaling with top-distance, using standard errors
#	23 Mar 2025.  Last modified 23 Mar 2025.
{
#	Check whether y is an EList with standard errors
	x <- y$E
	if(is.null(x)) stop("y should be an EList object")
	SE <- y$other$standard.error
	if(is.null(SE)) stop("y should contain standard errors")
	SE2 <- SE^2

#	Check dimensions
	nsamples <- ncol(x)
	if(nsamples < 3) stop(paste("Only",nsamples,"columns of data: need at least 3"))
	cn <- colnames(x)

#	Remove rows with missing or Inf values
	bad <- rowSums(is.finite(x)) < nsamples
	if(any(bad)) x <- x[!bad,,drop=FALSE]
	nprobes <- nrow(x)

#	Check top
	top <- min(top,nprobes)

#	Check labels and pch
	if(is.null(pch) & is.null(labels)) {
		labels <- colnames(x)
		if(is.null(labels)) labels <- 1:nsamples
	}
	if(!is.null(labels)) labels <- as.character(labels)

#	Check dim.plot
	dim.plot <- unique(as.integer(dim.plot))
	if(length(dim.plot) != 2L) stop("dim.plot must specify two dimensions to plot")

#	Check dim
	ndim <- max(dim.plot)
	if(ndim < 2L) stop("Need at least two dim.plot")
	if(nsamples < ndim) stop("ndim is greater than number of samples")
	if(nprobes < ndim) stop("ndim is greater than number of rows of data")

#	Distance matrix from pairwise leading fold changes
#	Distance measure is mean of top squared deviations for each pair of arrays
	dd <- matrix(0,nrow=nsamples,ncol=nsamples,dimnames=list(cn,cn))
	topindex <- nprobes-top+1L
	for (i in 2L:(nsamples))
	for (j in 1L:(i-1L)) {
		z2 <- (x[,i]-x[,j])^2 / (SE2[,i] + SE2[,j])
		dd[i,j] <- mean(sort.int(z2,partial=topindex)[topindex:nprobes])
	}
	axislabel <- "Leading z-statistic dim"

#	Multi-dimensional scaling
	dd <- dd + t(dd)
	rm <- rowMeans(dd)
	dd <- dd - rm
	dd <- t(dd) - (rm - mean(rm))
	mds <- eigen(-dd/2, symmetric=TRUE)
	names(mds) <- c("eigen.values","eigen.vectors")

#	Make MDS object
	lambda <- pmax(mds$eigen.values,0)
	mds$var.explained <- lambda / sum(lambda)
	mds$dim.plot <- dim.plot
	mds$distance.matrix.squared <- dd
	mds$top <- top
	mds$gene.selection <- "pairwise"
	mds$axislabel <- axislabel
	mds <- new("MDS",unclass(mds))

#	Add coordinates for plot
	i <- dim.plot[1]
	mds$x <- mds$eigen.vectors[,i] * sqrt(lambda[i])
	if(lambda[i] < 1e-13) warning("dimension ", i, " is degenerate or all zero")
	i <- dim.plot[2]
	mds$y <- mds$eigen.vectors[,i] * sqrt(lambda[i])
	if(lambda[i] < 1e-13) warning("dimension ", i, " is degenerate or all zero")

	if(plot)
		plotMDS(mds,labels=labels,pch=pch,cex=cex,xlab=xlab,ylab=ylab,var.explained=var.explained,...)
	else
		mds
}
