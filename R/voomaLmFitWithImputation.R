voomaLmFitWithImputation <- function(
	y, design=NULL, prior.weights=NULL, imputed=NULL,
	block=NULL,
	sample.weights=FALSE, var.design=NULL, var.group=NULL, prior.n=10,
	predictor=NULL, span=NULL, legacy.span=FALSE, plot=FALSE, save.plot=FALSE,
	keep.EList=TRUE
)
#	vooma+lmFit with iteration of sample weights and intrablock correlation.
#	Creates an MArrayLM object for entry to eBayes() etc in the limma pipeline.
#	Corrects for loss of residual df due to entirely imputed values in a group.
#	Mengbo Li and Gordon Smyth
#	Created 24 Nov 2023. Last modifed 10 Apr 2025.
{
	Block <- !is.null(block)
	PriorWeights <- !is.null(prior.weights)
	SampleWeights <- sample.weights || !is.null(var.design) || !is.null(var.group)
	
#	Can't specify prior weights and ask for sample weights to be estimated as well
	if(PriorWeights && SampleWeights) stop("Can't specify prior.weights and estimate sample weights")
	
#	Check y
	if(!is(y,"EList")) y <- new("EList",list(E=as.matrix(y)))
	narrays <- ncol(y)
	ngenes <- nrow(y)
	if(narrays < 2L) stop("Too few samples")
	if(ngenes < 2L) stop("Need multiple rows")
	A <- rowMeans(y$E,na.rm=TRUE)
	if(anyNA(A)) stop("y contains entirely NA rows")

#	Check imputed
	if(!is.null(imputed)) {
		if(!identical(nrow(imputed),ngenes) || !identical(ncol(imputed),narrays)) stop("imputed must have same dimensions as y")
		if(anyNA(y$E)) stop("y is imputed but also has NAs")
		NImputed <- rowSums(imputed)
		if(max(NImputed) <= 1L) imputed <- NULL
	}	

#	Check predictor
	if(!is.null(predictor)) {
		predictor <- as.matrix(predictor)
		if(!identical(nrow(predictor),ngenes)) stop("predictor is of wrong dimension")
		if(identical(ncol(predictor),1L)) {
			predictor <- matrix(predictor,ngenes,narrays)
		} else {
			if(!identical(ncol(predictor),narrays)) stop("predictor is of wrong dimension")
		}
		if(is.infinite(min(predictor,na.rm=TRUE))) stop("predictor contains infinite values")
		if(is.infinite(max(predictor,na.rm=TRUE))) stop("predictor contains infinite values")
		if(anyNA(predictor)) {
			if(anyNA(y$E)) {
				if(anyNA( predictor[!is.na(y$E)] )) stop("All observed y values must have non-NA predictors")
			} else {
				stop("All observed y values must have non-NA predictors")
			}
		}
	}
	
#	Check design	
	if(is.null(design)) design <- y$design
	if(is.null(design)) {
		design <- matrix(1,narrays,1)
		rownames(design) <- colnames(y)
		colnames(design) <- "GrandMean"
	}

#	Expand prior.weights if necessary
	if(!is.null(prior.weights)) prior.weights <- asMatrixWeights(prior.weights,dim(y))
	
#	Fit linear model
	fit <- lmFit(y, design, weights = prior.weights)
	
#	Find individual fitted values
	if(fit$rank < ncol(design)) {
		j <- fit$pivot[1:fit$rank]
		fitted.values <- fit$coefficients[,j,drop=FALSE] %*% t(fit$design[,j,drop=FALSE])
	} else {
		fitted.values <- fit$coefficients %*% t(fit$design)
	}

#	Correct for loss of residual degrees of freedom from imputed values
	DFLoss <- FALSE
	if(!is.null(imputed)) {
		h <- 1-hat(design)
		DFImputed <- imputed %*% (1-h)
		HasImp <- which(DFImputed > 0.9999)
		if(length(HasImp)) {
			EHasImp <- y$E[HasImp,,drop=FALSE]
			IHasImp <- imputed[HasImp,,drop=FALSE]
			fit1 <- lm.fit(design,t(EHasImp))
			EHasImp2 <- EHasImp
			EHasImp2[!IHasImp] <- EHasImp2[!IHasImp] + ncol(y) + 1
			fit2 <- lm.fit(design,t(EHasImp2))
			NoChange <- t(abs(fit2$fitted.values - fit1$fitted.values) < 0.01)
			AllImpFitVal <- IHasImp & NoChange
			HasAllImpFitVal <- which(rowSums(AllImpFitVal) > 0)
			DFLoss <- length(HasAllImpFitVal) > 0
		}
		if(DFLoss) {
			if(length(HasAllImpFitVal) < length(HasImp)) {
				HasImp <- HasImp[HasAllImpFitVal]
				EHasImp <- EHasImp[HasAllImpFitVal,,drop=FALSE]
				AllImpFitVal <- AllImpFitVal[HasAllImpFitVal,,drop=FALSE]
			}
			EHasImpNA <- EHasImp
			EHasImpNA[AllImpFitVal] <- NA
			fitNA <- suppressWarnings(lmFit(EHasImpNA,design,weights=prior.weights[HasImp,,drop=FALSE]))
			if(min(fitNA$df.residual) < 0.5) {
				NoDF <- fitNA$df.residual < 0.5
				fitNA$sigma[NoDF] <- median(fit$sigma,na.rm=TRUE)
				fitNA$Amean[NoDF] <- median(fit$Amean,na.rm=TRUE)
			}
			fit$sigma[HasImp] <- fitNA$sigma
			fit$df.residual[HasImp] <- fitNA$df.residual
			A[HasImp] <- fitNA$Amean
		}
	}

#	Prepare to fit lowess trend
	sx <- A
	sy <- sqrt(fit$sigma)
	mu <- fitted.values
	
#	Optionally combine ave log intensity with precision predictor
	if(!is.null(predictor)) {
		sxc <- rowMeans(predictor, na.rm = TRUE)
		if(DFLoss) {
			predictorNA <- predictor[HasImp,,drop=FALSE]
			predictorNA[AllImpFitVal] <- NA
			sxc[HasImp] <- rowMeans(predictorNA,na.rm=TRUE)
			if(anyNA(sxc)) sxc[is.na(sxc)] <- median(sxc,na.rm=TRUE)
		}
		vartrend <- lm.fit(cbind(1,sx,sxc),sy)
		beta <- coef(vartrend)
		sx <- vartrend$fitted.values
		mu <- beta[1] + beta[2]*mu + beta[3]*predictor
		xlab <- "Combined predictor"
		main.title <- "vooma variance trend"
	} else {
		xlab <- "Average log2 expression"
		main.title <- "vooma mean-variance trend"
	}

#	Choose span based on the number of genes
	if(is.null(span))
		if(legacy.span)
			span <- chooseLowessSpan(ngenes,small.n=10,min.span=0.3,power=0.5) 
		else
			span <- chooseLowessSpan(ngenes,small.n=50,min.span=0.3,power=1/3)
	
#	Fit lowess trend
	if(DFLoss) {
		lf <- loessFit(y=sy,x=sx,span=span,weights=fit$df.residual)
		o <- order(sx)
		l <- list(x=sx[o],y=lf$fitted[o])
	} else {
		l <- lowess(sx,sy,f=span)
	}
	if(plot) {
		plot(sx,sy,xlab=xlab,ylab="Sqrt( standard deviation )",pch=16,cex=0.25)
		title(main.title)
		lty <- ifelse(Block || SampleWeights,2,1)
		lines(l,col="red",lty=lty)
	}
	
#	Make interpolating rule
	f <- approxfun(l, rule=2, ties=list("ordered",mean))
	
#	Apply trend to individual observations to get vooma weights
	w <- 1/f(mu)^4
	dim(w) <- dim(y)
	colnames(w) <- colnames(y)
	rownames(w) <- rownames(y)
	
#	Add voom weights to prior weights
	if(PriorWeights) {
		weights <- w * prior.weights
		attr(weights,"arrayweights") <- NULL
	} else {
		weights <- w
	}
	
#	Estimate sample weights?
	if(SampleWeights) {
		sw <- arrayWeights(y,design,weights=weights,var.design=var.design,var.group=var.group,prior.n=prior.n)
		message("First sample weights (min/max) ", paste(format(range(sw)),collapse="/") )
		if(Block) weights <- t(sw * t(weights))
	}
	
#	Estimate correlation?
	if(Block) {
		dc <- suppressWarnings(duplicateCorrelation(y,design,block=block,weights=weights))
		correlation <- dc$consensus.correlation
		if(is.na(correlation)) correlation <- 0
		message("First intra-block correlation  ",format(correlation))
	} else {
		correlation <- NULL
	}
	
#	Second iteration to refine intra-block correlation or sample weights
	if(Block || SampleWeights) {
#		Rerun voom weights with new correlation and sample weights
		if(SampleWeights) {
			weights <- asMatrixWeights(sw,dim(y))
		} else {
			weights <- prior.weights
		}
		fit <- lmFit(y,design,block=block,correlation=correlation,weights=weights)
#		Find individual fitted values
		if(fit$rank < ncol(design)) {
			j <- fit$pivot[1:fit$rank]
			fitted.values <- fit$coefficients[,j,drop=FALSE] %*% t(fit$design[,j,drop=FALSE])
		} else {
			fitted.values <- fit$coefficients %*% t(fit$design)
		}
#		Correct for loss of residual degrees of freedom from imputed values
		if(DFLoss) {
			fitNA <- suppressWarnings(lmFit(EHasImpNA,design,block=block,correlation=correlation,weights=weights[HasImp,,drop=FALSE]))
			if(min(fitNA$df.residual) < 0.5) {
				fitNA$sigma[fitNA$df.residual < 0.5] <- median(fit$sigma,na.rm=TRUE)
			}
			fit$sigma[HasImp] <- fitNA$sigma
			fit$df.residual[HasImp] <- fitNA$df.residual
		}
#		Prepare to fit NEW lowess trend
		sx <- A
		sy <- sqrt(fit$sigma)
		mu <- fitted.values
#		Optionally combine ave log intensity with precision predictor for the NEW trend
		if(!is.null(predictor)) {
			vartrend <- lm.fit(cbind(1,sx,sxc),sy)
			beta <- coef(vartrend)
			sx <- vartrend$fitted.values
			mu <- beta[1] + beta[2]*mu + beta[3]*predictor
		}
		if(DFLoss) {
			lf <- loessFit(y=sy,x=sx,span=span,weights=fit$df.residual)
			o <- order(sx)
			l <- list(x=sx[o],y=lf$fitted[o])
		} else {
			l <- lowess(sx,sy,f=span)
		}
		if(plot) {
			lines(l,col="red")
			legend("topleft",lty=c(2,1),col="red",legend=c("First","Final"))
		}
#		Make NEW interpolating rule
		f <- approxfun(l, rule=2, ties=list("ordered",mean))
#		Apply NEW trend to individual observations to get vooma weights
		w <- 1/f(mu)^4
		dim(w) <- dim(y)
		colnames(w) <- colnames(y)
		rownames(w) <- rownames(y)
#		Add voom weights to prior weights
		if(PriorWeights) {
			weights <- w * prior.weights
			attr(weights,"arrayweights") <- NULL
		} else {
			weights <- w
		}
		if(SampleWeights) {
			sw <- arrayWeights(y,design,weights=weights,var.design=var.design,var.group=var.group,prior.n=prior.n)
			message("Final sample weights (min/max) ", paste(format(range(sw)),collapse="/") )
			weights <- t(sw * t(weights))
		}
		if(Block) {
			dc <- suppressWarnings(duplicateCorrelation(y,design,block=block,weights=weights))
			correlation <- dc$consensus.correlation
			if(is.na(correlation)) correlation <- 0
			message("Final intra-block correlation  ",format(correlation))
		}
	}
	
#	Final linear model fit with voom weights
	fit <- lmFit(y,design,block=block,correlation=correlation,weights=weights)
	if(DFLoss) {
		fitNA <- suppressWarnings(lmFit(EHasImpNA,design,block=block,correlation=correlation,weights=weights[HasImp,,drop=FALSE]))
		if(min(fitNA$df.residual) < 0.5) {
			fitNA$sigma[fitNA$df.residual < 0.5] <- median(fit$sigma,na.rm=TRUE)
		}
		fit$sigma[HasImp] <- fitNA$sigma
		fit$df.residual[HasImp] <- fitNA$df.residual
	}

#	Output
	if(!is.null(y$genes)) {
		fit$genes <- y$genes
	}
	if(SampleWeights) {
		fit$targets <- y$targets
		if(is.null(fit$targets)) {
			fit$targets <- data.frame(sample.weight=sw)
			row.names(fit$targets) <- colnames(y)
		} else {
			fit$targets$sample.weight <- sw
		}
	}
	fit$span <- span
	if(save.plot) {
		fit$voom.xy <- list(x=sx,y=sy,xlab=xlab,ylab="Sqrt( standard deviation )",pch=16,cex=0.25)
		fit$voom.line <- l
		fit$voom.line$col <- "red"
	}
	if(keep.EList) {
		y$weights <- weights
		fit$EList <- y
	}
	fit
}
