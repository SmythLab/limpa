dpcImpute <- function(y, ...)
  UseMethod("dpcImpute")

dpcImpute.default <- function(y, dpc=NULL, dpc.slope=0.8, verbose=TRUE, chunk=1000L, ...)
# Use the DPC to quantify protein expression values by maximum posterior.
# Created 25 Feb 2025. Last modified 25 Feb 2025.
{
# Check y
  y <- as.matrix(y)

# Construct EList and pass to EList method
  z <- new("EList",list(E=y))
  dpcImpute(z,dpc=dpc,dpc.slope=dpc.slope,verbose=verbose,chunk=chunk,...)
}

dpcImpute.EList <- function(y, dpc=NULL, dpc.slope=0.8, verbose=TRUE, chunk=1000L, ...)
# Use the DPC to quantify protein expression values by maximum posterior.
# Created 27 Dec 2024. Last modified 9 Jan 2024.
{
# Check dpc
  if(is.list(dpc)) dpc <- dpc$dpc
  if(is.null(dpc)) {
    beta0 <- estimateDPCIntercept(y, dpc.slope=dpc.slope)
    if(verbose) message("DPC intercept estimated as ", formatC(beta0))
    dpc <- c(beta0=beta0, beta1=dpc.slope)
  }

# Estimate Bayes hyperparameters
  if(verbose) message("Estimating hyperparameters ...")
  h <- dpcImputeHyperparam(y, dpc.slope=dpc[2], ...)

# Summarize peptides to proteins by maximum posterior estimation
  if(verbose) message("Imputing peptides ...")
  protein.id <- seq_len(nrow(y))
  y.protein <- peptides2Proteins(y,
                          protein.id = protein.id,
                          dpc = dpc,
                          sigma = 0.01,
                          prior.mean = h$prior.mean,
                          prior.sd = h$prior.sd,
                          prior.logFC = h$prior.logFC,
                          standard.errors = TRUE,
                          verbose = verbose,
                          chunk = chunk)

# Add back original original annotation
  if(!is.null(y$genes)) y.protein$genes <- data.frame(y$genes,y.protein$genes)

# Output
  rownames(y.protein) <- rownames(y)
  y.protein$targets <- y$targets
  y.protein$dpc <- dpc
  y.protein
}
