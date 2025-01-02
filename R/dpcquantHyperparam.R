dpcQuantHyperparam <- function(y, protein.id, dpc.slope=0.7, sd.quantile.for.logFC = 0.9, robust=FALSE, ...)
# Estimate hyperparameters needed for dpcQuant
# Mengbo Li and Gordon Smyth
# Created 12 Sept 2024. Last modified 29 Dec 2024.
{
  # Check input
  y <- as.matrix(y)
  nsamples <- ncol(y)

  # Check that peptides are in protein order
  if( !identical(protein.id,sort(protein.id)) ) {
    o <- order(protein.id)
    protein.id <- protein.id[o]
    y <- y[o,,drop=FALSE]
  }

  # Impute by exponential tilting
  yComplete <- imputeByExpTilt(y, dpc.slope=dpc.slope, ...)

  # Prior mean
  RowMeans <- rowMeans(yComplete)
  prior.mean <- mean(RowMeans)

  # Prior standard deviation
  prior.sd <- sd(RowMeans)

  # Prior logfc
  RowVar <- rowSums((yComplete-RowMeans)^2)/(nsamples-1)
  prior.logFC <- sqrt(quantile(RowVar, probs=sd.quantile.for.logFC))
  names(prior.logFC) <- NULL

  # Fit an additive model for each protein and save residual variance
  ProtVar <- proteinResVarFromCompletePeptideData(yComplete, protein.id=protein.id, reorder=FALSE)
 
  # Posterior protein variances
  ModProtVar <- squeezeVar(ProtVar$s2, ProtVar$df.residual, robust=robust)

  list(prior.mean = prior.mean, 
       prior.sd = prior.sd, 
       prior.logFC = prior.logFC, 
       sigma = sqrt(ModProtVar$var.post))
}