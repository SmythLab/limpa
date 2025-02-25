dpcImputeHyperparam <- function(y, dpc.slope=0.7, sd.quantile.for.logFC = 0.9, robust=FALSE, ...)
# Estimate hyperparameters needed for dpcImpute
# Mengbo Li and Gordon Smyth
# Created 25 Feb 2025. Last modified 29 Dec 2024.
{
  # Check input
  y <- as.matrix(y)
  nsamples <- ncol(y)

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

  # Rowwise residual variance
#  nsv <- round(nsamples^1/3)
#  SVD <- svd(yComplete,nu=0,nv=nsv)
#  fit <- lm.fit(SVD$v,t(yComplete))
#  s2 <- colMeans(fit$effects[-(1:nsv),]^2)
 
  # Posterior protein variances
#  s2.post <- squeezeVar(s2, nsamples-nsv, robust=robust)$var.post

  list(prior.mean = prior.mean, 
       prior.sd = prior.sd, 
       prior.logFC = prior.logFC)
}