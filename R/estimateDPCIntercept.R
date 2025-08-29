estimateDPCIntercept <- function(y, dpc.slope=0.8, verbose=FALSE)
# For a preset DPC slope, estimate the intercept that gives the correct proportion of missing values overall.
# Created 2 Jan 2025. Last modified 28 Jun 2025.
{
  y <- as.matrix(y)
  IsObs <- as.integer(!is.na(y))

# If dpc.slope is zero, return overall proportion of detection
  if(abs(dpc.slope) < 1e-8) {
    return(qlogis(mean(IsObs)))
  }

# Impute to get putative complete matrix
  if(verbose) message("Imputing to get putative complete matrix ...")
  y <- as.vector(imputeByExpTilt(y, dpc.slope=dpc.slope))

# Estimate intercept of logistic regression
# If more than 10k observations, then aggregate y into equal intervals
  if(length(y) > 10000L) {
    if(verbose) message("Aggregating to reduce number of observations ...")
    cuty <- as.integer(cut(y,300))
    N <- rowsum(rep_len(1L,length(y)),cuty)
    if(identical(min(N),0L)) {
      cuty <- cuty[N > 0L]
      N <- N[N > 0L]
    }
    MeanY <- rowsum(y,cuty) / N
    NObs <- rowsum(IsObs,cuty)
    PropObs <- NObs / N
    o <- dpc.slope*MeanY
    X <- matrix(1,length(NObs),1)
    if(verbose) message("Running glm.fit ...")
    fit <- suppressWarnings(glm.fit(X,PropObs,offset=o,family=binomial(),weights=N,control=glm.control(trace=verbose),intercept=FALSE))
  } else {
    o <- dpc.slope*y
    n <- length(IsObs)
    X <- matrix(1,n,1)
    N <- rep_len(1,n)
    if(verbose) message("Running glm.fit ...")
    fit <- suppressWarnings(glm.fit(X,IsObs,offset=o,family=binomial(),weights=N,control=glm.control(trace=verbose),intercept=FALSE))
  }

# Output
  if(verbose) message("Intercept = ",fit$coef)
  fit$coef
}
