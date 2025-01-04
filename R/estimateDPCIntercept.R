estimateDPCIntercept <- function(y, dpc.slope=0.8, trace=FALSE)
# For a preset DPC slope, estimate the intercept that gives the correct proportion of missing values overall.
# Created 2 Jan 2025. Last modified 4 Jan 2025.
{
  y <- as.matrix(y)
  IsObs <- as.integer(!is.na(y))

# If dpc.slope is zero, return overall proportion of detection
  if(abs(dpc.slope) < 1e-8) {
    beta0 <- qlogis(mean(IsObs))
    names(beta0) <- "Intercept"
    return(beta0)
  }

# Impute to get putative complete matrix
  y <- as.vector(imputeByExpTilt(y, dpc.slope=dpc.slope))

# Estimate intercept of logistic regression
# If more than 10k observations, then aggregate y into equal intervals
  if(length(y) > 10000L) {
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
    fit <- glm(PropObs~offset(o), family=binomial(), weights=N, trace=trace)
  } else {
    o <- dpc.slope*y
    fit <- glm(IsObs~offset(o), family=binomial(), trace=trace)
  }

# Output
  beta0 <- coef(fit)[1]
  names(beta0) <- "Intercept"
  beta0
}
