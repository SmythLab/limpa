estimateDPCIntercept <- function(y, dpc.slope=0.8, trace=FALSE)
# For a preset DPC slope, estimate the intercept that gives the correct proportion of missing values overall.
# Created 2 Jan 2025.
{
  y <- as.matrix(y)
  IsObs <- as.integer(!is.na(y))
  y <- as.vector(imputeByExpTilt(y, dpc.slope=dpc.slope))
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
  beta0 <- coef(fit)[1]
  names(beta0) <- "Intercept"
  beta0
}
