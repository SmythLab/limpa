dpcDE <- function(y, design, plot=TRUE, ...)
# Using `y` output by dpcQuant(), calls voomLmFit() read for DE analysis.
# Created 27 Dec 2024. Last modified 6 Apr 2025.
{
  if(is.null(y$other$standard.error)) stop("standard errors not found. y should be an EList produced by dpcQuant().")
# Arbitrary small value to avoid taking log of zero
  eps <- 1e-6
  voomaLmFitWithImputation(y=y, imputed=!y$other$n.observations, design=design, predictor=log(y$other$standard.error+eps), plot=plot, ...)
}
