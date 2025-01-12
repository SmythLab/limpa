dpcDE <- function(y, design, plot=TRUE, ...)
# Using `y` output by dpcQuant(), calls voomLmFit() read for DE analysis.
# Created 27 Dec 2024. Last modified 13 Jan 2025.
{
  if(is.null(y$other$standard.error)) stop("standard errors not found. y should be an EList produced by dpcQuant().")
  voomaLmFitWithImputation(y=y, imputed=!y$other$n.observations, design=design, predictor=log(y$other$standard.error), plot=plot, ...)
}
