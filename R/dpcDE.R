dpcDE <- function(y, design, plot=TRUE, ...)
# Using `y` output by dpcQuant(), calls voomLmFit() read for DE analysis.
# Created 27 Dec 2024.
{
  if(is.null(y$other$standard.error)) stop("standard errors not found. y should be an EList produced by dpcQuant().")
  voomaLmFit(y=y, design=design, predictor=log(y$other$standard.error), plot=plot, ...)
}
