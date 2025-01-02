removeNARows <- function(y, ...)
  UseMethod("removeNARows")

removeNARows.default <- function(y, nobs.min=1, ...)
# Remove entirely missing rows from a matrix
# Created 23 Dec 2024. Last modified 31 Dec 2024.
{
  y <- as.matrix(y)
  nobs <- ncol(y) - rowSums(is.na(y))
  if(min(nobs) < nobs.min) y <- y[nobs >= nobs.min,,drop=FALSE]
  y
}

removeNARows.EList <- function(y, nobs.min=1, ...)
# Remove entirely missing rows from an EList object.
# Created 31 Dec 2024. Last modified 31 Dec 2024.
{
  nobs <- ncol(y) - rowSums(is.na(y$E))
  if(min(nobs) < nobs.min) y <- y[nobs >= nobs.min,,]
  y
}
