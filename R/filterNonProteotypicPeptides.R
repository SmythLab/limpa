filterNonProteotypicPeptides <- function(y, ...)
  UseMethod("filterNonProteotypicPeptides")

filterNonProteotypicPeptides.default <- function(y, proteotypic, ...)
# Remove peptides that are not proteotypic.
# Default method for matrices.
# Created 11 Sep 2023. Last modified 26 Dec 2024.
{
  # Check for NAs
  if(anyNA(proteotypic)) stop("Need proteotypic information for each peptide!")

  # proteotypic is expected to be a 0/1 vector or a TRUE/FALSE vector.
  if(is.character(proteotypic)) proteotypic <- as.integer(proteotypic)
  if(is.factor(proteotypic)) proteotypic <- as.integer(as.character(proteotypic))
  if(length(unique(proteotypic)) > 2L) stop("proteotypic should contain only 0/1 or TRUE/FALSE values")

  # Filter
  y[as.logical(proteotypic),,drop=FALSE]
}

filterNonProteotypicPeptides.EList <- filterNonProteotypicPeptides.EListRaw <- function(y, proteotypic="Proteotypic", ...)
# Remove peptides that are not proteotypic.
# Method for EList and EListRaw objects.
# Created 11 Sep 2023. Last modified 26 Dec 2024.
{
  # Get proteotypic vector
  proteotypic <- as.character(proteotypic)
  if(identical(length(proteotypic),1L)) {
    ColName <- proteotypic
    proteotypic <- y$genes[[ColName]]
    if(is.null(proteotypic)) stop("Column \"",ColName,"\" not found in y$genes")
  } else {
    if(!identical(nrow(y),length(proteotypic))) stop("length(proteotypic) must match nrows(y) or be of length 1 as the column name for proteotypicity!")
    y$genes$proteotypic <- proteotypic
  }

  # Check for NAs
  if(anyNA(proteotypic)) stop("Need proteotypic information for each peptide!")

  # proteotypic is expected to be a 0/1 vector or a TRUE/FALSE vector.
  if(is.character(proteotypic)) proteotypic <- as.integer(proteotypic)
  if(is.factor(proteotypic)) proteotypic <- as.integer(as.character(proteotypic))
  if(length(unique(proteotypic)) > 2L) stop("proteotypic should contain only 0/1 or TRUE/FALSE values")

  # Filter
  y[as.logical(proteotypic),,drop=FALSE]
}
