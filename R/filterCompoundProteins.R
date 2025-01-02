filterCompoundProteins <- function(y, ...)
UseMethod("filterCompoundProteins")

filterCompoundProteins.default <- function(y, protein.group, ...)
# Remove compound proteins from matrix
# Created 24 July 2023. Last modified 10 August 2023.
{
# Check input
  if(anyNA(protein.group)) stop("NAs not allowed in protein.group")

# Remove rows with ";" join in protein.group
  i <- grep(";",protein.group)
  if (identical(length(i), 0L)) {
     return(y)
  } else {
     return(y[-i,,drop=FALSE])
  }

 }

filterCompoundProteins.EList <- filterCompoundProteins.EListRaw <- function(y, protein.group="Protein.Group", ...)
# Remove compound proteins from EList
# Created 24 July 2023. Last modified 27 December 2024.
{
# Get protein.group vector
  protein.group <- as.character(protein.group)
  if(identical(length(protein.group),1L)) {
    ColName <- protein.group
    protein.group <- y$genes[[protein.group]]
    if(is.null(protein.group)) stop("Column \"",ColName,"\" not found in y$genes")
  } else {
    if(!identical(nrow(y),length(protein.group))) stop("length(protein.group) must match nrows(y)")
  }

# Remove rows with ";" join in protein.group
  i <- grep(";",protein.group)
  if (identical(length(i), 0L)) {
     return(y)
  } else {
     return(y[-i,,drop=FALSE])
  }

}
