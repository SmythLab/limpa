dpcQuant <- function(y, ...)
  UseMethod("dpcQuant")

dpcQuant.default <- function(y, protein.id, dpc=NULL, dpc.slope=0.8, verbose=TRUE, ...)
# Use the DPC to quantify protein expression values by maximum posterior.
# Created 31 Dec 2024. Last modified 4 Jan 2025.
{
# Check y
  y <- as.matrix(y)

# Check protein.id
  if(missing(protein.id)) stop("Need protein ID")
  if(!identical(nrow(y),length(protein.id))) stop("length(protein.id) must match nrows(y)")

# Construct EList and pass to EList method
  z <- new("EList",list(E=y,genes=data.frame(Protein=protein.id)))
  dpcQuant(z,protein.id=protein.id,dpc=dpc,dpc.slope=dpc.slope,verbose=verbose,...)
}

dpcQuant.EList <- function(y, protein.id="Protein.Group", dpc=NULL, dpc.slope=0.8, verbose=TRUE, ...)
# Use the DPC to quantify protein expression values by maximum posterior.
# Created 27 Dec 2024. Last modified 4 Jan 2024.
{
# Check dpc
  if(is.list(dpc)) dpc <- dpc$dpc
  if(is.null(dpc)) {
    beta0 <- estimateDPCIntercept(y, dpc.slope=dpc.slope)
    if(verbose) message("DPC intercept estimated as ", formatC(beta0))
    dpc <- c(beta0=beta0, beta1=dpc.slope)
  }

# Get vector of protein IDs
# protein.id can be either an annotation column name or a vector of IDs.
  protein.id <- as.character(protein.id)
  if(identical(length(protein.id),1L)) {
    ColName <- protein.id
    protein.id <- y$genes[[ColName]]
    if(is.null(protein.id)) stop("Column \"",ColName,"\" not found in y$genes")
  } else {
    if(!identical(nrow(y),length(protein.id))) stop("length(protein.id) must match nrows(y)")
  }

# Sort peptides in protein order
  o <- order(protein.id)
  protein.id <- protein.id[o]
  y <- y[o,]

# Estimate Bayes hyperparameters
  if(verbose) message("Estimating hyperparameters ...")
  h <- dpcQuantHyperparam(y, protein.id=protein.id, dpc.slope=dpc[2], ...)

# Summarize peptides to proteins by maximum posterior estimation
  if(verbose) message("Quantifying proteins ...")
  y.protein <- peptides2Proteins(y,
                          protein.id = protein.id,
                          dpc = dpc,
                          sigma = h$sigma,
                          prior.mean = h$prior.mean,
                          prior.sd = h$prior.sd,
                          prior.logFC = h$prior.logFC,
                          standard.errors = TRUE,
                          verbose = verbose)

# Add back original original annotation
  d <- !duplicated(protein.id)
  genes <- y$genes[d,]
  row.names(genes) <- row.names(y.protein)
  y.protein$genes <- data.frame(genes,y.protein$genes)
  y.protein$targets <- y$targets
  y.protein$dpc <- dpc
  y.protein
}
