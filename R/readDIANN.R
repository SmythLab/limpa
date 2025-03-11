readDIANN <- function(file="Report.tsv", path=NULL, format = "tsv", sep="\t", log=TRUE, q.columns = c("Global.Q.Value", "Lib.Q.Value"), q.cutoffs = c(0.01, 0.01))
  # Read Report.tsv from DIA-NN output
  # Gordon Smyth and Mengbo Li
  # Created 3 July 2023. Last modified 3 March 2025.
{
  # Read DIA-NN report file
  if (!is.null(path)) file <- file.path(path, file)
  Select <- c("Run", "Protein.Group", "Protein.Names", "Genes", "Precursor.Id", "Proteotypic", "Precursor.Normalised", q.columns)
  format <- match.arg(format, choices = c("tsv", "parquet"))

  if (identical(format, "tsv")) {
    Report <- fread(file, sep = "\t", select = Select)
  }
  
  #	Use arrow package for reading
  if (identical(format, "parquet")) {
    suppressPackageStartupMessages(OK <- requireNamespace("arrow",quietly = TRUE))
	  if(!OK) stop("arrow package required but is not installed (or can't be loaded)")
    Report <- arrow::read_parquet(file)
    Report <- as.data.frame(Report[, Select])
  }
  
  # Filter by q-values
  if (length(q.columns) > 0L) {
    if (!identical(length(q.cutoffs), length(q.columns))) {
      q.cutoffs <- rep_len(q.cutoffs[1], length(q.columns))
      message("Length of q-value columns does not match with length of q-value cutoffs. Use q.cutoffs[1] for all columns.")
    }
    kp <- rep_len(TRUE, nrow(Report))
    for (qcol in seq_along(q.columns)) {
      kp <- kp & Report[[q.columns[qcol]]] <= q.cutoffs[qcol]
    }
    Report <- Report[kp, ]
  }
  
  # Convert intensities to wide format
  Samples <- unique(Report$Run)
  Precursors <- unique(Report$Precursor.Id)
  y <- matrix(0, length(Precursors), length(Samples))
  mSample <- match(Report$Run, Samples)
  mPrecursor <- match(Report$Precursor.Id, Precursors)
  i <- mPrecursor + (mSample - 1L) * length(Precursors)
  y[i] <- Report$Precursor.Normalised
  colnames(y) <- Samples
  rownames(y) <- Precursors
  
  # Precursor annotation in wide format
  d <- duplicated(Report$Precursor.Id)
  Genes <- data.frame(Report[!d, 2:6])
  row.names(Genes) <- Precursors
  
  # Output either unlogged EListRaw (with zeros) or logged Elist (with NAs)
  if (log) {
    y[y < 1e-8] <- NA
    y <- log2(y)
    new("EList", list(E = y, genes = Genes))
  } else {
    new("EListRaw", list(E = y, genes = Genes))
  }
}
