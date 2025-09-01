readDIANN <- function(
  file="Report.tsv", path=NULL, format = "tsv", sep="\t", log=TRUE, 
  run.column = "Run",
  precursor.column = "Precursor.Id",
  qty.column = "Precursor.Normalised",
  q.columns = c("Global.Q.Value", "Lib.Q.Value"), q.cutoffs = c(0.01, 0.01), 
  extra.columns = c("Protein.Group", "Protein.Names", "Genes", "Proteotypic")
  )
# Read Report.tsv from DIA-NN output
# Gordon Smyth and Mengbo Li
# Created 3 July 2023. Last modified 2 September 2025.
{
  # Check arguments
  if (!is.null(path)) file <- file.path(path, file)
  Select <- c(run.column, precursor.column, qty.column, q.columns, extra.columns)
  format <- match.arg(format, choices = c("tsv", "parquet"))

  # Read DIA-NN report file
  if (identical(format, "tsv")) {
    Report <- fread(file, sep = "\t", select = Select)
  } else {
  #	Use arrow package to read Parquet format file
    suppressPackageStartupMessages(OK <- requireNamespace("arrow",quietly = TRUE))
	  if(!OK) stop("arrow package required but is not installed (or can't be loaded)")
    Report <- arrow::read_parquet(file)
    Report <- as.data.frame(Report[, Select])
  }

  all.columns <- colnames(Report)
  if (any(!(Select %in% all.columns))) {
    no.in.Select <- setdiff(Select, all.columns)
    message(paste("Columns", paste(no.in.Select, collapse = ","), "not in data!", sep = " "))
    message("Reading the rest of the columns only.")
  }
  Select <- intersect(Select, all.columns)
  extra.columns <- intersect(extra.columns, all.columns)
  colnames(Report)[which(colnames(Report) == run.column)] <- "Run"
  colnames(Report)[which(colnames(Report) == precursor.column)] <- "Precursor.Id"
  colnames(Report)[which(colnames(Report) == qty.column)] <- "Intensity"
  
  # Filter by q-values
  if (length(q.columns) > 0L) {
    if (!identical(length(q.cutoffs), length(q.columns))) {
      q.cutoffs <- rep_len(q.cutoffs[1], length(q.columns))
      message("Length of q-value columns does not match with length of q-value cutoffs. Use q.cutoffs[1] for all columns.")
    }
    kp <- rep_len(TRUE, nrow(Report))
    for (qcol in seq_along(q.columns)) {
      kp[ Report[[ q.columns[qcol] ]] > q.cutoffs[qcol] ] <- FALSE
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
  y[i] <- Report$Intensity
  colnames(y) <- Samples
  rownames(y) <- Precursors
  
  # Precursor annotation in wide format
  d <- duplicated(Report$Precursor.Id)
  Genes <- data.frame(Report[!d, ..extra.columns])
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
