readSpectronaut <- function(
  file="Report.tsv", path=NULL, sep="\t", log=TRUE,
  run.column = "R.Raw File Name",
  precursor.column = "EG.PrecursorId",
  qty.column = "EG.TotalQuantity (Settings)",
  q.columns = c("EG.Qvalue", "PG.Qvalue"), q.cutoffs = 0.01,
  extra.columns = c("PG.ProteinAccessions", "EG.IsImputed")
)
# Read normal table output from Spectronaut.
# Gordon Smyth and Mengbo Li
# Created 18 December 2023. Last modified 7 April 2025.
{
  # Read Spectronaut report file
  if (!is.null(path)) file <- file.path(path, file)

  # Read column names
  all.columns <- fread(file, sep = "\t", nrows = 0L)
  all.columns <- colnames(all.columns)

  Select <- c(run.column, precursor.column, qty.column, q.columns, extra.columns)
  if (any(!(Select %in% all.columns))) {
    no.in.Select <- setdiff(Select, all.columns)
    message(paste("Columns", paste(no.in.Select, collapse = ","), "not in data!", sep = " "))
    message("Reading the rest of the columns only.")
  }
  Select <- intersect(Select, all.columns)
  extra.columns <- intersect(extra.columns, all.columns)
  Report <- fread(file, sep = "\t", select = Select)
  colnames(Report)[which(colnames(Report) == run.column)] <- "Run"
  colnames(Report)[which(colnames(Report) == precursor.column)] <- "Precursor.Id"
  colnames(Report)[which(colnames(Report) == qty.column)] <- "Intensity"

  # Filter by imputed
  if ("EG.IsImputed" %in% colnames(Report)) {
    message("Filtering out imputed values according to `EG.IsImputed`.")
    Report <- Report[Report$EG.IsImputed == FALSE, ]
  }

  # Filter by q-values
  if (length(q.columns) > 0L) {
    q.columns <- q.columns[q.columns %in% colnames(Report)]
    if (!identical(length(q.cutoffs), length(q.columns))) {
      q.cutoffs <- rep(q.cutoffs[1], length(q.columns))
      message("Length of q-value columns does not match with length of q-value cutoffs.
      Use q.cutoffs[1] for all columns.")
    }
    kp <- rep(TRUE, nrow(Report))
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
  y[i] <- Report$Intensity
  colnames(y) <- Samples
  rownames(y) <- Precursors

  # Precursor annotation in wide format
  d <- duplicated(Report$Precursor.Id)
  Genes <- data.frame(Report[!d, ])[, setdiff(extra.columns, "EG.IsImputed")]
  row.names(Genes) <- Precursors

  # Output either unlogged EListRaw (with zeros) or logged Elist (with NAs)
  if (log) {
    y[y <= 1L] <- NA
  # Remove rows that are missing in all samples
    ind <- rowMeans(is.na(y)) < 1
    y <- y[ind, , drop = FALSE]
    Genes <- Genes[ind, ]
  # Log
    y <- log2(y)
    new("EList", list(E = y, genes = Genes))
  } else {
    new("EListRaw", list(E = y, genes = Genes))
  }
}
