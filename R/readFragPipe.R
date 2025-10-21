readFragPipe <- function(
    file = "combined_ion.tsv", path = NULL, sep = "\t",
    log = TRUE,
    peptide.column = c("Modified Sequence", "Charge"),
    qty.column = NULL,
    qty.column.key = " Intensity",
    extra.columns = c("Protein", "Protein ID", "Gene", "Protein Description", "Mapped Proteins"),
    match.type.key = NULL,
    maxlfq = FALSE
)
# Read combined_ion.tsv, combined_peptide.tsv, or combined_modified_peptide.tsv file from FragPipe output
# Mengbo Li, Pedro Baldoni and Gordon Smyth
# Created 8 Oct 2025. Last modified 15 Oct 2025.
{
  # Read file
  if (!is.null(path)) file <- file.path(path, file)

  # Read column names
  all.columns <- fread(file, sep = sep, nrows = 0L, showProgress = FALSE)
  all.columns <- colnames(all.columns)

  # get qty.columns by key
  if (!is.null(qty.column.key)) {
    message("Searching for column names matching to `qty.column.key`.")
    qty.column <- all.columns[grepl(qty.column.key, all.columns)]
  }
  # get rid of MaxLFQ columns
  if (!maxlfq) {
    qty.column <- qty.column[setdiff(seq_len(length(qty.column)), grep("MaxLFQ", qty.column))]
  }

  Select <- c(peptide.column, qty.column, extra.columns)
  if (any(!(Select %in% all.columns))) {
    no.in.Select <- setdiff(Select, all.columns)
    message(paste("Columns", paste(no.in.Select, collapse = ","), "not in data!", sep = " "))
    message("Reading the rest of the columns only.")
  }
  Select <- intersect(Select, all.columns)
  extra.columns <- intersect(extra.columns, all.columns)

  # Read in all required columns
  Report <- fread(file, sep = sep, select = Select,
  data.table = FALSE, showProgress = FALSE)

  # if length of peptide.column > 1, concatenate them
  if(length(peptide.column) > 1L) {
    message("Multiple columns are input to `peptide.column`, concatenating them to make a unique precursor Id.")
    Report$Precursor.Id <- vapply(seq_len(nrow(Report)),
    function(i) paste(unlist(Report[i, peptide.column, drop = FALSE]), collapse = "_"),
    character(1))
    if(any(duplicated(Report$Precursor.Id))) {
      stop("Duplicated ion ids from `peptide.column`, check input!")
    } else {
      peptide.column <- "Precursor.Id"
      Report <- Report[, c(peptide.column, qty.column, extra.columns), drop = FALSE]
    }
  } else {
    if (any(duplicated(unlist(Report[, peptide.column, drop = FALSE])))) {
      stop("Duplicated ion ids from `peptide.column`, check input!")
    }
  }

  # Get the intensity matrix
  y <- Report[, qty.column, drop = FALSE]
  row.names(y) <- Report[,peptide.column]
  if (!is.null(qty.column.key)) {
    colnames(y) <- gsub(qty.column.key, "", colnames(y), fixed = TRUE)
  }

  # Precursor annotation in wide format
  Genes <- Report[, extra.columns, drop = FALSE]
  colnames(Genes) <- extra.columns
  row.names(Genes) <- Report[,peptide.column]

  # Grab MBR annotation if requested
  if (!is.null(match.type.key)) {
    match.type.cols <- all.columns[grepl(match.type.key, all.columns)]
    match.type <- fread(file, sep = sep, select = match.type.cols,
    data.table = TRUE, showProgress = FALSE)
    match.type <- as.data.frame(match.type)
    colnames(match.type) <- gsub(match.type.key, "", colnames(match.type), fixed = TRUE)
    row.names(match.type) <- Report$peptide.column
  }

  # Output either unlogged EListRaw (with zeros) or logged Elist (with NAs)
  if (log) {
    y[y == 0L] <- NA
    # Remove rows that are missing in all samples
    ind <- rowMeans(is.na(y)) < 1
    y <- y[ind, , drop = FALSE]
    Genes <- Genes[ind, ]
    # Log
    y <- log2(y)
    if (!is.null(match.type.key)) {
      out <- new("EList", list(E = y, genes = Genes))
      out$other$match.type <- match.type[ind, , drop = FALSE]
    } else {
      out <- new("EList", list(E = y, genes = Genes))
    }
  } else {
    if (!is.null(match.type.key)) {
      out <- new("EListRaw", list(E = y, genes = Genes))
      out$other$match.type <- match.type[ind, , drop = FALSE]
    } else {
      out <- new("EListRaw", list(E = y, genes = Genes))
    }
  }
  return(out)
}
