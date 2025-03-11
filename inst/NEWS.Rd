\name{limpanews}
\title{limpa News}
\encoding{UTF-8}

 
\section{Version 0.99.6 (2025-03-11)}{\itemize{
\item
expTiltByColumns() is now more robust for datasets with very few rows, especially when the column variance is exactly zero by chance.
It originally assumed that the number of rows would always be large.
This also improves the robustness of imputeByExpTilt().

\item
peptides2Proteins() now runs imputeByExpTilt() on the whole dataset in order to get starting values for each protein.
This improves speed slightly and makes the starting values more robust for large datasets with many NAs.

\item
readDIANN() now supports DIA-NN output in either tsv or Parquet formats.
}}


\section{Version 0.99.5 (2025-02-25)}{\itemize{
\item
New functions dpcImpute() and dpcImputeHyperparam() to facilitate imputation and differential expression analysis at the peptide level.
This may be needed for phosphoproteomics and analogous applications.
}}


\section{Version 0.99.4 (2025-02-15)}{\itemize{
\item
Package title changed to "Quantification and Differential Analysis of Proteomics Data", and package description expanded.

\item
New function plotProtein(), to plot the quantifications for a given protein with error bars.

\item
Cleanup annotation output from dpcQuant() to remove peptide-level annotation columns

\item
dpc() now outputs `n.detected` for each peptide.

\item
dpc() now checks for peptides that are entirely missing.

\item
plotDPC() now accepts input from either dpc() or dpcCN().
It now estimates jittering amount automatically.

\item
Add readDIANN() code to vignette.
}}


\section{Version 0.99.3 (2025-02-15)}{\itemize{
\item
Add peptide filtering to vignette.

\item
New function voomaLmFitWithImputation(), which is the same as limma::voomaLmFit() except that it allows for loss of degrees of freedom when all the values for a particular protein in a particular group are imputed and identical.
dpcDE() now calls thew new function instead of limma::voomaLmFit(), making the differential expression analysis more robust to high rates of missing values.
}}


\section{Version 0.99.2 (2025-01-13)}{\itemize{
\item
Fix dpcQuant() `genes` output when the input `y` is a plain matrix instead of an `EList` object.
}}


\section{Version 0.99.1 (2025-01-04)}{\itemize{
\item
New argument `dpc.slope` for dpcQuant(), to allow the use of a preset DPC slope.

\item
estimateDPCIntercept() now allows `dpc.slope=0`, corresponding to the case that missing values are uninformative.
}}


\section{Version 0.99.0 (2025-01-03)}{\itemize{
\item
Initial public release of the limpa package ("Linear Models for Proteomics Data")
at \url{https://github.com/SmythLab/limpa} with the following exported functions:
completeMomentsON, dpc, dpcCN, dpcDE, dpcQuant, dpcQuant.default,
dpcQuant.EList, dpcQuantHyperparam, dztbinom, estimateDPCIntercept, expTiltByColumns,
expTiltByRows, filterCompoundProteins, filterCompoundProteins.default,
filterCompoundProteins.EList, filterCompoundProteins.EListRaw,
filterNonProteotypicPeptides, filterNonProteotypicPeptides.default,
filterNonProteotypicPeptides.EList, filterNonProteotypicPeptides.EListRaw,
filterSingletonPeptides, filterSingletonPeptides.default,
filterSingletonPeptides.EList, filterSingletonPeptides.EListRaw,
fitZTLogit, imputeByExpTilt, imputeByExpTilt.default,
imputeByExpTilt.EList, imputeByExpTilt.EListRaw, observedMomentsCN,
peptides2ProteinBFGS, peptides2ProteinNewton, peptides2Proteins,
peptides2ProteinWithoutNAs, plotDPC, proteinResVarFromCompletePeptideData,
pztbinom, readDIANN, readSpectronaut, removeNARows, removeNARows.default,
removeNARows.EList, simCompleteDataCN, simCompleteDataON, simProteinDataSet.
}}
