\name{limpanews}
\title{limpa News}
\encoding{UTF-8}


\section{Version 0.99.12 (2025-04-07)}{\itemize{
\item
dpcQuant() and peptides2Proteins() now return essentially equal values if there are 
multiple entirely imputed values for the same protein.

\item
Update readSpectronaut() to allow more flexibility regarding which columns are read.

\item
Bug fix to voomaLmFitWithImputation() for proteins
that are almost entirely imputed and have no remaining
degrees of freedom for estimating the residual variance.
Previously the degrees of freedom were not always being recognized as zero
because of floating point errors.

\item
New calculation of lost degrees of freedom in voomaLmFitWithImputation().
The new method always returns integer residual degrees of freedom and no longer assumes that entirely imputed values are all equal for the same gene.
}}


\section{Version 0.99.11 (2025-04-05)}{\itemize{
\item
Remove reference to DIA in vignette, because limpa works for both DIA and DDA.
}}


\section{Version 0.99.10 (2025-03-25)}{\itemize{
\item
Use new function plotMDSUsingSEs() in vignette.

\item
Remove Author and Maintainer fields from DESCRIPTION again.
}}


\section{Version 0.99.9 (2025-03-23)}{\itemize{
\item
New function plotMDSUsingSEs().
}}


\section{Version 0.99.8 (2025-03-22)}{\itemize{
\item
Add Mengbo Li as joint maintainer.
Add Author and Maintainer fields to DESCRIPTION.
}}


\section{Version 0.99.7 (2025-03-21)}{\itemize{
\item
New argument `chunk` for dpcQuant() and dpcImpute().

\item
dpcImpute() now passes the `verbose` argument to peptides2Proteins().

\item
Add more background information to the vignette, including references for DIA-NN, Spectronaut and MaxQuant.
Add mention of voomaLmFitWithImputation().

\item
Bug fix to starting values in peptides2Proteins() for single-peptide proteins.
}}


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
