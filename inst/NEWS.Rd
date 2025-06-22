\name{limpanews}
\title{limpa News}
\encoding{UTF-8}


\section{Version 1.0.5 (2025-06-22)}{\itemize{
\item
New function plotPeptides() to plot the peptide-level log-intensities for one protein.

\item
New argument `logit.p` for dztbinom(), which allows the function to avoid subtractive cancellation for success probabilities close to 1.
dpc() now uses the new argument when calling dztbinom().

\item
Set the default value of `b1.upper` back to 1 in dpc().

\item
Revisions to dpcCN() to improve convergence.
To improve speed, the function now randomly selects a subset of rows of data to process, with the number of rows given by the new argument `subset`.
The default number of outer iterations is decreased to 2.

\item
Add sections on dpcImpute() and imputeByExpTilt() to the package vignette.
}}


\section{Version 1.0.4 (2025-06-19)}{\itemize{
\item
dpcQuant() now calls dpcImpute() automatically if all proteins have just one peptide.
}}


\section{Version 1.0.3 (2025-05-26)}{\itemize{
\item
Fix simProteinDataSet() processing of prop.missing.

\item
Increase default value of `b1.upper` to 2 in dpc().

\item
Set row.names in dpcImpute().
}}


\section{Version 1.0.2 (2025-04-30)}{\itemize{
\item
Cite Li et al (bioRxiv 2025)
}}


\section{Version 1.0.1 (2025-04-20)}{\itemize{
\item
Expand Quantification section in the vignette to give more math details.
}}


\section{Version 1.0.0 (2025-04-16)}{\itemize{
\item
dpcQuant() and peptides2Proteins() now return essentially equal values if there are 
multiple entirely imputed values for the same protein.

\item
Update readSpectronaut() to allow more flexibility regarding which columns are read.

\item
Remove reference to DIA in vignette, because limpa works for both DIA and DDA.

\item
New function plotMDSUsingSEs(), which accounts for standard errors when computing between-sample distances for the plot.

\item
New argument `chunk` for dpcQuant() and dpcImpute().

\item
dpcImpute() now passes the `verbose` argument to peptides2Proteins().

\item
Add more background information to the vignette, including references for DIA-NN, Spectronaut and MaxQuant.
Add mention of voomaLmFitWithImputation().

\item
Bug fix to starting values in peptides2Proteins() for single-peptide proteins.

\item
expTiltByColumns() is now more robust for datasets with very few rows, especially when the column variance is exactly zero by chance.
It originally assumed that the number of rows would always be large.
This also improves the robustness of imputeByExpTilt().

\item
peptides2Proteins() now runs imputeByExpTilt() on the whole dataset in order to get starting values for each protein.
This improves speed slightly and makes the starting values more robust for large datasets with many NAs.

\item
readDIANN() now supports DIA-NN output in either tsv or Parquet formats.

\item
New functions dpcImpute() and dpcImputeHyperparam() to facilitate imputation and differential expression analysis at the peptide level.
This may be needed for phosphoproteomics and analogous applications.

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
It now estimates the jittering amount automatically.

\item
Add readDIANN() code to vignette.

\item
Add peptide filtering to vignette.

\item
New function voomaLmFitWithImputation(), which is the same as limma::voomaLmFit() except that it allows for loss of degrees of freedom when all the values for a particular protein in a particular group are imputed.
dpcDE() now calls this new function instead of limma::voomaLmFit(), making the differential expression analysis more robust to high rates of missing values.

\item
Fix dpcQuant() `genes` output when the input `y` is a plain matrix instead of an `EList` object.

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
