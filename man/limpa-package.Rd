\name{limpa-package}
\docType{package}
\alias{limpa}
\alias{limpa-package}
\title{Linear Models for Proteomics Data (Accounting for Missing Values)}

\description{
This package implements a pipeline for quantification and differential expression analysis of mass-spectrometry-based proteomics data.

Mass spectrometry (MS)-based proteomics is a powerful tool in biomedical research. 
Recent advancements in label-free methods and MS instruments have enabled the quantitative characterisation of large-scale complex biological samples with the increasingly deeper coverage of the proteome. 
However, missing values are still ubiquitous in MS-based proteomics data. 
We observe from a wide range of real datasets that missingness in label-free data is intensity-dependent, so that the missing values are missing not at random or, in other words, are non-ignorable. 

This package implements statistical and computational methods for analysing MS-based label-free proteomics data with non-ignorable missing values. 
The package use the observed proteomics data to estimate the detection probability curve (DPC), which provides a formal probabilistic model for the intensity-dependent missingness.
Based on exponential tilting, the DPC estimates the detection probabilities given the underlying intensity of each observation, observed or unobserved. 
Importantly, the DPC evaluates how much statistical information can or cannot be recovered from the missing value pattern, and can be used to inform downstream analyses such as differential expression (DE) analysis. 

Next, the package implements a novel protein quantification method, called DPC-quant, where missing values are represented by the DPC. 
An empirical Bayes scheme is employed to borrow information across the tens of thousands of peptides measured in a typical experiment. 
A multivariate normal prior is estimated empirically from data to describe the variability in log-intensities across the samples and across the peptides.

Finally, quantification uncertainty is incorporated into the differential expression analysis using precision weights.
Leveraging the limma package, a new variance modelling approach with multiple predictors is used, which allows the DPC-quant precisions to be propagated to the differential expression analysis while simultaneously assuming a mean-variance relationship.
The new differential expression pipeline has been implemented in the limma R package in the \code{vooma()} function. 

The limpa package is fully compatible with limma pipelines, allowing any arbitrarily complex experimental design and other downstream tasks such as the gene ontology or pathway analysis. 
}

\author{
Mengbo Li and Gordon K Smyth
}

\references{
Li M (2024).
Linear Models and Empirical Bayes Methods for Mass Spectrometry-based Proteomics Data.
PhD Thesis, University of Melbourne.
\url{http://hdl.handle.net/11343/351600}

Li M, Smyth GK (2023).
Neither random nor censored: estimating intensity-dependent probabilities for missing values in label-free proteomics.
\emph{Bioinformatics} 39(5), btad200.
\doi{10.1093/bioinformatics/btad200}

Li M, Cobbold SA, Smyth GK (2025).
Quantification and differential analysis of mass spectrometry proteomics data with probabilistic recovery of information from missing values.
\emph{bioRxiv} 2025/651125.
\doi{10.1101/2025.04.28.651125}
}

\concept{Documentation}
