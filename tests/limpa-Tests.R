library(limpa)
options(warnPartialMatchArgs=TRUE,warnPartialMatchAttr=TRUE,warnPartialMatchDollar=TRUE,width=120)

set.seed(0); u <- runif(100)

y.peptide <- simProteinDataSet()
names(y.peptide)
colSums(is.na(y.peptide$E))
summary(y.peptide$E[,1])
names(y.peptide)

dpcfit <- dpc(y.peptide)
dpcfit$dpc

y.protein <- dpcQuant(y.peptide,protein.id="Protein",dpcfit)
summary(y.protein$E[,1])

Group <- factor(y.peptide$targets$Group)
design <- model.matrix(~Group)
fit <- dpcDE(y.protein, design, plot=FALSE)
summary(fit$coefficients)
summary(fit$sigma)
