\name{fitZTLogit}
\alias{fitZTLogit}
\title{Fit Capped Logistic Regression To Zero-Truncated Binomial Data}
\description{Estimate a logistic regression, with optionally capped probabilities, by maximum likelihood with zero-truncated data.}

\usage{
fitZTLogit(n.successes, n.trials, X = NULL, capped = FALSE, beta.start = NULL, alpha.start = 0.95)
}

\arguments{
  \item{n.successes}{number of binomial successes (numeric vector).
    Should be bounded below by 1 and bounded above by \code{n.trials}.}
  \item{n.trials}{number of binomial trials (numeric vector).}
  \item{X}{the regression design matrix. Number of rows should match \code{length(n.successes)}.}
  \item{capped}{if \code{TRUE}, then probability of a success will be capped at alpha < 1, where alpha is to be estimated.}
  \item{beta.start}{starting values for the regression coefficients. Of same length as \code{ncol(X)}.}
  \item{alpha.start}{starting value for alpha.}
}

\details{
Estimates a logistic regression equation for zero-truncated binomial observations.
Optionally estimates a limiting value for the probabilities that may be less than one.

The function maximizes the zero-truncated binomial likelihood using the \code{optim} function with \code{method="BFGS"}.
The fitted probabilities are equal to \code{alpha * plogis(X \%*\% beta)}.
}

\value{
A list with components
  \item{beta}{linear predictor coefficients.}
  \item{alpha}{capping parameter, maximum or asymptotic value for the probabilities.}
  \item{p}{fitted probabilities.}
  \item{deviance}{minus twice the maximized log-likelihood.}
  \item{calls}{number of function calls used in the optimization.}
}

\references{
Li M, Smyth GK (2023).
Neither random nor censored: estimating intensity-dependent probabilities for missing values in label-free proteomics.
\emph{Bioinformatics} 39(5), btad200.
\url{10.1093/bioinformatics/btad200}
}

\examples{
# Generate binomial data
n <- 30
n.trials <- rep(4,n)
x <- seq(from=3, to=9, length.out=n)
X <- model.matrix(~x)
beta <- c(-4,0.7)
p <- plogis(X \%*\% beta)
n.successes <- rbinom(n, size=n.trials, prob=p)

# Zero truncation
is.pos <- (n.successes > 0)
n.successes <- n.successes[is.pos]
n.trials <- n.trials[is.pos]
x <- x[is.pos]
X <- X[is.pos,]

# Zero-truncated regression
fit <- fitZTLogit(n.successes, n.trials, X)
p.observed <- n.successes / n.trials
plot(x, p.observed)
lines(x, fit$p)
}
