\name{PMLcfaefa}
\alias{PMLcfaefa}
\title{Fit a PML CFA-EFA model through a compiled Julia library}
\usage{
PMLcfaefa(X, p = ncol(X), q, K, c = max(X),
  ics = matrix(1, nrow(X), 1), tol = 1e-6,
  maxit = 1000L, maxNumRestart = 20L, replix = 10L,
  fixthrehsolds = 0, L_str = FALSE,
  latentmeans = 1, parametr = FALSE, verbose = FALSE, fast = FALSE)
}
\arguments{
\item{X}{An \code{n x p} matrix of integer category responses, coded from 1 to \code{c}.}
\item{p}{The number of observed items. Must match \code{ncol(X)}.}
\item{q}{The number of CFA factors.}
\item{K}{The number of EFA factors.}
\item{c}{The number of ordinal response categories.}
\item{ics}{A numeric matrix of covariates. Supply an intercept column explicitly.}
\item{tol}{Convergence tolerance.}
\item{maxit}{Maximum number of EM iterations.}
\item{maxNumRestart}{Maximum number of restarts.}
\item{replix}{Number of initial value attempts.}
\item{fixthrehsolds}{Original Julia threshold mode (original spelling preserved).}
\item{L_str}{FALSE for the default CFA loading structure, or a numeric \code{p x q} matrix.}
\item{latentmeans}{Original Julia latent mean mode (0, 1, or 2).}
\item{parametr}{Not supported in the initial native bridge; keep FALSE.}
\item{verbose}{Logical. If TRUE, print the composite log-likelihood and the actual stopping differences at every EM iteration. FALSE by default. The primary stopping difference is absolute for latentmeans = 1 and fixthrehsolds = 0, and relative in the other supported branches.}
\item{fast}{Logical. If TRUE, aggregate the fixed EM weights by the c squared ordinal response pairs for the M-step objective and analytic gradient. Supported only when latentmeans=1 and fixthrehsolds=0. The default FALSE uses the original algorithm. Floating-point results may differ in the last digits due to summation order.}
}
\value{An object of class \code{pmlcfaefa_fit} with elements \code{raw},
\code{status_code}, \code{restart_count}, \code{n}, \code{p}, \code{q},
\code{K} and \code{c}. The \code{raw} element preserves the original Julia
numeric result; its parameter segmentation depends on model settings.}
\description{Calls the Julia PMLcfaefa algorithm from R through a compiled library created with PackageCompiler.}
