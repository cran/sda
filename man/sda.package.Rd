\encoding{latin1}
\name{sda-package}
\alias{sda-package}
\docType{package}
\title{The sda package}

\description{This package performs linear discriminant analysis (LDA)
and diagonal discriminant analysis (DDA)
with variable selection using correlation-adjusted t (CAT) scores.

The classifier is trained using James-Stein-type shrinkage estimators.
Variable selection is based on ranking predictors by CAT scores (LDA) 
or t-scores (DDA).  A cutoff is chosen by false non-discovery 
rate (FNDR) or  higher criticism (HC) thresholding.

This approach is particularly suited for high-dimensional classification
with correlation among predictors.
For details see Zuber and Strimmer (2009) and Ahdesm\"aki and Strimmer (2010).

Typically the functions in this package are applied in three steps:
\itemize{
\item{Step 1:}{feature selection with \code{\link{sda.ranking}},} 
\item{Step 2:}{training the classifier with \code{\link{sda}}, and}
\item{Step 3:}{classification using \code{\link{predict.sda}}.} 
}
The accompanying web site (see below) provides example R scripts to illustrate
the functionality of this package.
}

\references{
Ahdesm\"aki, A., and K. Strimmer. 2010.  Feature selection in omics prediction problems 
using cat scores and false non-discovery rate control. Ann. Appl. Stat. 4: 503-519.
Preprint available from \url{http://arxiv.org/abs/0903.2003}.

Zuber, V., and K. Strimmer. 2009.  Gene ranking and biomarker discovery under correlation.
Bioinformatics 25: 2700-2707.
Preprint available from \url{http://arxiv.org/abs/0902.0751}.

See website: \url{http://strimmerlab.org/software/sda/}
}

\author{Miika Ahdesm\"aki, Verena Zuber and Korbinian Strimmer (\url{http://strimmerlab.org/})}

\keyword{multivariate}
\seealso{
 \code{\link{catscore}}, \code{\link{sda.ranking}}, \code{\link{sda}}, \code{\link{predict.sda}}. 
}
