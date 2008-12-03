\encoding{latin1}
\name{sda-package}
\alias{sda-package}
\docType{package}
\title{The sda package}

\description{This package performs linear and diagonal discriminant analysis
where the classifier is trained using Stein-type shrinkage estimators.

Specifically, three separate shrinkage estimators are employed for the

\item{class frequencies}{following Hausser and Strimmer (2008),} 
\item{variances}{following Opgen-Rhein and Strimmer (2007), and}
\item{correlations}{following Sch\"afer and Strimmer (2005).} 

The estimates are plugged into the LDA and DDA discriminant
function for prediction.  Note that the three corresponding regularization parameters
are obtained analytically without resorting to computer intensive resampling.

This approach is particularly suited for high-dimensional classification.
}

\author{Korbinian Strimmer (\url{http://strimmerlab.org/})}
\references{
  See website: \url{http://strimmerlab.org/software/sda/}
}
\keyword{multivariate}
\seealso{
\code{\link{sda}, \link{predict.sda}}
}
