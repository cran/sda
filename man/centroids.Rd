\name{centroids}
\alias{centroids}

\title{Group Centroids, (Pooled) Variances, and Inverse Pooled Correlation Matrix}

\usage{
centroids(x, L, var.pooled=TRUE, var.groups=FALSE, invcor.pooled=FALSE, shrink=FALSE, verbose=TRUE)
}
\description{
  \code{centroids} computes group centroids and optionally the pooled variance, the group specific variances, and
                    the inverse of the pooled correlation matrix.
}

\arguments{
  \item{x}{A matrix  containing the  data set. Note that 
                the rows are sample observations and the columns
                are variables.}
  \item{L}{A factor with the group labels. }
  \item{var.pooled}{Estimate the pooled variances.}
  \item{var.groups}{Estimate all group-specific variances.}
  \item{invcor.pooled}{Estimate inverse pooled correlation matrix.}
  \item{shrink}{Use empirical or shrinkage estimator.}
  \item{verbose}{Provide some messages while computing.}
}
\details{
If option \code{shrink=TRUE} then the shrinkage estimators 
 \code{\link[corpcor:cov.shrink]{var.shrink}} from Opgen-Rhein and Strimmer (2007)
and  \code{\link[corpcor:invcov.shrink]{invcor.shrink}} from Sch\"afer and Strimmer (2005) are used.
}
\value{
  \code{centroids} returns a list
  with the following components:

  \item{samples}{a vector containing the samples sizes in each group,}
   \item{means}{the empirical group means,}
   \item{var.pooled}{a vector containing the pooled variances,}
   \item{var.groups}{a matrix containing the group-specific variances, and}
   \item{invcor.pooled}{a matrix containing the inverse pooled correlation matrix 
                        (if all correlations are zero a vector of 1s is returned to save space). }
}


\author{
  Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\seealso{
\code{\link[corpcor:cov.shrink]{var.shrink}}, 
\code{\link[corpcor:invcov.shrink]{invcor.shrink}}.}

\examples{
library("sda")

## prepare data set
data(iris) # good old iris data
X = as.matrix(iris[,1:4])
Y = iris[,5]

## estimate centroids and empirical pooled variances
centroids(X, Y)

## include also group-specific variances
centroids(X, Y, var.groups=TRUE)

## and inverse pooled correlation
centroids(X, Y, var.groups=TRUE, invcor.pooled=TRUE)

## use shrinkage estimator for variances and correlations
centroids(X, Y, var.groups=TRUE, invcor.pooled=TRUE, shrink=TRUE)
}
\keyword{multivariate}
