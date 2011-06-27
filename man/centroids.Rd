\name{centroids}
\alias{centroids}

\title{Group Centroids and (Pooled) Variances}

\usage{
centroids(x, L, var.groups=FALSE, centered.data=FALSE, shrink=FALSE, verbose=TRUE)
}
\description{
  \code{centroids} computes group centroids, the pooled mean 
   and pooled variance, and optionally the group specific variances.
}

\arguments{
  \item{x}{A matrix  containing the  data set. Note that 
                the rows are sample observations and the columns
                are variables.}
  \item{L}{A factor with the group labels. }
  \item{var.groups}{Estimate group-specific variances.}
  \item{centered.data}{Return column-centered data matrix.}
  \item{shrink}{Use empirical estimates or a shrinkage estimator for the variances.}
  \item{verbose}{Provide some messages while computing.}
}
\details{
If option \code{shrink=TRUE} then the shrinkage estimators 
 \code{\link[corpcor:cov.shrink]{var.shrink}} from Opgen-Rhein and Strimmer (2007)
and  \code{\link[corpcor:cov.shrink]{cor.shrink}} from Sch\"afer and Strimmer (2005) are used.

Details on the algorithm for efficiently computing the power of the shrinkage correlation 
matrix are given in Zuber and Strimmer (2009).
}
\value{
  \code{centroids} returns a list
  with the following components:
  \item{samples}{a vector containing the samples sizes in each group,}
   \item{means}{the group means and the pooled mean,}
   \item{variances}{the group-specific and the pooled variances, and}
   \item{centered.data}{a matrix containing the centered data.}
}


\author{
  Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\seealso{
\code{\link[corpcor:cov.shrink]{var.shrink}}, 
\code{\link[corpcor]{powcor.shrink}}.}

\examples{
# load sda library
library("sda")

## prepare data set
data(iris) # good old iris data
X = as.matrix(iris[,1:4])
Y = iris[,5]

## estimate centroids and empirical pooled variances
centroids(X, Y)
          
## also compute group-specific variances
centroids(X, Y, var.groups=TRUE)
   
## use shrinkage estimator for the variances
centroids(X, Y, var.groups=TRUE, shrink=TRUE)

## return centered data
xc = centroids(X, Y, centered.data=TRUE)$centered.data
apply(xc, 2, mean)

## compute pooled inverse correlation matrix
powcor.shrink(xc, alpha=-1)
}
\keyword{multivariate}
