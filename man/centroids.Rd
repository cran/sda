\name{centroids}
\alias{centroids}

\title{Group Centroids and (Pooled) Variances}

\usage{
centroids(x, L, lambda.var, lambda.freqs, var.groups=FALSE, 
  centered.data=FALSE, verbose=TRUE)
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
  \item{lambda.var}{Shrinkage intensity for the variances. If not specified it is 
    estimated from the data, see details below. \code{lambda.var=0} implies no shrinkage
    and \code{lambda.var=1} complete shrinkage. }
  \item{lambda.freqs}{Shrinkage intensity for the frequencies. If not specified it is 
    estimated from the data. \code{lambda.freqs=0} implies no shrinkage (i.e. empirical frequencies)
    and \code{lambda.freqs=1} complete shrinkage (i.e. uniform frequencies). }
  \item{var.groups}{Estimate group-specific variances.}
  \item{centered.data}{Return column-centered data matrix.}
  \item{verbose}{Provide some messages while computing.}
}
\details{
As estimator of the variance we employ 
 \code{\link[corpcor:cov.shrink]{var.shrink}} as described in Opgen-Rhein and Strimmer (2007).
For the estimates of frequencies we rely on 
\code{\link[entropy:freqs.shrink]{freqs.shrink}} as described in Hausser and Strimmer (2009).
Note that the pooled mean is computed using the estimated frequencies.
}
\value{
  \code{centroids} returns a list
  with the following components:
  \item{samples}{a vector containing the samples sizes in each group,}
  \item{freqs}{a vector containing the estimated frequency in each group,}
   \item{means}{the group means and the pooled mean,}
   \item{variances}{the group-specific and the pooled variances, and}
   \item{centered.data}{a matrix containing the centered data.}
}


\author{
  Korbinian Strimmer (\url{https://strimmerlab.github.io}).
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
centroids(X, Y, lambda.var=0)
          
## also compute group-specific variances
centroids(X, Y, var.groups=TRUE, lambda.var=0)
   
## use shrinkage estimator for the variances
centroids(X, Y, var.groups=TRUE)

## return centered data
xc = centroids(X, Y, centered.data=TRUE)$centered.data
apply(xc, 2, mean)

## useful, e.g., to compute the inverse pooled correlation matrix
powcor.shrink(xc, alpha=-1)
}
\keyword{multivariate}
