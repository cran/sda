\name{sda}
\alias{sda}
\alias{predict.sda}

\title{Shrinkage discriminant analysis}

\usage{
sda(Xtrain, L, diagonal=FALSE, verbose=TRUE)
\method{predict}{sda}(object, Xtest, ...)
}
\description{
  \code{rda} trains a LDA or DDA classifier using Stein-type shrinkage estimation.
 
  \code{predict.sda} performs the corresponding class prediction.
}

\arguments{
  \item{Xtrain}{A matrix  containing the training data set. Note that 
                the rows are sample observations and the columns
                are variables.}
  \item{L}{A factor with the class labels of the traning samples. }
  \item{diagonal}{Chooses between LDA (default, \code{diagonal=FALSE}) and DDA (\code{diagonal=TRUE}).}
  \item{verbose}{Report shrinkage intensities.}
  \item{object}{An \code{sda} fit object obtained from the function \code{sda}.}
  \item{Xtest}{A matrix containing the test data set..}
  \item{...}{Additional arguments for generic predict.}
}
\details{
In order to train the LDA or DDA classifier, three separate shrinkage estimators are employed:

\item{class frequencies}{ the  estimator \code{\link[entropy:entropy.shrink]{freqs.shrink}} 
           from   Hausser and Strimmer (2008),}
\item{variances}{the estimator \code{\link[corpcor:cov.shrink]{var.shrink}} from Opgen-Rhein and Strimmer (2007),  }

\item{correlations}{the estimator \code{\link[corpcor:invcov.shrink]{invcor.shrink}} from Sch\"afer and Strimmer (2005).  }

These estimates are plugged into the LDA and DDA discriminant
function for prediction.  Note that the three corresponding regularization parameters
are obtained analytically without resorting to computer intensive resampling.

This approach is particularly suited for high-dimensional classification.
}
\value{
  \code{sda} trains the classifier and returns an \code{sda} object
  with the following components:

  \item{regularization}{a vector containing the three estimated shrinkage intensities,}
   \item{prior}{the estimated class frequencies,}
   \item{means}{a matrix containing the group means (centroids), and}
   \item{invcov}{the inverse of the estimated pooled covariance matrix (this is a matrix in the LDA and
    a vector in the DDA case).}

   \code{predict.sda} predicts class probabilities for each test sample and returns
   a list with two components:
  \item{yhat}{a factor with the the most likely class assignment, and}
   \item{probs}{a matrix containing the class probabilities for each test sample.}

}


\author{
  Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\seealso{\code{\link[entropy:entropy.shrink]{freqs.shrink}}, 
\code{\link[corpcor:cov.shrink]{var.shrink}}, 
\code{\link[corpcor:invcov.shrink]{invcor.shrink}}.}

\examples{
library("sda")

## prepare data set
data(iris) # good old iris data
X = as.matrix(iris[,1:4])
Y = iris[,5]

# divide into test and traing data set
tr.index = sample(1:length(Y), 2/3*length(Y))
train.x = X[tr.index,]                     
train.y = Y[tr.index]
test.x = X[-tr.index,] 
test.y = Y[-tr.index]


## shrinkage LDA
sda.fit = sda(Xtrain=train.x, L=train.y)
sda.fit

predict(sda.fit, test.x)

ynew = predict(sda.fit, test.x)$yhat
sum(ynew != test.y)


## shrinkage DDA
sda.fit = sda(Xtrain=train.x, L=train.y, diagonal=TRUE)
sda.fit

predict(sda.fit, test.x)

ynew = predict(sda.fit, test.x)$yhat
sum(ynew != test.y)
}
\keyword{univar}
