\name{sda}
\alias{sda}
\alias{predict.sda}

\title{Shrinkage Discriminant Analysis}

\description{
  \code{sda} trains a LDA or DDA classifier using Stein-type shrinkage estimation.
 
  \code{predict.sda} performs the corresponding class prediction.
}

\usage{
sda(Xtrain, L, diagonal=FALSE, fdr=FALSE, plot.fdr=FALSE, verbose=TRUE)
\method{predict}{sda}(object, Xtest, feature.idx, verbose=TRUE, ...)
}

\arguments{
  \item{Xtrain}{A matrix  containing the training data set. Note that 
                the rows are sample observations and the columns
                are variables.}
  \item{L}{A factor with the class labels of the training samples. }
  \item{diagonal}{Chooses between LDA (default, \code{diagonal=FALSE}) and DDA (\code{diagonal=TRUE}).}
  \item{fdr}{compute FDR values for each feature.}
  \item{plot.fdr}{show plot with estimated FDR values.}
  \item{verbose}{Report shrinkage intensities (sda) and number of used features (predict.sda).}
  \item{object}{An \code{sda} fit object obtained from the function \code{sda}.}
  \item{Xtest}{A matrix containing the test data set.}
  \item{feature.idx}{A vector indicating which features to employ for prediction (if unspecified all
                     features will be used).}
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
}

\value{
  \code{sda} trains the classifier and returns an \code{sda} object
  with the following components needed for the subsequent prediction:

  \item{regularization}{a vector containing the three estimated shrinkage intensities,}
   \item{prior}{the estimated class frequencies,}
   \item{predcoef}{matrix containing the coefficients used for prediction, and}
   \item{ranking}{matrix containing the ``correlation-adjusted t scores'' for each feature
                 and group.  The overall ranking of a feature is determined by the 
                 sum of the squared cat scores across all groups.}
   \code{predict.sda} predicts class probabilities for each test sample and returns
   a list with two components:
  \item{class}{a factor with the most probable class assignment, and}
   \item{posterior}{a matrix containing the class posterior probabilities for each test sample.}
}


\author{
  Miiika Ahdesm\"aki and Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\seealso{\code{\link[entropy:entropy.shrink]{freqs.shrink}}, 
\code{\link[corpcor:cov.shrink]{var.shrink}}, 
\code{\link[corpcor:invcov.shrink]{invcor.shrink}}.}

\examples{
# load sda library
library("sda")

# load full Khan et al (2001) data set
data(khan2001)
dim(khan2001$x)
levels(khan2001$y)

# create data set containing only the SRBCT samples
del.idx = which( khan2001$y == "non-SRBCT" )
srbct.x = khan2001$x[-del.idx,]
srbct.y = factor(khan2001$y[-del.idx])
dim(srbct.x)
levels(srbct.y)

# divide into training and test data
train.x = srbct.x[1:63,]
train.y = srbct.y[1:63]
test.x = srbct.x[64:83,]
test.y = srbct.y[64:83]

###################################################
# classification with correlation (shrinkage LDA) #
###################################################

sda.fit = sda(train.x, train.y)
ynew = predict(sda.fit, test.x)$class # using all 2308 features
sum(ynew != test.y) # 0

sda.fit$ranking[1:20,]
fidx = sda.fit$ranking[1:20,"idx"]
ynew = predict(sda.fit, test.x, feature.idx = fidx)$class # using the top 20 features
sum(ynew != test.y) # 1

###########################################################
# classification with diagonal covariance (shrinkage DDA) #
###########################################################

sda.fit = sda(train.x, train.y, diagonal=TRUE)
ynew = predict(sda.fit, test.x)$class # using all 2308 features
sum(ynew != test.y) # 4

sda.fit$ranking[1:20,]
fidx = sda.fit$ranking[1:20,"idx"]
ynew = predict(sda.fit, test.x, feature.idx = fidx)$class # using the top 20 features
sum(ynew != test.y) # 2
}
\keyword{multivariate}
