\name{sda}
\alias{sda}

\title{Shrinkage Discriminant Analysis 2: Training Step}

\description{
  \code{sda} trains a LDA or DDA classifier using James-Stein-type shrinkage estimation.
}

\usage{
sda(Xtrain, L, lambda, lambda.var, lambda.freqs, diagonal=FALSE, verbose=TRUE)
}

\arguments{
  \item{Xtrain}{A matrix  containing the training data set. Note that 
                the rows correspond to observations and the columns
                to variables.}
  \item{L}{A factor with the class labels of the training samples. }
  \item{lambda}{Shrinkage intensity for the correlation matrix. If not specified it is 
    estimated from the data. \code{lambda=0} implies no shrinkage
    and \code{lambda=1} complete shrinkage. }
  \item{lambda.var}{Shrinkage intensity for the variances. If not specified it is 
    estimated from the data. \code{lambda.var=0} implies no shrinkage
    and \code{lambda.var=1} complete shrinkage. }
  \item{lambda.freqs}{Shrinkage intensity for the frequencies. If not specified it is 
    estimated from the data. \code{lambda.freqs=0} implies no shrinkage (i.e. empirical frequencies)
    and \code{lambda.freqs=1} complete shrinkage (i.e. uniform frequencies). }
 
  \item{diagonal}{Chooses between LDA (default, \code{diagonal=FALSE}) and DDA (\code{diagonal=TRUE}).}
  \item{verbose}{Print out some info while computing.}
}
\details{
In order to train the LDA or DDA classifier, three separate shrinkage estimators are employed:

\describe{
 \item{class frequencies:}{ the  estimator \code{\link[entropy:entropy.shrink]{freqs.shrink}} 
           from   Hausser and Strimmer (2008),}
 \item{variances:}{ the estimator \code{\link[corpcor:cov.shrink]{var.shrink}} from Opgen-Rhein and Strimmer (2007),  }

 \item{correlations:}{ the estimator \code{\link[corpcor:cov.shrink]{cor.shrink}} from Sch\"afer and Strimmer (2005).  }
}
Note that the three corresponding regularization parameters
are obtained analytically without resorting to computer intensive resampling.
}

\value{

   \code{sda} trains the classifier and returns an \code{sda} object
  with the following components needed for the subsequent prediction:

  \item{regularization}{a vector containing the three estimated shrinkage intensities,}
   \item{freqs}{the estimated class frequencies,}
   \item{alpha}{vector containing the intercepts used for prediction,}
   \item{beta}{matrix containing the coefficients used for prediction.}
}

\references{
Ahdesm\"aki, A., and K. Strimmer. 2010.  Feature selection in omics prediction problems 
using cat scores and false non-discovery rate control. Ann. Appl. Stat. 4: 503-519.
<DOI:10.1214/09-AOAS277>
}

\author{
  Miika Ahdesm\"aki and Korbinian Strimmer (\url{https://strimmerlab.github.io}).
}
\seealso{\code{\link{predict.sda}}, \code{\link{sda.ranking}},
\code{\link[entropy:entropy.shrink]{freqs.shrink}}, 
\code{\link[corpcor:cov.shrink]{var.shrink}}, 
\code{\link[corpcor:invcov.shrink]{invcor.shrink}}.}

\examples{
# load sda library
library("sda")

########################## 
# training and test data #
##########################

# data set containing the SRBCT samples
get.srbct = function()
{
  data(khan2001)
  idx = which( khan2001$y == "non-SRBCT" )
  x = khan2001$x[-idx,]
  y = factor(khan2001$y[-idx])
  descr = khan2001$descr[-idx]

  list(x=x, y=y, descr=descr)
}
srbct = get.srbct()

# training data
Xtrain = srbct$x[1:63,]
Ytrain = srbct$y[1:63]
Xtest = srbct$x[64:83,]
Ytest = srbct$y[64:83]


###################################################
# classification with correlation (shrinkage LDA) #
###################################################

sda.fit = sda(Xtrain, Ytrain)
ynew = predict(sda.fit, Xtest)$class # using all 2308 features
sum(ynew != Ytest) 

###########################################################
# classification with diagonal covariance (shrinkage DDA) #
###########################################################

sda.fit = sda(Xtrain, Ytrain, diagonal=TRUE)
ynew = predict(sda.fit, Xtest)$class # using all 2308 features
sum(ynew != Ytest) 

#################################################################
# for complete example scripts illustrating classification with #
# feature selection visit https://strimmerlab.github.io/software/sda/  #
#################################################################
}
\keyword{multivariate}
