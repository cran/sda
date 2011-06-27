\name{catscore}
\alias{catscore}


\title{Estimate CAT scores and t-scores}

\description{
  \code{catscore} computes  CAT scores
   (correlation-adjusted t-scores)
   between the group centroids and the pooled mean. 
}

\usage{
catscore(Xtrain, L, diagonal=FALSE, shrink=FALSE, verbose=TRUE)
}

\arguments{
  \item{Xtrain}{A matrix  containing the training data set. Note that 
                the rows correspond to observations and the columns
                to variables.}
  \item{L}{A factor with the class labels of the training samples. }
  \item{diagonal}{for \code{diagonal=FALSE} (the default) CAT scores are computed;
                otherwise with \code{diagonal=TRUE} t-scores.}
  \item{shrink}{Use empirical estimates or a shrinkage estimator for the CAT score.}
  \item{verbose}{Print out some info while computing.}
}
\details{CAT scores generalize conventional t-scores 
to account for correlation among predictors (Zuber and Strimmer 2009).
If there is no correlation then CAR scores reduce to t-scores.  The squared
CAR scores provide a decomposition of Hotelling's T^2 statistic.

CAT scores for two classes are described in Zuber and Strimmer (2009),
for the multi-class case see Ahdesm\"aki and Strimmer (2010).
}

\value{

   \code{catscore} returns a matrix containing the cat score (or t-score) between 
   each group centroid and the pooled mean for each feature. 
}
\references{
Ahdesm\"aki, A., and K. Strimmer. 2010.  Feature selection in omics prediction problems 
using cat scores and false non-discovery rate control. Ann. Appl. Stat. 4: 503-519.
Preprint available from \url{http://arxiv.org/abs/0903.2003}.

Zuber, V., and K. Strimmer. 2009.  Gene ranking and biomarker discovery under correlation.
Bioinformatics 25: 2700-2707.
Preprint available from \url{http://arxiv.org/abs/0902.0751}.
}

\author{
 Verena Zuber, Miika Ahdesm\"aki and Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\seealso{\code{\link{sda.ranking}}, \code{\link[care]{carscore}},.}

\examples{
# load sda library
library("sda")

################# 
# training data #
#################

# prostate cancer set
data(singh2002)

# training data
Xtrain = singh2002$x
Ytrain = singh2002$y
dim(Xtrain)


####################################################
# shrinkage t-score (DDA setting - no correlation) #
####################################################

tstat = catscore(Xtrain, Ytrain, diagonal=TRUE, shrink=TRUE)
dim(tstat)
tstat[1:10,]


########################################################
# shrinkage CAT score (LDA setting - with correlation) #
########################################################

cat = catscore(Xtrain, Ytrain, diagonal=FALSE, shrink=TRUE)
dim(cat)
cat[1:10,]

}
\keyword{multivariate}
