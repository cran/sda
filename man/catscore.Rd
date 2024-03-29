\name{catscore}
\alias{catscore}


\title{Estimate CAT Scores and t-Scores}

\description{
  \code{catscore} computes  CAT scores
   (correlation-adjusted t-scores)
   between the group centroids and the pooled mean. 
}

\usage{
catscore(Xtrain, L, lambda, lambda.var, lambda.freqs, diagonal=FALSE, verbose=TRUE)
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
  \item{diagonal}{for \code{diagonal=FALSE} (the default) CAT scores are computed;
                otherwise with \code{diagonal=TRUE} t-scores.}
  \item{verbose}{Print out some info while computing.}
}
\details{CAT scores generalize conventional t-scores 
to account for correlation among predictors (Zuber and Strimmer 2009).
If there is no correlation then CAR scores reduce to t-scores.  The squared
CAR scores provide a decomposition of Hotelling's T^2 statistic.

CAT scores for two classes are described in Zuber and Strimmer (2009),
for the multi-class case see Ahdesm\"aki and Strimmer (2010).

The scale factors for t-scores and CAT-scores are computed from the estimated frequencies
(for empirical scale factors set  \code{lambda.freqs=0}).
}

\value{

   \code{catscore} returns a matrix containing the cat score (or t-score) between 
   each group centroid and the pooled mean for each feature. 
}
\references{
Ahdesm\"aki, A., and K. Strimmer. 2010.  Feature selection in omics prediction problems 
using cat scores and false non-discovery rate control. Ann. Appl. Stat. 4: 503-519.
<DOI:10.1214/09-AOAS277>

Zuber, V., and K. Strimmer. 2009.  Gene ranking and biomarker discovery under correlation.
Bioinformatics 25: 2700-2707.
<DOI:10.1093/bioinformatics/btp460>
}

\author{
 Verena Zuber, Miika Ahdesm\"aki and Korbinian Strimmer (\url{https://strimmerlab.github.io}).
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

tstat = catscore(Xtrain, Ytrain, diagonal=TRUE)
dim(tstat)
tstat[1:10,]


########################################################
# shrinkage CAT score (LDA setting - with correlation) #
########################################################

cat = catscore(Xtrain, Ytrain, diagonal=FALSE)
dim(cat)
cat[1:10,]

}
\keyword{multivariate}
