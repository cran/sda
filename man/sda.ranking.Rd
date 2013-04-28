\name{sda.ranking}
\alias{sda.ranking}
\alias{plot.sda.ranking}

\title{Shrinkage Discriminant Analysis 1: Predictor Ranking}

\description{
  \code{sda.ranking} determines a ranking of predictors by computing  CAT scores
   (correlation-adjusted t-scores)
   between the group centroids and the pooled mean.

  \code{plot.sda.ranking} provides a graphical visualization of the top ranking features..
}

\usage{
sda.ranking(Xtrain, L, lambda, lambda.var, diagonal=FALSE, fdr=TRUE, plot.fdr=FALSE, verbose=TRUE)
\method{plot}{sda.ranking}(x, top=40, arrow.col="blue", ...)
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
  \item{diagonal}{Chooses between LDA (default, \code{diagonal=FALSE}) and DDA (\code{diagonal=TRUE}).}
  \item{fdr}{compute FDR values and HC scores for each feature.}
  \item{plot.fdr}{Show plot with estimated FDR values.}
  \item{verbose}{Print out some info while computing.}
  \item{x}{An "sda.ranking" object -- this is produced by the sda.ranking() function.}
  \item{top}{The number of top-ranking features shown in the plot (default: 40).}
  \item{arrow.col}{Color of the arrows in the plot (default is \code{"blue"}).}
   \item{...}{Additional arguments for generic plot.}
}
\details{
For each predictor variable and centroid a shrinkage CAT scores of the mean versus
the pooled mean is computed.  The overall ranking of a feature is determined
by the sum of the squared cat scores across all centroids.
For the diagonal case (LDA) the (shrinkage) CAT score reduces to the (shrinkage) t-score.  Thus in the two-class diagonal case the feature are simply ranked according to the
(shrinkage) t-scores.

Calling \code{sda.ranking} is step 1 in a classification analysis with the
sda package.  Steps 2 and 3 are 
\code{\link{sda}} and \code{\link{predict.sda}}

See Ahdesm\"aki and Strimmer (2010) for details on multi-class CAT scores,
Zuber and Strimmer (2009) for CAT scores in general.  For shrinkage t scores
see Opgen-Rhein and Strimmer (2007). 
}

\value{

   \code{sda.ranking} returns a matrix with the following columns:

  \item{idx}{original feature number}
   \item{score}{sum of the squared CAT scores across groups - this determines the overall ranking of a feature}
   \item{cat}{for each group and feature the cat score of the centroid versus the pooled mean}

  If \code{fdr=TRUE} then additionally local false discovery rate (FDR) values 
  as well as higher criticism (HC) scores are computed for each feature
  (using \code{\link[fdrtool]{fdrtool}}).

}
\references{
Ahdesm\"aki, A., and K. Strimmer. 2010.  Feature selection in omics prediction problems 
using cat scores and false non-discovery rate control. Ann. Appl. Stat. 4: 503-519.
Preprint available from \url{http://arxiv.org/abs/0903.2003}.

Opgen-Rhein, R., and K. Strimmer. 2007. Accurate ranking of differentially expressed 
genes by a distribution-free shrinkage approach.  
 Statist. Appl. Genet. Mol. Biol. 6:9.

Zuber, V., and K. Strimmer. 2009.  Gene ranking and biomarker discovery under correlation.
Bioinformatics 25: 2700-2707.
Preprint available from \url{http://arxiv.org/abs/0902.0751}.}

\author{
  Miika Ahdesm\"aki, Verena Zuber, Sebastian Gibb, and Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\seealso{\code{\link{catscore}}, \code{\link{sda}}, \code{\link{predict.sda}}.}

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

######################################### 
# feature ranking (diagonal covariance) #
#########################################

# ranking using t-scores (DDA)
ranking.DDA = sda.ranking(Xtrain, Ytrain, diagonal=TRUE)
ranking.DDA[1:10,]

# plot t-scores for the top 40 genes
plot(ranking.DDA, top=40) 

# number of features with local FDR < 0.8 
# (i.e. features useful for prediction)
sum(ranking.DDA[,"lfdr"] < 0.8)

# number of features with local FDR < 0.2 
# (i.e. significant non-null features)
sum(ranking.DDA[,"lfdr"] < 0.2)

# optimal feature set according to HC score
plot(ranking.DDA[,"HC"], type="l")
which.max( ranking.DDA[1:1000,"HC"] ) 


##################################### 
# feature ranking (full covariance) #
#####################################

# ranking using CAT-scores (LDA)
ranking.LDA = sda.ranking(Xtrain, Ytrain, diagonal=FALSE)
ranking.LDA[1:10,]

# plot t-scores for the top 40 genes
plot(ranking.LDA, top=40) 

# number of features with local FDR < 0.8 
# (i.e. features useful for prediction)
sum(ranking.LDA[,"lfdr"] < 0.8)

# number of features with local FDR < 0.2 
# (i.e. significant non-null features)
sum(ranking.LDA[,"lfdr"] < 0.2)

# optimal feature set according to HC score
plot(ranking.LDA[,"HC"], type="l")
which.max( ranking.LDA[1:1000,"HC"] ) 

}
\keyword{multivariate}
