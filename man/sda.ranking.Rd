\name{sda.ranking}
\alias{sda.ranking}
\alias{plot.sda.ranking}

\title{Shrinkage Discriminant Analysis 1: Predictor Ranking}

\description{
  \code{sda.ranking} determines a ranking of predictors by computing  CAT scores
   (correlation-adjusted t-scores)
   between the group centroids and the pooled mean.

  \code{plot.sda.ranking} provides a graphical visualization of the top ranking features.
}

\usage{
sda.ranking(Xtrain, L, lambda, lambda.var, lambda.freqs,
  ranking.score=c("entropy", "avg", "max"), 
  diagonal=FALSE, fdr=TRUE, plot.fdr=FALSE, verbose=TRUE)
\method{plot}{sda.ranking}(x, top=40, arrow.col="blue", zeroaxis.col="red",
   ylab="Features", main, ...)
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
  \item{ranking.score}{how to compute the summary score for each variable from the CAT scores of all classes - see Details.}
  \item{fdr}{compute FDR values and HC scores for each feature.}
  \item{plot.fdr}{Show plot with estimated FDR values.}
  \item{verbose}{Print out some info while computing.}
  \item{x}{An "sda.ranking" object -- this is produced by the sda.ranking() function.}
  \item{top}{The number of top-ranking features shown in the plot (default: 40).}
  \item{arrow.col}{Color of the arrows in the plot (default is \code{"blue"}).}
  \item{zeroaxis.col}{Color for the center zero axis (default is \code{"red"}).}
  \item{ylab}{Label written next to feature list (default is \code{"Features"}).}
  \item{main}{Main title (if missing, \code{"The", top, "Top Ranking Features"} is used).}
  \item{...}{Other options passed on to generic plot().}
}
\details{
For each predictor variable and centroid a shrinkage CAT scores of the mean versus
the pooled mean is computed.   If there are only two classes the CAT score vs. the pooled
mean reduces to the CAT score between the two group means.
Moreover, in the diagonal case (LDA) the (shrinkage) CAT score reduces to the (shrinkage) t-score.

The overall ranking of a feature is determine by computing a summary score from the CAT scores.
This is controlled by the option \code{ranking.score}.  The default setting
(\code{ranking.score="entropy"}) uses mutual information
between the response and the respective predictors (\code{ranking.score}) for ranking.  This is equivalent to
a weighted sum of squared CAT scores across the classes.  Another possibility is to employ
the average of the squared CAT scores for ranking (as suggested in Ahdesm\"aki and Strimmer 2010) 
by setting \code{ranking.score="avg"}.  A third option is to use the maximum of the squared CAT scores across groups (similarly as in the PAM algorithm) via setting \code{ranking.score="max"}.
Note that in the case of two classes all three options are equivalent and
 lead to identical scores.  Thus, the choice of \code{ranking.score} is important only
in the multi-class setting.  In the two-class case the features are simply ranked according to the
(shrinkage) squared CAT-scores (or t-scores, if there is no correlation among predictors).

The current default approach is to use ranking by mutual information (i.e. relative entropy
between full model vs. model without predictor) and to use shrinkage estimators of frequencies.
In order to reproduce exactly the ranking computed by previous versions (1.1.0 to 1.3.0) of the \code{sda} package set the options  \code{ranking.score="avg"} and \code{lambda.freqs=0}.

Calling \code{sda.ranking} is step 1 in a classification analysis with the
sda package.  Steps 2 and 3 are 
\code{\link{sda}} and \code{\link{predict.sda}}

See Zuber and Strimmer (2009) for CAT scores in general, and 
Ahdesm\"aki and Strimmer (2010) for details on multi-class CAT scores.
For shrinkage t scores see Opgen-Rhein and Strimmer (2007). 
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
<DOI:10.1214/09-AOAS277>

Opgen-Rhein, R., and K. Strimmer. 2007. Accurate ranking of differentially expressed 
genes by a distribution-free shrinkage approach.  
 Statist. Appl. Genet. Mol. Biol. 6:9.
<DOI:10.2202/1544-6115.1252>

Zuber, V., and K. Strimmer. 2009.  Gene ranking and biomarker discovery under correlation.
Bioinformatics 25: 2700-2707. <DOI:10.1093/bioinformatics/btp460>
}

\author{
  Miika Ahdesm\"aki, Verena Zuber, Sebastian Gibb, and Korbinian Strimmer (\url{https://strimmerlab.github.io}).
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
