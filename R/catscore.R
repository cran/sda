### catscore.R  (2011-06-26)
###
###    Estimate CAT scores and t-scores
###
### Copyright 2008-11 Verena Zuber, Miika Ahdesmaki and Korbinian Strimmer
###
###
### This file is part of the `sda' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 3, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA


catscore = function(Xtrain, L, diagonal=FALSE, shrink=FALSE, verbose=TRUE)
{
  if (!is.matrix(Xtrain)) stop("Training data must be given as matrix!")
  if (missing(L)) stop("Class labels are missing!")

  if(verbose) 
  {
    if (shrink) 
      type="shrinkage"
    else 
      type="empirical"

    if(diagonal)
      cat("Computing", type, "t-scores (centroid vs. pooled mean) for feature ranking\n\n")
    else
      cat("Computing", type, "cat scores (centroid vs. pooled mean) for feature ranking\n\n")
  }

  tmp = centroids(Xtrain, L, var.groups=FALSE, centered.data=TRUE, shrink=shrink, verbose=verbose)

  cl.count = ncol(tmp$means)-1    # number of classes
  
  mu = tmp$means[,1:cl.count]
  mup = tmp$means[,cl.count+1]
  s2 = tmp$variances[,1]
  sc = sqrt(s2)
  if(shrink) 
  {
    lambda.var = attr(tmp$variances,"lambda.var")[1]
    lambda.var.estimated = TRUE
  }

  nk = tmp$samples       # samples in class k
  n = sum(nk)            # number of samples
  p = nrow(mu)           # number of features
  
  xc = tmp$centered.data # to compute pooled correlation matrix

  rm(tmp)


  ############################################################# 
  # compute coefficients for feature ranking
  #############################################################

  cat = array(0, dim=c(p, cl.count))
  if(diagonal)
    colnames(cat) = paste("t.", colnames(mu), sep="")
  else
    colnames(cat) = paste("cat.", colnames(mu), sep="")
  rownames(cat) = rownames(mu)
  
  # first compute t-scores (centroid vs. pooled mean)
  m = sqrt(1/nk - 1/n) # note the minus sign!
  for (k in 1:cl.count)
  {
    diff = mu[,k]-mup
    cat[,k] = diff/(m[k]*sc)   # t-scores
  }
  if (shrink) 
  {
    attr(cat, "lambda.var") = lambda.var
    attr(cat, "lambda.var.estimated") = lambda.var.estimated
    class(cat) = "shrinkage"
  }

  # then compute cat scores
  if (!diagonal)
  {

    if(verbose) cat("Computing the square root of the inverse pooled", type, "correlation matrix\n")
    if (shrink)
    {
      cat = crossprod.powcor.shrink(xc, cat, alpha=-1/2, verbose=FALSE)
      if(verbose) cat("Estimating optimal shrinkage intensity lambda (correlation matrix):", 
         round(attr(cat, "lambda"), 4), "\n")
      attr(cat, "lambda.var") = lambda.var
      attr(cat, "lambda.var.estimated") = lambda.var.estimated
    }
    else
    {
      cat = crossprod.powcor.shrink(xc, cat, alpha=-1/2, lambda=0, verbose=FALSE)
    }
  }

  return(cat)
}

