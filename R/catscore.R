### catscore.R  (2012-11-24)
###
###    Estimate CAT scores and t-scores
###
### Copyright 2008-12 Verena Zuber, Miika Ahdesmaki and Korbinian Strimmer
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


catscore = function(Xtrain, L, lambda, lambda.var, diagonal=FALSE, verbose=TRUE)
{
  if (!is.matrix(Xtrain)) stop("Training data must be given as matrix!")
  if (missing(L)) stop("Class labels are missing!")

  if(verbose) 
  {
    if(diagonal)
      cat("Computing t-scores (centroid vs. pooled mean) for feature ranking\n\n")
    else
      cat("Computing cat scores (centroid vs. pooled mean) for feature ranking\n\n")
  }

  tmp = centroids(Xtrain, L, lambda.var=lambda.var, var.groups=FALSE, centered.data=TRUE, verbose=verbose)

  cl.count = ncol(tmp$means)-1    # number of classes
  
  mu = tmp$means[,1:cl.count]
  mup = tmp$means[,cl.count+1]
  s2 = tmp$variances[,1]
  sc = sqrt(s2)
  lambda.var = attr(tmp$variances,"lambda.var")[1]
  lambda.var.estimated = attr(tmp$variances,"lambda.var.estimated")
 
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
  attr(cat, "lambda.var") = lambda.var
  attr(cat, "lambda.var.estimated") = lambda.var.estimated
  class(cat) = "shrinkage"

  # then compute cat scores
  if (!diagonal)
  {
    if(verbose) cat("Computing the square root of the inverse pooled correlation matrix\n")     
    cat = crossprod.powcor.shrink(xc, cat, alpha=-1/2, lambda=lambda, verbose=FALSE)
    attr(cat, "lambda.var") = lambda.var
    attr(cat, "lambda.var.estimated") = lambda.var.estimated 

    if(verbose)
    {
       if(attr(cat, "lambda.estimated") )
         cat("Estimating optimal shrinkage intensity lambda (correlation matrix):", 
         round(attr(cat, "lambda"), 4), "\n")
       else
         cat("Specified shrinkage intensity lambda (correlation matrix):", 
         round(attr(cat, "lambda"), 4), "\n")
    }

  }

  return(cat)
}

