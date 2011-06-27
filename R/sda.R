### sda.R  (2011-06-26)
###
###    Shrinkage discriminant analysis (training the classifier)
###
### Copyright 2008-11 Miika Ahdesmaki and Korbinian Strimmer
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


sda = function(Xtrain, L, diagonal=FALSE, verbose=TRUE)
{
  if (!is.matrix(Xtrain)) stop("Training data must be given as matrix!")
  if (missing(L)) stop("Class labels are missing!")

  # shrinkage intensities
  regularization = rep(NA, 3)
  names(regularization) = c("lambda.freqs", "lambda.var", "lambda")
  regularization[3] = 1 # for diagonal=TRUE


  tmp = centroids(Xtrain, L, var.groups=FALSE, centered.data=TRUE,
    shrink=TRUE, verbose=verbose)
  
  cl.count = ncol(tmp$means)-1    # number of classes
 
  mu = tmp$means[,1:cl.count]
  mup = tmp$means[,cl.count+1]
  s2 = tmp$variances[,1]
  sc = sqrt(s2)  
  regularization[2] = attr(tmp$variances, "lambda.var")[1]

  nk = tmp$samples       # samples in class k
  n = sum(nk)            # number of samples
  p = nrow(mu)           # number of features
 
  xc = tmp$centered.data # to compute inverse pooled correlation matrix

  rm(tmp)

  
  ############################################################# 
  # compute coefficients for prediction 
  #############################################################

  # class frequencies
  prior = freqs.shrink( nk, verbose=verbose )
  regularization[1] = attr(prior, "lambda.freqs")
  attr(prior, "lambda.freqs") = NULL

  # reference means
  ref = array(0, dim=c(p, cl.count))
  colnames(ref) = paste("ref.", colnames(mu), sep="")
  rownames(ref) = rownames(mu)
  for (k in 1:cl.count)
  {
    ref[,k] = (mu[,k]+mup)/2
  }

  # prediction weights
  pw = array(0, dim=c(p, cl.count) )
  colnames(pw) = paste("pw.", colnames(mu), sep="")
  rownames(pw) = rownames(mu)

  for (k in 1:cl.count)
  {
    diff = mu[,k]-mup  
    pw[,k] = diff/sc
  }


  if(!diagonal)
  {
    if(verbose) cat("\nComputing inverse correlation matrix (pooled across classes)\n")
    pw = crossprod.powcor.shrink(xc, pw, alpha=-1, verbose=FALSE)
    regularization[3] = attr(pw, "lambda")
    attr(pw, "lambda") = NULL
    if(verbose) cat("Estimating optimal shrinkage intensity lambda (correlation matrix):", 
                    round(regularization[3], 4), "\n")
    
  }

  for (k in 1:cl.count)
  {
    pw[,k] = pw[,k]/sc
  }


  ############################################################# 

  out = list(regularization=regularization, prior=prior, 
             predcoef=cbind(ref, pw))
  class(out)="sda"

  return (out)
}

