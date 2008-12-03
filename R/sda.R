### sda.R  (2008-12-03)
###
###    Shrinkage discriminant analysis (training the classifier)
###
### Copyright 2008 Korbinian Strimmer
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
   if (missing(L)) stop("Class labels are missing!")

   # shrinkage intensities
   regularization = rep(NA, 3)
   names(regularization) = c("lambda.freqs", "lambda.var", "lambda")

   # centroids, variances and inverse correlation matrix
   if (diagonal) # no correlations in this case
   {
     tmp = centroids(Xtrain, L, mean.pooled=TRUE, var.pooled=TRUE, var.groups=FALSE, 
            powcor.pooled=FALSE, shrink=TRUE, verbose=verbose)
     regularization[2] = attr(tmp$var.pooled, "lambda.var")
     regularization[3] = 1
   }
   else
   {
     tmp = centroids(Xtrain, L, mean.pooled=TRUE, var.pooled=TRUE, var.groups=FALSE, 
            powcor.pooled=TRUE, alpha=-1/2, shrink=TRUE, verbose=verbose)
     regularization[2] = attr(tmp$var.pooled, "lambda.var")
     regularization[3] = attr(tmp$powcor.pooled, "lambda")

     if(verbose) cat("Computing inverse correlation matrix (pooled across classes)\n\n")
     ic = centroids(Xtrain, L, mean.pooled=FALSE, var.pooled=FALSE, var.groups=FALSE, 
            powcor.pooled=TRUE, alpha=-1, shrink=TRUE, verbose=FALSE)$powcor.pooled
   }

   n = sum(tmp$samples)        # number of samples
   p = nrow(tmp$means)         # number of features
   cl.count = ncol(tmp$means)  # number of classes
 
   # class frequencies
   prior = freqs.shrink( tmp$samples, verbose=verbose )
   regularization[1] = attr(prior, "lambda.freqs")
   attr(prior, "lambda.freqs") = NULL
   
   # correlation adjusted two-sample t-statistic centroid versus pooled mean
   cat = array(0, dim=c(p, cl.count))
   coef = array(0, dim=c(p, cl.count))
   colnames(cat) = colnames(tmp$means)
   rownames(cat) = rownames(tmp$means)
   colnames(coef) = colnames(tmp$means)
   rownames(coef) = rownames(tmp$means)
   s = sqrt(tmp$var.pooled)
   attr(s, "lambda.var") = NULL
   m = sqrt(1/tmp$samples - 1/n) # note the minus sign!
   for (i in 1:cl.count)
   {
      diff = tmp$means[,i]-tmp$mean.pooled
      cat[,i] = diff/(m[i]*s)           # t score for feature selection
      coef[,i] = diff/tmp$var.pooled    # coefficient for prediction
   }
   if (diagonal==FALSE) # adjust for correlations between features
   {
     if (!is.null(dim(tmp$powcor.pooled)))
     {
       cat = crossprod(tmp$powcor.pooled, cat) # decorrelate
       coef = crossprod(ic, coef)              # multiply by inverse correlation 
     }
     rm(ic)
   } 
   
   # reference means
   ref = array(0, dim=c(p, cl.count))
   colnames(ref) = colnames(tmp$means)
   rownames(ref) = rownames(tmp$means)
   for (i in 1:cl.count)
   {
      ref[,i] = (tmp$means[,i]+tmp$mean.pooled)/2
   }
   rm(tmp)

   out = list(regularization=regularization, prior=prior, ref=ref, coef=coef, cat=cat)
   class(out)="sda"

   return (out)
}


