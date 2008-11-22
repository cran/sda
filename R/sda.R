### sda.R  (2008-11-20)
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
   tmp = centroids(Xtrain, L, var.pooled=TRUE, var.groups=FALSE, 
            invcor.pooled=TRUE, shrink=TRUE, verbose=verbose)

   # class frequencies
   prior = freqs.shrink( tmp$samples, verbose=verbose )
   regularization[1] = attr(prior, "lambda.freqs")
   attr(prior, "lambda.freqs") = NULL
   
   # means
   mu = tmp$means
   
   if (diagonal)
   {
      regularization[2] = attr(tmp$var.pooled, "lambda.var")
      regularization[3] = 1

      invS  = 1/tmp$var.pooled

      attr(invS, "lambda.var") = NULL
   }
   else
   {
     sc = sqrt( tmp$var.pooled )
     invS = tmp$invcor.pooled

     regularization[2] = attr(sc, "lambda.var")
     regularization[3] = attr(invS, "lambda")

     if (is.null(dim(invS)))
       invS = invS/sc/sc
     else
       invS = sweep(sweep(invS, 1, 1/sc, "*"), 2, 1/sc, "*")

     attr(invS, "lambda") = NULL
   }
   
   out = list(regularization=regularization, prior=prior, means=mu, invcov=invS)
   class(out)="sda"

   return (out)
}


