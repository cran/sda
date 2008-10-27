### sda.R  (2008-10-26)
###
###    Shrinkage discriminant analysis (training)
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
   p = ncol(Xtrain)
   n = nrow(Xtrain)
   
   y = factor(L)
   if (length(y) != n) stop("for each sample there must be a class label!")
   cl = levels(y)
   clcount = length(cl)

   # shrinkage intensities
   regularization = rep(NA, 3)
   names(regularization) = c("lambda.freqs", "lambda.var", "lambda")
   

   # class frequencies
   f = freqs.shrink( table(y), verbose=verbose )
   prior = numeric(clcount)
   names(prior) = cl
   for (k in 1:clcount) prior[k] = f[cl[k]]
   regularization[1] = attr(f, "lambda.freqs")
   
   # means
   mu = array(NA, dim=c(p, clcount))
   colnames(mu) = cl
   rownames(mu) = colnames(Xtrain)
   Xc = array(0, dim=c(n,p))  # centered data
   colnames(Xc) = colnames(Xtrain)
   for (k in 1:clcount)
   {
      idx = ( y == cl[k] )
      Xk = Xtrain[idx,]
      mu[,k] = apply(Xk, 2, mean)
      Xc[idx,] = sweep(Xk, 2, mu[,k]) # center data
   }
 
   if (diagonal)
   {
      invvar = 1/var.shrink(Xc, verbose=verbose) *(n-clcount)/(n-1)
      regularization[2] = attr(invvar, "lambda.var")
      regularization[3] = 1
      invS = invvar
   }
   else
   {
     invS = invcov.shrink(Xc, verbose=verbose)*(n-clcount)/(n-1)
     regularization[2] = attr(invS, "lambda.var")
     regularization[3] = attr(invS, "lambda")
     attr(invS,"lambda") = NULL
     attr(invS,"lambda.estimated") = NULL
   }
   attr(invS,"lambda.var") = NULL
   attr(invS,"lambda.var.estimated") = NULL
   class(invS) = NULL
  
   out = list(regularization=regularization, prior=prior, means=mu, invcov=invS)
   class(out)="sda"

   return (out)
}


