### predict.sda.R  (2008-11-15)
###
###    Shrinkage discriminant analysis (prediction)
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




predict.sda = function(object, Xtest, ...)
{
  if ( missing(object) ) {
    stop("A sda fit object must be supplied.")
  }

  if ( missing(Xtest) ) {
    stop("A new data to predict must be supplied.")
  }

  m = object$means
  freq = object$prior
  iS = object$invcov
  nt = nrow(Xtest)
  clcount = length(freq)

  if( is.null(dim(iS)) )
    diagonal = TRUE
  else
    diagonal = FALSE

  probs = array(0, dim=c(nt, clcount) )
  score = numeric(clcount) 
  yhat = integer(nt)
  for (i in 1:nt)
  {
    x = t(Xtest[i,,drop=FALSE])  # p x 1
    for (k in 1:clcount)
    {
       if (diagonal)
       {
         tmp = m[,k]*iS
         score[k] = sum( tmp*x -0.5*tmp*m[,k] ) + log(freq[k])
       }
       else
       {
         tmp = crossprod(iS, m[,k])
         score[k] = crossprod(tmp, x) -0.5*crossprod(tmp, m[,k]) + log(freq[k])
       }
    }
    probs[i,] = score2prob(score)
    yhat[i] = which.max(score)
  }
  probs = zapsmall(probs)
  attr(yhat, "levels") = names(freq)
  class(yhat)= "factor"
  colnames(probs) = names(freq)
  rownames(probs) = rownames(Xtest)

  return(list(yhat=yhat, probs=probs) )
}

score2prob = function(x)
{
   x = x-max(x)
   x = exp(x)
   x = x/sum(x)

   return(x)
}

