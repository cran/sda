### centroids.R  (2008-11-20)
###
###    Group centroids, variances and inverse correlation matrix
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



centroids = function(x, L, var.pooled=TRUE, var.groups=FALSE, invcor.pooled=FALSE, 
  shrink=FALSE, verbose=TRUE)
{
  p = ncol(x)
  n = nrow(x)

  if (length(L) != n) stop("For each sample there must be a class label!")   
  
  groups = pvt.groups(L)
  samples = groups$samples
  cl.count = groups$cl.count
  cl.names = groups$cl.names

  if (verbose) cat("Number of class labels: ", cl.count, "\n")

  # means
  mu = array(NA, dim=c(p, cl.count))
  colnames(mu) = cl.names
  rownames(mu) = colnames(x)
 
  if(var.pooled)
  {
     xc = array(0, dim=c(n,p))  # storage for centered data
     #colnames(xc) = colnames(x)
  }

  if(var.groups)
  {
    # storage for variances
    v = array(0, dim=c(p, cl.count))
    colnames(v) = c(cl.names)
    rownames(v) = colnames(x)
    if (shrink) attr(v, "lambda.var") = numeric(cl.count)
  }
  else
  {
    v = NULL
  }

  for (k in 1:cl.count)
  {
     idx = groups$idx[,k]
     Xk = x[ idx, ,drop = FALSE]
     mu[,k] = colMeans(Xk)

     if(var.pooled)
       xc[idx,] = sweep(Xk, 2, mu[,k]) # center data

     if (var.groups)
     {
       if(verbose) cat("Estimating variances (class #", k, ")\n", sep="")
        if (shrink)
        {
          vs = var.shrink(Xk, verbose=verbose)
          v[,k] = as.vector(vs)
          attr(v, "lambda.var")[k] = attr(vs, "lambda.var")
        }
        else
          v[,k] = as.vector(var.shrink(Xk, lambda.var=0, verbose=FALSE))
     }
  }

  if (var.pooled)
  {
    if (verbose) cat("Estimating variances (pooled across classes)\n")
   
    if (shrink)
    {
      v.pool = var.shrink(xc, verbose=verbose)
    }
    else
    {
      v.pool = as.vector(var.shrink(xc, lambda.var=0, verbose=FALSE))
      attr(v.pool, "lambda.var") = NULL
    }
    attr(v.pool, "class") = NULL
    attr(v.pool, "lambda.var.estimated") = NULL
    
    v.pool = v.pool*(n-1)/(n-cl.count)
    names(v.pool) = colnames(x)
  }
  else
  {
    v.pool = NULL
  }
  
  if (invcor.pooled == TRUE)
  {
    if (verbose) cat("Estimating inverse correlation matrix (pooled across classes)\n")

    if (shrink)
    {
      ic = invcor.shrink(xc, collapse=TRUE, verbose=verbose)
    }
    else
    {
      ic = invcor.shrink(xc, lambda=0, collapse=TRUE, verbose=FALSE)
      attr(ic, "lambda") = NULL
    }
    attr(ic, "class") = NULL
    attr(ic, "lambda.estimated") = NULL      

    # note there is no correction factor for (inverse) correlation
    colnames(ic) = colnames(x)
    rownames(ic) = colnames(x)
  }
  else
  {
    ic = NULL
  }


  return( list(samples=samples, means=mu, var.pooled=v.pool, var.groups=v, invcor.pooled=ic))
}



## private function ##

pvt.groups = function(L)
{
   y = factor(L)  # note that this creates new levels (in contrast to as.factor)
     
   cl.names = levels(y)
   cl.count = length(cl.names)

   idx = array(FALSE, dim = c(length(y), cl.count))
   colnames(idx) = cl.names
   nn = integer(cl.count)
   names(nn) = cl.names
   
   for (k in 1:cl.count)
   {
      idx[,k] = ( y == cl.names[k] )
      nn[k] = sum(idx[,k])
   }

   return( list(idx=idx, samples=nn, cl.count=cl.count, cl.names=cl.names) )
}

