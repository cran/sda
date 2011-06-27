### centroids.R  (2011-06-26)
###
###    Group centroids, variances, and correlations
###
### Copyright 2008-2011 Korbinian Strimmer
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



centroids = function(x, L,  
  var.groups=FALSE, centered.data=FALSE, shrink=FALSE, verbose=TRUE)
{
  if (!is.matrix(x)) stop("Input x must be a matrix!")
  p = ncol(x)
  n = nrow(x)

  if (length(L) != n) stop("For each sample there must be a class label!")   
  
  groups = pvt.groups(L)
  samples = groups$samples
  cl.count = groups$cl.count
  cl.names = groups$cl.names

  if (verbose)
  {
    cat("Number of variables:", p, "\n")
    cat("Number of observations:", n, "\n")
    cat("Number of classes:", cl.count, "\n\n")
  }

  
  # setup arrays
  mu = array(0, dim=c(p, cl.count+1))
  colnames(mu) = c(cl.names, "(pooled)")
  rownames(mu) = colnames(x)

  xc = array(0, dim=c(n,p))  # storage for centered data
  rownames(xc) = rownames(x)
  colnames(xc) = colnames(x)
 
  if(var.groups)
  {
    v = array(0, dim=c(p, cl.count+1)) # storage for variances
    if(shrink) attr(v, "lambda.var") = numeric(cl.count+1)
    colnames(v) = c(cl.names, "(pooled)")
    rownames(v) = colnames(x)
  }
  else
  {
    v = array(0, dim=c(p, 1)) # store only pooled variances
    if(shrink) attr(v, "lambda.var") = numeric(1)
    colnames(v) = c("(pooled)")
    rownames(v) = colnames(x)
  }  

  # compute means and variance in each group
  for (k in 1:cl.count)
  {
     idx = groups$idx[,k]
     Xk = x[ idx, ,drop = FALSE]
     mu[,k] = colMeans(Xk)

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
  
  # compute pooled mean and variance
  mu[,cl.count+1] = colMeans(x) # pooled mean
 
  if (verbose) cat("Estimating variances (pooled across classes)\n")
  if (var.groups)
  {
    if (shrink)
      v.pool = var.shrink(xc, verbose=verbose)
    else
      v.pool = as.vector(var.shrink(xc, lambda.var=0, verbose=FALSE))
    v[,cl.count+1] = v.pool*(n-1)/(n-cl.count) # correction factor
    if (shrink) attr(v, "lambda.var")[cl.count+1] = attr(v.pool, "lambda.var")
  }  
  else
  {
    if (shrink)
      v.pool = var.shrink(xc, verbose=verbose)
    else
      v.pool = as.vector(var.shrink(xc, lambda.var=0, verbose=FALSE))
    v[,1] = v.pool*(n-1)/(n-cl.count) # correction factor
    if (shrink) attr(v, "lambda.var")[1] = attr(v.pool, "lambda.var")
  }

  if(centered.data == FALSE) xc=NULL

  ##

  return( list(samples=samples, means=mu, variances=v, 
     centered.data=xc))
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

