### sda.ranking.R  (2011-06-26)
###
###    Shrinkage discriminant analysis (feature ranking)
###
### Copyright 2008-11 Miika Ahdesmaki, Verena Zuber and Korbinian Strimmer
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


sda.ranking = function(Xtrain, L, diagonal=FALSE, fdr=TRUE, plot.fdr=FALSE, verbose=TRUE)
{
  cat = catscore(Xtrain, L, diagonal=diagonal, shrink=TRUE, verbose=verbose)

  cl.count = dim(cat)[2]

  score = apply(cat^2, 1, sum) # sum of squared CAT-scores
  names(score) = rownames(cat)
  idx = order(score, decreasing = TRUE)

  if (fdr)
  {
    if (verbose) cat("\nComputing false discovery rates and higher cricitism scores for each feature\n")

    if (cl.count == 2)
    {
      fdr.out = fdrtool(cat[,1], plot=plot.fdr, verbose=FALSE)
    }
    else
    {
      z = score^(1/3) # Wilson-Hilferty transformation to normality
     
      # center before feeding into fdrtool
      #offset = median(z)
      d = density(z)
      offset = d$x[which.max(d$y)]
      z = z-offset 
      fdr.out = fdrtool(z, plot=plot.fdr, verbose=FALSE)
    }
    lfdr = fdr.out$lfdr # local false discovery rates
    pval = fdr.out$pval # p-values

    # compute HC score for each p-value
    HC = pvt.hcobj(pval)

    ranking = cbind(idx, score[idx], cat[idx, , drop=FALSE], lfdr[idx], HC[idx])
    colnames(ranking) = c("idx", "score", colnames(cat), "lfdr", "HC")
    rm(fdr.out)
  }
  else
  {
    ranking = cbind(idx, score[idx], cat[idx, , drop=FALSE])
    colnames(ranking) = c("idx", "score", colnames(cat))
  }
  rm(cat)

  attr(ranking, "class") = "sda.ranking"
  attr(ranking, "diagonal") = diagonal
  attr(ranking, "cl.count") = cl.count

  return(ranking)
}


plot.sda.ranking = function(x, top=40, ...)
{
  if ( class(x) != "sda.ranking" )
    stop ("sda.ranking x needed as input!")

  cl.count = attr(x, "cl.count")
  diagonal = attr(x, "diagonal")
  if (diagonal) 
    xlab = "t-Scores (Centroid vs. Pooled Mean)"
  else 
    xlab = "Correlation-Adjusted t-Scores (Centroid vs. Pooled Mean)"

  top = min( nrow(x), top ) # just to be sure ...

  idx = 2+(1:cl.count)
  cn = colnames(x)[idx]
  if(diagonal)
    colnames(x)[idx] = substr(cn, 3, nchar(cn))
  else
    colnames(x)[idx] = substr(cn, 5, nchar(cn))

  if (is.null(rownames(x))) rownames(x) = x[, 1]

  score = x[1:top, 2]
  DATA = as.data.frame.table( x[1:top, idx] )

  require("lattice")
  dotplot(reorder(Var1,rep(score, cl.count )) ~ Freq | Var2, 
    data = DATA, origin = 0, type = c("p", "h"), 
    main = paste("The", top, "Top Ranking Features"), 
    xlab = xlab, 
    layout=c(cl.count,1), ...) 
}


