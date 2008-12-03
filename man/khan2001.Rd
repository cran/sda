\name{khan2001}
\alias{khan2001}

\title{Childhood Cancer Study of Khan et al. (2001)}


\usage{
data(khan2001)
}

\description{
Gene expression data (2308 genes for 88 samples) from the
microarray study of Khan et al. (2001).  
}

\details{ This data set contains measurements of the gene
expression of 2308 genes for 88 observations:  29 cases of Ewing 
sarcoma (EWS), 11 cases of Burkitt lymphoma (BL), 18 cases of 
neuroblastoma (NB), 25 cases of rhabdomyosarcoma (RMS), 
and 5 other (non-SRBCT) samples.
}

\format{
  \code{khan.data$x} is a 88 x 2308 matrix containing the expression levels. Note that
   rows correspond to samples, and columns to genes. The row names are the original
   image IDs, and the column names the orginal sample labels.

   \code{khan2001$y} is a factor containing the diagnosis for each sample ("BL", "EWS", "NB", "non-SRBCT", "RMS").
 
   \code{khan2001$descr} provides some annotation for each gene.
}

\source{
  The data are described in Khan et al. (2001) and can be
  obtained from 
  \url{http://research.nhgri.nih.gov/microarray/Supplement/}.
}

\references{
 Khan et al. 2001.
 Classification and diagnostic prediction of cancers using gene expression 
 profiling and artificial neural networks. Nature Medicine 7:673--679.
}


\examples{
# load sda library
library("sda")

# load full Khan et al (2001) data set
data(khan2001)
dim(khan2001$x)
khan2001$y

# create data set containing only the SRBCT samples
idx = which( khan2001$y == "non-SRBCT" )
srbct.x = khan2001$x[-idx,]
srbct.y = factor(khan2001$y[-idx])
srbct.descr = khan2001$descr[-idx]
dim(srbct.x)
srbct.y
}


\keyword{datasets}
