\name{ascrda-package}
\alias{ascrda-package}
\docType{package}

\title{Augmented shrunken centroid discriminant analysis}


\description{
This paper uses nearest neighbour autocovariates (Holmes and Adams, 2003)
together with regularized shrunken centroid discriminant analysis (Guo, Hastie and
Tibshirani, 2007) for discriminant analysis or class prediction.
This technique was described in the thesis by Islam (2008).
This method is especially useful in small-n-big-p problems as in microarrays
where n, the number of subjects is usually less than 100 and p, the number
of genes is over 10,000.
Other functions implemented in this package include nearest neighbor
autocovariates, the delete-d method for cross-validation and a method
for simulating non-linear microarray data.
An interesting microarray dataset with smokers (Spira et al., 2004) is provided.
For convenience, a consistent interfact is provided for other classification
methods including nearest neighbours, diagonal linear discriminant analysis and
support vector machines.
}

\details{
\tabular{ll}{
Package: \tab ascrda\cr
Type: \tab Package\cr
Version: \tab 1.15\cr
Date: \tab 2011-12-13\cr
License: \tab GPL (>= 2.10)\cr
LazyLoad: \tab yes\cr
LazyData: \tab yes\cr
}
}

\author{
M. S. Islam and A. I. McLeod
 <aimcleod@uwo.ca>
}

\references{
Guo Y, Hastie T, Tibshirani R (2007). 
Regularized linear discriminant analysis and its application in microarrays. 
Biostatistics 8, 86-100. 

Holmes, C. C. and Adams, N. M. (2003). 
Likelihood inference in nearest-neighbour classification models, 
Biometrika, 90(1), 99-112.
\url{http://biomet.oxfordjournals.org/cgi/content/abstract/90/1/99}

Islam, M. S. (2008). Periodicity, Change Detection and Prediction in Microarrays.
Ph.D. Thesis, The University of Western Ontario.

Avrum Spira, Jennifer Beane, Vishal Shah, Gang Liu, Frank Schembri, Xuemei Yang, John Palma, and Jerome S. Brody
Effects of cigarette smoke on the human airway epithelial cell transcriptome
PNAS, Jul 2004; 101: 10143 - 10148. 


}


\keyword{ package }

