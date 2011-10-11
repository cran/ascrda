\name{ascrda}
\alias{ascrda}

\title{Test error rate for ASCRDA and other methods}

\description{The misclassifcation rate on the test data is computed for the following methods:
ASCRDA -  ASCRDA - SCRDA with nearest neightbor autocovariates;
SCRDA - regular SCRDA (shrunken centroid regularized discriminant analysis;
SVM - support vector machine;
SCRDA/SCRDA - SCRDA using feature selection with SCRDA;
SVM/SCRDA - SVM using feature selection with SCRDA. }

\usage{
ascrda(X, y, Xt, yt, k = "default", alpha = seq(0, 0.99, 0.11), delta = seq(0, 3, 0.2), ...)
}

\arguments{
  \item{X}{training, expression matrix}
  \item{y}{training output}
  \item{Xt}{test expression matrix}
  \item{yt}{test output}
  \item{k}{"default" use MLE k. Otherwise set to specified value, k>0.}
  \item{alpha}{ alpha in rda function}
  \item{delta}{ delta in rda function}
  \item{\dots}{ optional arguments passed to rda}
}


\value{Vector of length 5, misclassification error rates for ASCRDA, SCRDA, SVM, SCRDA/SCRDA, SVM/SCRDA.}

\references{
Guo Y, Hastie T, Tibshirani R (2007). 
Regularized linear discriminant analysis and its application in microarrays. 
Biostatistics 8, 86-100. 

Islam, M. S. (2008). Periodicity, Change Detection and Prediction in Microarrays. Ph.D. thesis, The University of Western Ontario, Canada.
}

\author{M. S. Islam and A. I. McLeod, email: mislam22@uwo.ca}

\examples{
#G = 400 genes
#fE=0.05 mean 5% of genes are expressed at 0.5 rest at 0.0
#takes about 30 seconds
set.seed(132561) #may not work with other seeds, due to bug with rda
StartTime <- proc.time()[3]
M<-synma(n=c(50,50), nt=c(100,100), B=c(20,4), m=c(20,100), fE=0.05)
ANS<-ascrda(M$X, M$y, M$Xt, M$yt)
EndTime <- proc.time()[3]
TotalTime<-EndTime-StartTime
TotalTime
ANS
}
\keyword{classif}


