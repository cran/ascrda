\name{rdaSVM}
\alias{rdaSVM}
\title{Support Vector Machine after gene selection using SCRDA}
\description{
Fits Support Vector Machine after gene selection using SCRDA. It computes mis-classification rates on training and test data.
}
\usage{
rdaSVM(X, y, Xt, yt, alpha = seq(0, 0.99, 0.11), delta = seq(0, 3, 0.2), ...)
}

\arguments{
  \item{X}{training, expression matrix }
  \item{y}{response vector taking values 1, 2, ...}
  \item{Xt}{test, expression matrix }
  \item{yt}{ response vector for test data taking values 1, 2, ... }
  \item{alpha}{ tuning parameter }
  \item{delta}{ turning parameter }
  \item{...}{ optional arguments, passed to rda }
}

\value{
Error rates for training and test data
}


\references{ 
Guo Y, Hastie T, Tibshirani R (2007). 
Regularized linear discriminant analysis and its application in microarrays. 
Biostatistics 8, 86-100. 
}

\author{M. S. Islam and A. I. McLeod, email: mislam22@uwo.ca}

\examples{
## Example
#G = 400 genes
#fE=0.05 mean 5% of genes are expressed at 0.5 rest at 0.0
#Takes about 2 seconds
M<-synma(n=c(50,50), nt=c(100,100), B=c(20,4), m=c(20,100), fE=0.05)
rdaSVM(M$X, M$y, M$Xt, M$yt) 
}

\keyword{ classif }
