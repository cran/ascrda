\name{FitSVM}
\alias{FitSVM}
\title{Classification using Support Vector Machine }
\description{
Fits Support Vector Machine, computes mis-classification rates on training and test data.
}
\usage{
FitSVM(X, y, Xt, yt)
}

\arguments{
  \item{X}{training, expression matrix }
  \item{y}{response vector taking values 1, 2, ...}
  \item{Xt}{test, expression matrix }
  \item{yt}{ response vector for test data taking values 1, 2, ... }
  }

\value{
Error rates for training and test data
}


\references{ 
Hastie T, Tibshirani R, Frienman J: The Elements of Statistical Learning: Data Mining, Inference, and
Prediction. Springer-Verlag 2009.
}

\author{M. S. Islam and A. I. McLeod, email: mislam22@uwo.ca}

\examples{
#G = 400 genes
#fE=0.05 mean 5% of genes are expressed at 0.5 rest at 0.0
#Takes about 2 seconds
M<-synma(n=c(50,50), nt=c(100,100), B=c(20,4), m=c(20,100), fE=0.05)
FitSVM(M$X, M$y, M$Xt, M$yt)
}

\keyword{ classif }
