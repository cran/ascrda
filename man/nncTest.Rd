\name{nncTest}
\alias{nncTest}
\title{Nearest neighbor autocovariates for test data}

\description{The relative proportion of class j to class 1 based on the
training data is computed for each point in the test data.
This is useful in cross-validation with regularized discriminant
analysis and as well as with statistical models such as logistic
regression. 
}

\usage{nncTest(X, Y, XTest, k)}

\arguments{
  \item{X}{training inputs, n-by-p matrix, where n=sample size, p=number of input variables}
  \item{Y}{training outputs with K levels, length n}
  \item{XTest}{test inputs, m-by-p matrix, m=sample size}
  \item{k}{neighborhood size}
}

\value{vector of length m of nearest neighbor autocovariates corresponding to each
row of XTest}

\references{
Holmes, C. C. and Adams, N. M. (2003). 
Likelihood inference in nearest-neighbour classification models, 
Biometrika, 90(1), 99-112.
\url{http://biomet.oxfordjournals.org/cgi/content/abstract/90/1/99}
}

\author{
M. S. Islam and A. I. McLeod
Maintainer: <aimcleod@uwo.ca>
}

\seealso{
\code{\link{nnc}},
\code{\link[class]{knn}}, 
\code{\link[class]{knn.cv}} 
}

\examples{
#Example 1. With 2 classes.
library("MASS") #need synth.tr
X <- synth.tr[,1:2]
Y <- synth.tr[,3]
XTest <- synth.te[,1:2]
k <- 21
nncTest(X=X, Y=Y, XTest=XTest, k) 
}

\keyword{classif}


