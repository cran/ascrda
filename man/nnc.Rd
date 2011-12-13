\name{nnc}
\alias{nnc}
\title{Computes nearest neighbor autocovariates}

\description{
Given n training examples in inputs X and output Y with two classes,
the corresponding nearest neighbor autocovariates are computed.
When there are more than 2 classes, X and Y contain the inputs and outputs for
two of the classes and Xb contains inputs for all other classes.
}

\usage{nnc(X, Y, k)}

\arguments{
  \item{X}{inputs, matrix with columns corresponding to inputs and rows to examples}
  \item{Y}{output, vector of length n, where n is the number of examples. 
Vector may be numeric, factor or character.}
  \item{k}{neighborhood size}
}

\details{
The k-th nearest neighbor autovariate is defined between two classes which we
may denote by -1 and 1.
For a given observation, we find all the k-nearest neighbors, and sum
the number correspond to to each class.
Then take the difference and divide by k.
Symbolically,  
\deqn{\sum_{j in k_i} (I(y_j = 1) - I(y_j = -1))/k}
{sum {j in k_i} (I(y_j = 1) - I(y_j = -1))/k}
where eqn{k_i}{k_i} denotes the indicies of all observations
in the k-th neighbor of the i-th observation.

}

\value{vector of length n of nearest neighbor autocovariates corresponding to training data 
for given k}

\references{
Holmes, C. C. and Adams, N. M. (2003). 
Likelihood inference in nearest-neighbor classification models, 
Biometrika, 90(1), 99-112.
\url{http://biomet.oxfordjournals.org/cgi/content/abstract/90/1/99}
}

\author{
M. S. Islam and A. I. McLeod
Maintainer: <aimcleod@uwo.ca>
}

\seealso{
\code{\link[class]{knn}}, 
\code{\link[class]{knn.cv}} 
}

\examples{
#Example 1. With 2 classes.
library("MASS")
X <- synth.tr[,1:2]
Y <- synth.tr[,3]
# Nearest neighbor auto-covariates for synthetic data
sapply(1:10, FUN = function(k) nnc(X=X, Y=Y, k)) 

#Example 2. Iris dataset with 3 classes.
#Find both autocovariates
X <- iris[,1:4]
Y<- iris[,5]
k<-15
z<-nnc(X=X, Y=Y, k=k)
}

\keyword{classif}
