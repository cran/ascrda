\name{khat}
\alias{khat}

\title{Optimal k in kNN using pseudolikelihood}
\description{
Uses the profile pseudolikelihood to obtain the estimate for k, the number of
nearest neighbors parameter in kNN.
}

\usage{khat(X, Y, kmax = ceiling(length(Y) * 0.5), plot = TRUE)}

\arguments{
  \item{X}{An n-by-p matrix of covariates}
  \item{Y}{Outputs with Q classes}
  \item{kmax}{The maximum size of k}
  \item{plot}{if TRUE, plot the profile deviance otherwise no plot}
}

\details{
When Q=2, the glm algorithm is used to compute the profile
pseudologlikelihood and for Q>2, the function \code{multinom}
in \pkg{nnet} is used.
}

\value{
The estimate of k obtained by maximizing the pseudolikelihood
is returned.

The result is returned invisibly if plot is TRUE. 
}

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
\code{\link[nnet]{multinom}} 
}


\examples{

#Two classes example
library("MASS")#need synth.tr
X <- synth.tr[,1:2]
Y <- synth.tr[,3]
khat(X=X, Y=Y)

#Three classes example
library("MASS") #need lda
Y<- iris[,5]
X<- iris[,1:4]
kopt<-khat(X, Y)
kopt
#Mis-classification rates on training data.
#Of course FLDA does better in this case.
y <- factor(Y)
ans <-knn(train=X, test=X, k=kopt, cl=y)
etaKNN <- sum(ans!=y)/length(y)
iris.ldf <-lda(X, y)
yfitFLDA <- predict(iris.ldf, newdata=X, dimen=1)$class
etaFLDA <- sum(yfitFLDA!=y)/length(y)
eta<-c(etaFLDA, etaKNN)
names(eta)<-c("FLDA", "kNN")
eta
}

\keyword{ classif }

