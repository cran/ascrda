\name{FitPAM}
\alias{FitPAM}
\title{Shrunken Centroid (Predictive Analysis) DA }
\description{
Class prediction using nearest centroid using regularization.
}

\usage{
FitPAM(X, Y, Xt, Yt, PamMethod=c("default", "kNN1", "kNN"))
}

\arguments{
  \item{X}{ design matrix, training sample }
  \item{Y}{ output, training, consecutive integers}
  \item{Xt}{ design matrix, test sample }
  \item{Yt}{ output, test, consecutive integers}
  \item{PamMethod}{ one of "default", "kNN1", "kNN". Default is no autocovariate adjustment}
}


\details{
See references.
}

\value{
list with two components:
Err = misclassification rate
confusionMatrix = confusion matrix for test.
}

\references{
R. Tibshirani et al. (2003),
"Class Prediction by Nearest Shrunken Centroids, with Applications to DNA Microarrays"
Statistical Science, 18, 104-117.
}

\author{
A. I. McLeod, aimcleod@uwo.ca
}

\seealso{
\link[pamr]{pamr.train} 
}

\examples{
#takes about 0.5 sec
#Spira Dataset
y <- as.numeric(Spira$Patients$STATUS)
X <- Spira$X
genes <- Spira$genes
X <- X[,genes]
set.seed(132561)#for exact reproducibility
M <-  getTrainTest(X, y, d=0.25)
StartTime <- proc.time()[3]
ANS <- FitPAM(M$X, M$y, M$Xt, M$yt)
EndTime <- proc.time()[3]
TotalTime <-EndTime-StartTime
TotalTime
ANS

#Further examples (continuing from above)
FitPAM(M$X, M$y, M$Xt, M$yt, PamMethod="kNN1")
FitPAM(M$X, M$y, M$Xt, M$yt, PamMethod="kNN")

}

\keyword{ classification }
