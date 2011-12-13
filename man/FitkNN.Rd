\name{FitkNN}
\alias{FitkNN}
\title{ kNN classification }
\description{
kNN classification
}
\usage{
FitkNN(X, Y, Xt, Yt, k = 1)
}

\arguments{
  \item{X}{ design matrix, training sample }
  \item{Y}{ response, training sample }
  \item{Xt}{ design matrix, test sample }
  \item{Yt}{ response, test  }
  \item{k}{number of nearest neighbours}
}


\value{
list with two components:
misclassification rates
and confusion matricies for both the training and test samples.
}
\author{ AIM }

\examples{
#Example 1. Small example with 2 classes.
set.seed(132561)#for exact reproducibility
M<-synma(n=c(50,50), nt=c(100,100), B=c(20,4), m=c(20,100), fE=0.05)
StartTime <- proc.time()[3]
ANS <- FitkNN(M$X, M$y, M$Xt, M$yt)
EndTime <- proc.time()[3]
TotalTime <-EndTime-StartTime
out <- c(ANS$r[2], TotalTime)
names(out) <- c("MisclassificationRate", "TotalTime")
out
#
#Example 2.  large example
#training sample size = 200
#test sample size = 1000
#number of genes = 1000
\dontrun{ #FitkNN is fast but not synma
set.seed(132561)#for exact reproducibility
BScenario <- list(c(50, 10),c(10, 50), c(10, 10),  c(50, 50))
mScenario <- list(c(20,100),c(100, 20),c(100, 100),c(20, 20))
M <- synma(n=c(100,100), nt=c(500,500), rho = c(0.9, 0.9),
      B=BScenario[[1]], m=mScenario[[1]], fE=0.05)
StartTime <- proc.time()[3]
ANS <- FitkNN(M$X, M$y, M$Xt, M$yt)
EndTime <- proc.time()[3]
TotalTime <-EndTime-StartTime
ans <- c(ANS$r[2], TotalTime)
names(ans) <- c("MisclassificationRate", "TotalTime")
ans
}
##> ans
##MisclassificationRate             TotalTime
##                0.185                 2.660
#
#Example 3. With 3 classes.  Still fast - about 0.3 sec
#Spira Dataset
y <- as.numeric(Spira$Patients$STATUS)
X <- Spira$X
genes <- Spira$genes
X <- X[,genes]
set.seed(132561)#for exact reproducibility
M <-  getTrainTest(X, y, d=0.25)
StartTime <- proc.time()[3]
ANS <- FitkNN(M$X, M$y, M$Xt, M$yt)
EndTime <- proc.time()[3]
TotalTime <- EndTime-StartTime
out <- c(ANS$r[2], TotalTime)
names(out) <- c("MisclassificationRate", "TotalTime")
out
ANS$CMats[[2]]
}

\keyword{ classif }


