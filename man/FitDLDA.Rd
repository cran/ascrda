\name{FitDLDA}
\alias{FitDLDA}
\title{ Diagonal linear/quadratic discriminant analysis }

\description{ Diagnonal linear/quadratic DA}

\usage{
FitDLDA(X, Y, Xt, Yt, pool)
}

\arguments{
  \item{X}{ design matrix, training sample }
  \item{Y}{ output, training, consecutive integers}
  \item{Xt}{ design matrix, test sample }
  \item{Yt}{ output, test, consecutive integers}
  \item{pool}{ set to FALSE for quadratic}
}

\details{
FitDLDA is very fast. 
Often there is little difference between the linear and quadratic methods.
}

\value{
list with two components:
misclassification rate
and confusion matrix for test.
}

\author{ AIM }

\references{
R package: sfsmisc.
Utilities from Seminar fuer Statistik ETH Zurich.
Martin Maechler. 2010
}



\examples{
#Example 1. small example - takes about 1.23 sec
\dontrun{
set.seed(132561)#for exact reproducibility
M<-synma(n=c(50,50), nt=c(100,100), B=c(20,4), m=c(20,100), fE=0.05)
StartTime <- proc.time()[3]
ANS1 <- FitDLDA(M$X, M$y, M$Xt, M$yt)
EndTime1 <- proc.time()[3]
ANS2 <- FitDLDA(M$X, M$y, M$Xt, M$yt, pool=FALSE)
EndTime2 <- proc.time()[3]
TotalTime1 <-EndTime1-StartTime
TotalTime2 <-EndTime2-EndTime1
TotalTime <- c(TotalTime1, TotalTime2)
names(TotalTime) <- c("DLDA","DQDA")
TotalTime
ANS1
ANS2
}
#Example 2. large example
#training sample size = 200
#test sample size = 1000
#number of genes = 1000
\dontrun{set.seed(132561)#for exact reproducibility
BScenario <- list(c(50, 10),c(10, 50), c(10, 10),  c(50, 50))
mScenario <- list(c(20,100),c(100, 20),c(100, 100),c(20, 20))
M <- synma(n=c(100,100), nt=c(500,500), rho = c(0.9, 0.9),
  B=BScenario[[1]], m=mScenario[[1]], fE=0.05)
StartTime <- proc.time()[3]
ANS <- FitDLDA(M$X, M$y, M$Xt, M$yt)
EndTime <- proc.time()[3]
TotalTime <-EndTime-StartTime
ans <- c(ANS$Err, TotalTime)
names(ans) <- c("MisclassificationRate", "TotalTime")
ans
##> ans
##MisclassificationRate             TotalTime
##                0.312                 0.110
}

#Example 3. Spira Dataset
#Spira Dataset
y <- as.numeric(Spira$Patients$STATUS)
X <- Spira$X
genes <- Spira$genes
X <- X[,genes]
set.seed(132561)#for exact reproducibility
M <-  getTrainTest(X, y, d=0.25)
StartTime <- proc.time()[3]
ANS1 <- FitDLDA(M$X, M$y, M$Xt, M$yt)
EndTime1 <- proc.time()[3]
ANS2 <- FitDLDA(M$X, M$y, M$Xt, M$yt, pool=FALSE)
EndTime2 <- proc.time()[3]
TotalTime1 <-EndTime1-StartTime
TotalTime2 <-EndTime2-EndTime1
TotalTime <- c(TotalTime1, TotalTime2)
names(TotalTime) <- c("DLDA","DQDA")
TotalTime
ANS1
ANS2

}

\keyword{ classif }

