#Source: Synthetic.R
#
#number of simulations of each scenario. 
#So there are 4*NSIM simulations in total.
NSIM <- 4
#In our paper we report results with NSIM<-1000
#NSIM <- 1000 #takes about 23.3 hours
#
#Mac computer has 8 nodes. PC has 12 nodes.
nWorkers<-8
#load required libraries
require("ascrda")
require("parallel")
require("xtable")
require("lattice")
#The size of the job is quite different depending on the scenario.
#By default, we randomly allocated the jobs.
#Use randomization?
RandomizeQ <- TRUE
#If you set RandomizeQ <- FALSE, then the randomization feature is not
# used. Equivalent but not exactly equal results are obtained.
#
#Set output file base name.
#This is system and user dependent.
#Windows/Mac version:
if (.Platform$OS.type=="windows") BASE <- "d:/r/2011/ascrda/OUTTable1/" else 
	BASE <- "/Users/aim/R/2011/ascrda/OUTTable1/"
#
#file to save matrix of mis-classification rates in text form
Outtb  <- paste(BASE, "tb.txt", sep="")
OuttbMed  <- paste(BASE, "tbMed.txt", sep="")
#file for latex output
Outtbl <- paste(BASE, "tb.tex", sep="")
OuttblMed <- paste(BASE, "tbMed.tex", sep="")
#raw output
OUTRaw <- paste(BASE, "OUT.Rdata", sep="")
WS <- paste(BASE, "Table1.Rdata", sep="")
#
#4 scenarios, Regular Case
# c(50, 10, 20, 100)
# c(10, 50, 100, 20)
# c(10, 10, 100, 100)
# c(50, 50, 20, 20)
BScenario <- list(c(50, 10),c(10, 50), c(10, 10),  c(50, 50))
mScenario <- list(c(20,100),c(100, 20),c(100, 100),c(20, 20))
#
# 4 scenarios, Tiny case. Use this only for testing.
# Comment out the next two lines, to produce actual example.
#BScenario <- list(c(5, 10),c(10, 5), c(10, 10),  c(50, 50))
#mScenario <- list(c(20,10),c(10, 20),c(10, 10),  c(2, 2))
#
ind <- as.list(rep(1:4, rep(NSIM, 4)))
#randomize
NumInd <- length(ind)
if (RandomizeQ) {
	RInd <- sample(1:NumInd, size=NumInd, replace=FALSE)
	OInd <- order(RInd)
	} else { RInd <- 1:NumInd
	OInd <- 1:NumInd
	}
ind <- ind[RInd]
#
#initialize functions
DAFNS <- function(X, y, Xt, yt){
 a1 <- ascrda(X=X, Y=y, Xt=Xt, Yt=yt)[[1]]
 a2 <- ascrda(X=X, Y=y, Xt=Xt, Yt=yt, SCRDAmethod = "ASCRDA1")[[1]]
 a3 <- ascrda(X=X, Y=y, Xt=Xt, Yt=yt, SCRDAmethod = "SCRDA")[[1]]
 a4 <- FitPAM(X=X, Y=y, Xt=Xt, Yt=yt)$Err
 a5 <- FitSvm(X=X, Y=y, Xt=Xt, Yt=yt)$r[2]
 a6 <- FitDLDA(X=X, Y=y, Xt=Xt, Yt=yt)$Err
 eta <- c(a1,a2,a3,a4,a5,a6)
 names(eta) <- c("ASCRDA","ASCRDA1","SCRDA","PAM","SVM","DLDA")
 eta
    }
cNames<-c("ASCRDA","ASCRDA1","SCRDA","PAM","SVM","DLDA")
nNames<-length(cNames)
#d=0.25 means that 25% of data used for test sample
OneIt <- function(iScenario){
    B <- BScenario[[iScenario]]
    m <- mScenario[[iScenario]]
    M <- synma(n=c(100,100), nt=c(500,500), rho = c(0.9, 0.9), B=B, m=m, fE=0.05)
    tryCatch(DAFNS(M$X, M$y, M$Xt, M$yt),error=function(e) rep(NA,6))
}

#
cl <- makeCluster(getOption("cl.cores", nWorkers))
clusterSetRNGStream(cl, 6789123)
clusterEvalQ(cl, require(ascrda))
clusterExport(cl, list("OneIt", "DAFNS", "BScenario","mScenario", "NSIM"))
startTime <- proc.time()[3]
date()
OUT <- parLapply(cl, X=ind, fun=OneIt)
stopCluster(cl)
date()
endTime <- proc.time()[3]
totalTime <- endTime - startTime
totalTime
#save workspace
save.image(WS)
#re-order output
OUT <-  OUT[OInd]
ind <- unlist(ind[OInd])
#list which ones are NA
WhichNA<-(1:length(OUT))[unlist(lapply(OUT, function(x) any(is.na(x))))]
NumNA <- length(WhichNA)
if (NumNA > 0) 
	cat(paste("number of NA's:", NumNA, fill=TRUE))
#WhichNA
#save raw output
save(OUT, file=OUTRaw)
#adjust raw output
rNames <- c("Scenario 1", "Scenario 2", "Scenario 3", "Scenario 4")
ans <- t(matrix(unlist(OUT), nrow=nNames))
colnames(ans) <- cNames
rownames(ans) <- ind
tb <- tbMed <- matrix(numeric(4*6), ncol=6)
dimnames(tbMed)<-dimnames(tb)<-list(rNames,cNames)
for (i in 1:4){ 
  tb[i,] <- 100*colMeans(ans[ind==i,],, na.rm=TRUE)
  tbMed[i,] <- 100*apply(X=ans[ind==i,], MARGIN=2, FUN=median, na.rm=TRUE)
}
tb
tbMed
#
write.table(tb, file=Outtb, row.names=FALSE)
#make latex table and output it
tbl <- xtable(tb,  digits=rep(1, ncol(tb)+1)) 	
print(tbl, include.rownames=FALSE, file=Outtbl)
#
write.table(tbMed, file=OuttbMed, row.names=FALSE)
#make latex table and output it
tblMed <- xtable(tbMed,  digits=rep(1, ncol(tb)+1)) 	
print(tblMed, include.rownames=FALSE, file=OuttblMed)
#
err <- as.vector(ans)
scenarios <- ordered(rep(paste("Scenario",ind), 6))
methods <-  factor(rep(cNames, rep(nrow(ans), 6)))
methods <- ordered(methods, levels=cNames)
ans.df <- data.frame(err=err, scenarios=scenarios, methods=methods)
bwplot(methods~err|scenarios, xlab=expression(eta), 
       panel=function(x,y){
       panel.grid(v=-5, h=0)
       panel.bwplot(x,y, notch=TRUE, fill="blue", lwd=2, coef=1.5, pch="|", cex=0.25) 
       			} 
       )
