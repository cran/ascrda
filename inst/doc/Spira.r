#Source: Spira.R
#Note: this is a large dataset, the full dataset can be run PC with 24 GB RAM
#  but with less memory it may crash. The bottleneck is with the SVM code.
#comparing: 
#           "ASCRDA","ASCRDA1","SCRDA","kNN", "SVM", "DLDA"
#load required libraries
require("ascrda")
require("parallel")
require("xtable")
#
#System Requirements
#requires package snow for multicore computing
#multicore computer. See below CPU CUSTOMIZATION - number of cores
#
#DATASET CUSTOMIZATION
#customize for dataset. Will work for other X and y
X <- Spira$X
y <- as.numeric(Spira$Patients$STATUS)
X <- t(scale(t(X)))  #row standardization

#Spira DATA. smokers: current, former and never
#Summary of performance
#> tb
#       ASCRDA ASCRDA1 SCRDA  PAM  SVM DLDA
#mean     37.2    37.5  36.9 36.9 47.6 53.5
#median   38.9    38.9  38.9 38.9 44.4 55.6 
#NumIter = 500 took about 10.2 hrs on unit.uwo.ca
#13 evaluations failed


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
OneIt <- function(it){
    M <- getTrainTest(X, y, d=0.25)
    tryCatch(DAFNS(M$X, M$y, M$Xt, M$yt),error=function(e) rep(NA,6))
}
#
#Optional: time one iteration
#about 37 sec
#startTime <- proc.time()[3]
#OneIt(1)
#endTime <- proc.time()[3]
#endTime-startTime

TextOutput <- "TableSpira3Class.txt" 
TexOutput  <- "TableSpira3Class.tex"
ansName    <- "ansSpira3Class.Rdata"
Title <- "Spira dataset with 3 classes"
WS <- "Spira.Rdata"
#
#where n is the number of samples and G is the number of genes

#OUTPUT FILES. 
#User dependent, where to put the output?
#Set base name - must end in /
#This is system dependent - PC, Mac or linux
#Windows version:
if (.Platform$OS.type=="windows") BASE <- "d:/r/2011/ascrda/OUTSpira/" else BASE <- "/Users/aim/R/2011/ascrda/ascrda/inst/doc/"
#Names of files to save - choose different names for different datasets!
OutFileTxt  <- paste(BASE, TextOutput, sep="")
OutFileTex  <- paste(BASE, TexOutput, sep="")
ansFile  <- paste(BASE, ansName, sep="")
WS <- paste(BASE, WS, sep="")
#
#CPU CUSTOMIZATION
#Using intel i7 with 6 cores and 12 compute nodes
#takes 1237 sec with 12 nodes & 12 NumIter
nWorkers<-12
#
#Number of cross-validation iteratios, recommended NumIter<-1000
NumIter <- 500
Iter <- as.list(1:NumIter)
#
cl <- makeCluster(getOption("cl.cores", nWorkers))
clusterSetRNGStream(cl, 6789123)
clusterEvalQ(cl, require(ascrda))
clusterExport(cl, list("X", "y", "OneIt", "DAFNS"))
startTime <- proc.time()[3]
date()
OUT <- parLapply(cl, X=Iter, fun=OneIt)
stopCluster(cl)
date()
endTime <- proc.time()[3]
totalTime <- endTime - startTime
totalTime

#save ws
save.image(WS)
#normally this does not happen!
if (any(is.na(unlist(OUT))) ) {
	OUT2<-vector("list", NumIter)
	cat("Some evaluations failed - removing them..", fill=TRUE)
	j <- 0
	for (i in 1:length(OUT)) {
		if (any(is.na(OUT[[i]]))) next else{
			j <- j+1
			OUT2[[j]] <- OUT[[i]]
			}
		}
	OUT <- OUT2[1:j]
	cat(paste("Removed:", NumIter-j),fill=TRUE)
}
#
#adjust raw output and make table
ans <- t(matrix(unlist(OUT), nrow=nNames))
dimnames(ans)[[2]]<-cNames
save(ans, file=ansFile)
tb <- matrix(colMeans(ans), nrow=1)
tb <- round(100*tb,1) # use percentages
tbMed <- round(100*matrix(apply(ans, 2, median), nrow=1),1)
tb<-rbind(tb, tbMed)
dimnames(tb)[[1]]<- c("mean", "median")
dimnames(tb)[[2]]<- cNames
tb
#
write.table(tb, file=OutFileTxt, row.names=FALSE)
#latex table
tbl <- xtable(tb,  digits=rep(1, ncol(tb)+1), caption=Title) 
print(tbl, include.rownames=FALSE, file=OutFileTex)
#boxplot
require(lattice)
methodsFactor <- rep(dimnames(tb)[[2]], rep(length(OUT), nNames))
methodsFactor <- ordered(methodsFactor, levels=cNames)
ans.df <- data.frame(err=as.vector(ans), method=methodsFactor)
graphics.off()
win.graph()
bwplot(method~err, data=ans.df, xlab=expression(eta), panel=function(x,y){
	panel.bwplot(x, y, horizontal=TRUE, fill="lightblue")
	panel.grid(v=-1, h=0)
})


