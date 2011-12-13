ascrda <-
function(X, Y, Xt, Yt, alpha=seq(0,0.99,0.11), delta=seq(0,3,0.2), 
	SCRDAmethod=c("ASCRDAk","ASCRDA1","SCRDA"), ...){
#
#Misclassification rate for train and test data using RDA
#
#INPUT
# X        : Training sample, data matrix, n-by-G
# Xt       : Test sample
# Y        : Class in training data: y should be 1, 2...
# Yt       : Class in test data. 
# alpha    : alpha vector in function 'rda'
# delta    : delta vector in function 'rda'
#
#OUTPUT
#mis-classification rate
#
#ensure classes are 1, 2, ... etc
if (is.data.frame(Y)) y<-as.matrix.data.frame(Y) else y<-Y
if (is.data.frame(Yt)) yt<-as.matrix.data.frame(Yt) else yt<-Yt
y<-as.numeric(ordered(y))
yt<-as.numeric(ordered(yt))
#
#input validation
    stopifnot(is.vector(y), is.vector(yt), 
                nrow(X)==length(y), 
                nrow(Xt)==length(yt),
                ncol(Xt)==ncol(X)
		)
whichMethod <- match.arg(SCRDAmethod)
res <- switch(whichMethod,
	SCRDA = {
		Err <- FitRda(X,y, Xt,yt)[2]
		list(Err = Err, k=0)
	},
	ASCRDA1 = {  	
		XAC <- nnc(X, y, k=1)
  		XZ <-  cbind(X, XAC)
		XACt <- nncTest(X, y, Xt, k=1)
		XZt<- cbind(Xt,XACt)
		Err <- FitRda(XZ,y, XZt,yt)[2]
		list(Err = Err, k = 1)
	},
	ASCRDAk = {
		KHAT <- khat(X, y, plot=FALSE) 
		XAC<-nnc(X, y, k=KHAT)
		XZ<-cbind(X, XAC)
		XACt <- nncTest(X, y, Xt, k=KHAT)
		XZt<-cbind(Xt,XACt)
		Err <-FitRda(XZ,y, XZt,yt)[2]
		list(Err = Err, k = KHAT)
	}
)
#OUTPUT
resNames <- c(as.character(whichMethod), "K")
names(res) <- resNames
res      
} 


