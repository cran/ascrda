rdaSVM <-
function(X, y, Xt, yt, alpha=seq(0,0.99,0.11), delta=seq(0,3,0.2), ...){
# We select genes using RDA before using SVM
# Misclassification rate for train and test data using SVMnnc
#
# X        : Training sample, expression matrix
# Xt       : Test sample, expression matrix
# y        : Class in training data: y should be 1, 2...
# yt       : Class in test data.
# k        : "default",  use pseudolikelihood, k>0 set k in kNN
# alpha    : alpha vector in function 'rda'
# delta    : delta vector in function 'rda'
#
	stopifnot(is.vector(y), is.vector(yt),
				nrow(X)==length(y),
				nrow(Xt)==length(yt),
	 			ncol(Xt)==ncol(X))

	fit <- rda(t(X), y)
	sink("Junk.txt")
	tryCatch(
     		fit.cv <- rda.cv(fit,t(X),y,alpha=alpha,delta=delta,...),
		    finally = {
					 sink()
					 unlink("Junk.txt")
					 })
	a<-fit.cv$cv.err
	b<-fit.cv$ngene
	La<-length(alpha)
	Ld<-length(delta)
	IndDel<-rep(delta,each=La)
	IndAlp<-rep(alpha,Ld)
# find ALL the (alpha, delta) that produce(s) minimum cv.err
	Ia <- 1:length(a)
	NInd<-Ia[a==min(a)]
	Ngene<-b[NInd]
	OptInd<-NInd[Ngene==min(Ngene)]
	Alpha<-IndAlp[OptInd]
	Delta<-IndDel[OptInd]
	nAlpha<-length(Alpha)

  glist<-as.numeric(genelist.rda(t(X), y, alpha=Alpha[1], delta=Delta[1]))
  Xnew<-X[,glist]
  Xtnew<-Xt[,glist]
  FitSVM(Xnew,y, Xtnew,yt)
}
