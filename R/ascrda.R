ascrda <-
function(X, y, Xt, yt, k="default", alpha=seq(0,0.99,0.11), delta=seq(0,3,0.2), ...){
#
#Misclassification rate for train and test data using RDA
#
#INPUT
# X        : Training sample, expression matrix
# Xt       : Test sample, expression matrix
# y        : Class in training data: y should be 1, 2...
# yt       : Class in test data. 
# k        : "default",  use pseudolikelihood, k>0 set k in kNN
# alpha    : alpha vector in function 'rda'
# delta    : delta vector in function 'rda'
#
#OUTPUT
#mis-classification rates (test data only) for
#ASCRDA, SCRDA, SVM, SCRDA/SCRDA, SVM/SCRDA
#
#input validation
    stopifnot(is.vector(y), is.vector(yt), 
                nrow(X)==length(y), 
                nrow(Xt)==length(yt),
                ncol(Xt)==ncol(X),
                is.character(k)|is.numeric(k))
    is.wholenumber <-
        function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
    if (is.numeric(k))
        stopifnot(is.wholenumber(k), k>=1)
    #
    #output vector labels
    resNames <- c("ASCRDA", "SCRDA", "SVM", "SCRDA/SCRDA", "SVM/SCRDA")
    #
    #SCRDA
    #we don't use Fitrda since we need gene selection
    fit <- rda(t(X), y) 
    outRDA <- try(capture.output(
                 fit.cv <-rda.cv(fit,t(X),y,alpha=alpha,delta=delta,...)),TRUE)
    if (class(outRDA) == "try-error") return(-1)
    a<-fit.cv$cv.err
    b<-fit.cv$ngene
    La<-length(alpha)
    Ld<-length(delta)
    IndDel<-rep(delta,each=La)
    IndAlp<-rep(alpha,Ld)
# find ALL the (alpha, delta) that produce(s) minimum cv.err
    Ia <- 1:length(a)
    NInd <- Ia[a==min(a)] 
    Ngene <- b[NInd]  
    OptInd <- NInd[Ngene==min(Ngene)]
    Alpha <- IndAlp[OptInd]
    Delta <- IndDel[OptInd] 
    nAlpha <- length(Alpha)
  ytp<- sapply(1:nAlpha, function(i) predict(fit, x=t(X), y=y, xnew=t(Xt), alpha=Alpha[i], delta=Delta[i]))   
  if(length(yt)==1) 
    {MisclassTest<- ytp!=yt}  
  else
    {MisclassTest<-apply(ytp,2,function(x) sum(yt!=x)/length(yt))} # we do so in case test sample has one row
  TestRDA <- mean(MisclassTest) 
#
#feature selection methods  
glist<-as.numeric(genelist.rda(t(X), y, alpha=Alpha[1], delta=Delta[1]))
Xnew<-X[,glist]
Xtnew<-Xt[,glist]
#
#SCRDA/SCRDA
TestSCRDA2 <- FitRda(Xnew, y, Xtnew, yt)[2]
#  
#SVM and SVM/SCRDA
  SVMSCRDATest <- FitSVM(Xnew, y, Xtnew, yt)[2]
  SVMTestError <- FitSVM(X, y, Xt, yt)[2]
  if(is.character(k)) 
      K <- khat(Xnew, y, plot=FALSE) 
    else
      K <- k
  XZ<-cbind(Xnew,nnc(Xnew,y, K))
  XZt<-cbind(Xtnew,nncTest(Xnew,y,Xtnew, K))
  Testnnc<-FitRda(XZ,y, XZt,yt)[2]
# cat(" optimal k =", K, fill=TRUE) 
# resNames <- c("ASCRDA", "SCRDA", "SVM", "SCRDA/SCRDA", "SVM/SCRDA")
  res <- c(Testnnc, TestRDA, SVMTestError,  TestSCRDA2, SVMSCRDATest)
  names(res) <- resNames
  res      
} 


