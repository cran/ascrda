FitPAM <-
function(X, Y, Xt, Yt, PamMethod=c("default", "kNN1", "kNN")){
ok1 <- is.matrix(X)&&is.matrix(Xt)
ok2 <- nrow(X)==length(Y)&&nrow(Xt)==length(Yt) #check in design matrix form
uy <- sort(unique(Y)) #sorted unique y's
uyt <- sort(unique(Yt)) #sorted unique y's
ok3 <- all(uy==seq(length(uy)))
ok4 <- all(uyt==seq(length(uyt)))
ok <- ok1&&ok2&&ok3&&ok4
if (!ok) stop("FitPAM - Error code:", c("ok1","ok2","ok3","ok4")[!c(ok1,ok2,ok3,ok4)])
##
    whichMethod <- match.arg(PamMethod)
    y <- Y
    yt <- Yt
#using expression matrix form
    X <- t(X)
    Xt <- t(Xt)
    M <- switch(whichMethod,
	   default = {
		list(X=X, Xt=Xt)
	           },
	   kNN1 = {
            XAC <- nnc(t(X), y, k=1)
            x.adj <- pamr.decorrelate(X, as.data.frame.matrix(XAC))$x.adj
            XACt <- nncTest(t(X), y, t(Xt), k=1)
            xt.adj <- pamr.decorrelate(Xt,  as.data.frame.matrix(XACt))$x.adj
            list(X=x.adj, Xt=xt.adj)
	       },
	   kNN = {
            KHAT <- khat(t(X), y, plot=FALSE)
            XAC <- nnc(t(X), y, k=KHAT)
            x.adj <- pamr.decorrelate(X, as.data.frame.matrix(XAC))$x.adj
            XACt <- nncTest(t(X), y, t(Xt), k=KHAT)
            xt.adj <- pamr.decorrelate(Xt, as.data.frame.matrix(XACt))$x.adj
            list(X=x.adj, Xt=xt.adj)
	       } )
    X <- M$X
    Xt <- M$Xt
    capture.output(
        {mydata <<- list(x=X, y=as.factor(y))
         mytrain <-   pamr.train(mydata)
         new.scales <- pamr.adaptthresh(mytrain)
         mytrain2 <- pamr.train(mydata, threshold.scale=new.scales)
         }
    )
    ytp <- pamr.predict(mytrain2, Xt, threshold=new.scales)
    ytp <- as.numeric(ytp)
    cmatTest <- table(yt, ytp, dnn = c("Observed", "Predicted"))
    Err <- 1 - sum(diag(cmatTest))/sum(cmatTest)
    list(Err = Err, confusionMatrix = cmatTest)
    }
