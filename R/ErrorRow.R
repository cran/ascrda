ErrorRow <-
function(r, Data, cl, FUN, ...){
 #
 # Finding training and test test error rate where the test data is
 # r-th row. Number of variables selected beforehand.
 #
 # r    : The row numbers which correspond to test data
 # Data : Combination of training and test data
 # cl   : Combination of training and test class
 # P    : Number of genes to be selected initially
 # FUN  : Methods to be implemented
 #
            Xtr<-Data[-r,]
            ytr<-cl[-r]
            yte<-cl[r]
            rownames(Xtr)<-1:nrow(Xtr)
            Xte = as.matrix(Data[r,], nrow=length(r))
            rownames(Xte)<-1:nrow(Xte)
            colnames(Xte)<-colnames(Xtr)
            out <- try( FUN(Xtr, ytr, Xte, yte, ...), TRUE)
            if(class(out)=="try-error") out <- -1
            out
 }

