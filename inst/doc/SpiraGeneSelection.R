library("ascrda")
X <- Spira$X
y <- as.numeric(Spira$Patients$STATUS)
X <- scale(t(X))  #row standardization
glist<-as.numeric(genelist.rda(X, y, alpha=0.3, delta=0.1))
length(glist)
