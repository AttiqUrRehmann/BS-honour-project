# Comparison of the proposed, shrinkage of Schafer & Strimmer (2005) and
# maximum likelihood methods sum of absolute errors in estimated eigenvalues.
# Here the true covariance matrix and target is echangeable.

library(MASS)
library(corpcor)

# t.exch function

t.exch <-  function(x) { # The function that estimate "t" in case of
  n <- nrow(x)           # exchangeable covariance structure given in
  p <- ncol(x)           # equation 3.13. We will mention it in a 
  cc <- cor(x)           # comment wherever it is required
  sum <- 0
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      sum <- sum + cc[i,j]
    }
  }
  hat.t <- 2*sum/(p*(p-1))
}

p <- c(30, 50 , 100)
n <- 50
t <- 0.5

gamma <- rep(NA,length(p))
ALL.EIGEN <- matrix(NA, length(p), length(p))

res <- replicate(1000,{
  
  for(i in 1:length(p)) {
    
    sigma <-  matrix(t, p[i], p[i])
    diag(sigma) = 1
    e.tre <- eigen(sigma)$values
    x <- mvrnorm(n=n, mu= rep(0,p[i]), Sigma=sigma)
    # Exchangeable as a target matrix
    # compute target matrix
    t.hat <- t.exch(x)
    TAR.EXCH <- matrix(t.hat, p[i], p[i])
    diag(TAR.EXCH) = 1
    S <- cor(x)
    
    # Gamma values for the proposed method
    gamma[i] <- 1/(1+sum(abs(S - TAR.EXCH))/(p[i] + (p[i]*(p[i]-1)*t.hat)))
    # Compute proposed estimator
    sigma.gamma <- gamma[i]*S + (1-gamma[i])*TAR.EXCH
    
    # Compute eigenvalues of all three competing estimators
    MLA.eigen_values <- eigen(sigma.gamma)$values
    shrink.eigen_values <- eigen(cor.shrink(x,verbose = FALSE))$values
    MLE.eigen_values <- eigen(S)$values
    
    # Compute sum of absolute errors in estimated eigenvalues of all three competing estimators
    sum.eigen.MLA <-  sum(abs(MLA.eigen_values - e.tre ))/sum(e.tre)
    sum.eigen.shrink <- sum(abs(shrink.eigen_values - e.tre ))/sum(e.tre)
    sum.eigen.MLE <- sum(abs(MLE.eigen_values - e.tre ))/sum(e.tre)
    
    ALL.EIGEN[,i] <- c(sum.eigen.MLA, sum.eigen.shrink, sum.eigen.MLE)
  }
  ALL.EIGEN
})

vec <- as.vector(res)
arr <- as.vector(array(c(1:9),dim = c(3,3,1000)))

pdf(file="Ch04_Fig01(c).pdf")
par(mar=c(6,7,2,1), mgp = c(4,1,0))
boxplot(vec~arr, outline=FALSE, ylab= expression(sum(abs(hat(lambda[i]) - lambda[i]))/sum(lambda[i])),
        ylim=c(0,1), las=2, xlab = "" ,xaxt="n",at=c(1,2,3,5,6,7,9,10,11),
        cex.axis=2, cex.lab=2, col=c("red","green", "blue"))
mtext("(c)", side = 3, line = 0.5, cex = 2)
mtext("30",at=2, line = -29, cex= 2)
mtext("50",at=6, line = -29, cex= 2)
mtext("100",at=10, line = -29, cex= 2)
mtext("p",at=6, line = -31, cex= 2)

legend("topleft", legend = c("MLA","Shrinkage","MLE"),cex = 1.5,box.lwd = 2,box.lty = 1,
       text.font = 3,col = c("red","green", "blue"), pch = 15)
dev.off()