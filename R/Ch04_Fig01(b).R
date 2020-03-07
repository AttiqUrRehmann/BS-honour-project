# Comparison of the proposed, shrinkage of Schafer & Strimmer (2005) and 
# maximum likelihood methods sum of absolute errors in estimated eigenvalues.
# Here the true covariance matrix is Identity matrix and target is AR(1)

library(MASS)
library(corpcor)

# t.ar1 function 

t.ar1 <-  function(x) { # The function that estimate "t" in case of AR(1)
  n <- nrow(x)          # covariance structure given in equation 3.19. 
  p <- ncol(x)          # We will mention it in a comment wherever it is
  cor.psi <- cor(x)     # required
  sum <- 0
  for(t in 2:p){
    sum <- sum + cor.psi[t,t-1]
  }
  hat.t <- sum/(p-1)
}

p <- c(30, 50 , 100)
n <- 50
t <- 0.5

gamma <- rep(NA,length(p))
ALL.EIGEN <- matrix(NA, length(p), length(p))

res <- replicate(1000,{
  for(i in 1:length(p)) {
    sigma <- diag(p[i]) # Identity as a true covariance matrix
    e.tre <- eigen(sigma)$values
    x <- mvrnorm(n=n, mu= rep(0,p[i]), Sigma=sigma)
    # AR(1) as a target matrix
    # compute target matrix
    t.hat <- t.ar1(x)
    TAR.AR1 <- t.hat ^ outer(1:p[i], 1:p[i], function(aa, bb) abs(aa - bb))
    S <- cor(x)
    k <- seq(0,p[i]-1)
    # Gamma values for the proposed method
    gamma[i] <- 1/(1+sum(abs(S - TAR.AR1))/(p[i]+sum(2*k*(t.hat)^(p[i]-k))))
    # Compute proposed estimator
    sigma.gamma <- gamma[i]*cor(x) + (1-gamma[i])*TAR.AR1
    
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

pdf(file="Ch04_Fig01(b).pdf")
par(mar=c(6,7,2,1), mgp = c(4,1,0))
boxplot(vec~arr, outline=FALSE, ylab= expression(sum(abs(hat(lambda[i]) - lambda[i]))/sum(lambda[i])),
        ylim=c(0,1.2), las=2, xlab = "", xaxt="n",at=c(1,2,3,5,6,7,9,10,11),
        cex.axis=2, cex.lab=2, col=c("red","green", "blue"))
mtext("(b)", side = 3, line = 0.5, cex = 2)
mtext("30",at=2, line = -29, cex= 2)
mtext("50",at=6, line = -29, cex= 2)
mtext("100",at=10, line = -29, cex= 2)
mtext("p",at=6, line = -31, cex= 2)

legend("topleft", legend = c("MLA","Shrinkage","MLE"),cex = 1.5,box.lwd = 2,box.lty = 1,
       text.font = 3,col = c("red","green", "blue"), pch = 15)
dev.off()