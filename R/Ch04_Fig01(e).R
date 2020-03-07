# Comparison of the proposed, shrinkage of Schafer & Strimmer (2005) and
# maximum likelihood methods sum of absolute errors in estimated eigenvalues.
# Here Ture covariance matrix is random and target is identity

library(MASS)
library(corpcor)
# Function of random covariance structure 

un.str.cov <- function(n.var, prop.non.zero, const){ # The function to compute the
  AA <- matrix(0, n.var, n.var)      # random covariance structure of Schafer &         
  SS <- matrix(0, n.var, n.var)      # Strimmmer (2004) We will mention it in a
  # comment wherever it is required  
  AA[upper.tri(AA)] <- c(rbinom(n.var*(n.var-1)/2,1,prop.non.zero)) # proportion of non-zero elements
  
  BB <- AA + t(AA)
  for (i in 1:(n.var-1)){
    for (j in (i+1):n.var){
      if (AA[i,j]==1){
        SS[i,j] <- runif(1,-1,1) # Replace 1s with random values from the uniform distribution
      }
    }
  }
  
  SS <- SS + t(SS)
  ABS.SS <- abs(SS)
  ColSum <- apply(ABS.SS, 2, sum)
  
  for (i in 1:n.var){
    SS[i,i] <- ColSum[i] + const # fill diagonal entries with column sum plus small positive constant
  }
  
  cov2cor(SS)
}

p <- n.var <- c(30, 50 , 100)
n <- 50

gamma <- rep(NA,length(p))
ALL.EIGEN <- matrix(NA, length(p), length(p))

res <- replicate(1000,{
  for(i in 1:length(p)) {
    prop.non.zero <- 0.3
    const <- 0.01
    # random covariance structure as a ture covariance matrix 
    sigma <- un.str.cov(n.var[i], prop.non.zero, const)
    e.tre <- eigen(sigma)$values
    x <- mvrnorm(n=n, mu= rep(0,p[i]), Sigma=sigma)
    S <- cor(x)
    
    # Compute gamma values using Identity as a target
    gamma[i] <- 1/(1+sum(abs(S - diag(p[i])))/(p[i]))
    # Compute proposed estimator
    sigma.gamma <- gamma[i]*S + (1-gamma[i])*diag(p[i])
    
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

pdf(file="Ch04_Fig01(e).pdf")
par(mar=c(6,7,2,1), mgp = c(4,1,0))
boxplot(vec~arr, outline=FALSE, ylab= expression(sum(abs(hat(lambda[i]) - lambda[i]))/sum(lambda[i])),
        ylim=c(0,1), las=2,xlab = "" , xaxt="n",at=c(1,2,3,5,6,7,9,10,11),
        cex.axis=2, cex.lab=2, col=c("red","green", "blue"))
mtext("(e)", side = 3, line = 0.5, cex = 2)
mtext("30",at=2, line = -29, cex= 2)
mtext("50",at=6, line = -29, cex= 2)
mtext("100",at=10, line = -29, cex= 2)
mtext("p",at=6, line = -31, cex= 2)

legend("topleft", legend = c("MLA","Shrinkage","MLE"),cex = 1.5,box.lwd = 2,box.lty = 1,
       text.font = 3,col = c("red","green", "blue"), pch = 15)
dev.off()