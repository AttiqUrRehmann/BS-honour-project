# Mean square error (MSE) Comparison of the proposed, shrinkage of
# Schafer & Strimmer (2005) and
# maximum likelihood methods.
# Here Ture covariance matrix is random and target is Identity marix.

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

p <- n.var <- 10
n <- c(10, 20,30,40,50,60,70,90,120,200)

gamma <- rep(NA, length(n))
MSE <- matrix(NA, 3, length(n))

prop.non.zero <- 0.3
const <- 0.01
sigma <- un.str.cov(n.var, prop.non.zero, const) # Random covariance structure as a true covariance matrix

res <- replicate(1000,{
  for(i in 1:length(n)) {
    x <- mvrnorm(n=n[i], mu= rep(0,p), Sigma=sigma)
    S <- cor(x)
    
    # Calculate gamma values for the proposed estimator
    gamma[i] <- 1/(1+sum(abs(S - diag(p)))/(p))
    # Calculate proposed and shrinkage estimator 
    sigma.gamma <- gamma[i]*S + (1-gamma[i])*diag(p)
    sigma.shrink <- cor.shrink(x,verbose = FALSE)
    
    # Calculate MSE of all three competing procedures
    MSE_OF_MLE <- sum((S - sigma)^2)
    MSE_OF_SIGMA.GAMMA <- sum((sigma.gamma - sigma)^2)
    MSE_OF_SIGMA.SHRINK <- sum((sigma.shrink - sigma)^2)
    
    MSE[,i] <-  c(MSE_OF_MLE,MSE_OF_SIGMA.GAMMA,MSE_OF_SIGMA.SHRINK)
    
  }
  MSE
})
# Average the resulting MSE over 1000 simulations
ave_MSE <- apply(res, c(1,2), mean)
y.min <- min(ave_MSE)
y.max <- max(ave_MSE)

pdf(file="Ch04_Fig02(e).pdf")
par(mar=c(6,7,2,1), mgp = c(4,1,0))
plot(ave_MSE[1,], type="b",lwd=3,pch=1,ylim = c(y.min,y.max),
     cex.lab=2, ylab = "MSE", xlab = "",las=2, xaxt='n',
     cex.axis=2, col=20)
axis(1, at=c(1,2,3,4,5,6,7,8,9,10),
     labels = c(10, 20,30,40,50,60,70,90,120,200),
     cex.axis=1.2)
mtext("n",at=6, line = -31, cex= 2)
points(ave_MSE[2,], type="b", lwd=3,pch=1 ,col= "red")
points(ave_MSE[3,], type="b", lwd=3,pch=1 ,col= "green")
mtext("(e)", side = 3, line = 0.5, cex = 2)

legend("topright", legend = c("MLA", "Shrinkage", "MLE"),
       cex = 1.5,box.lwd = 1,box.lty = 1,
       text.font = 3,col = c("red", "green", 20),
       lwd=3, lty = 1, pch = 1)
dev.off()