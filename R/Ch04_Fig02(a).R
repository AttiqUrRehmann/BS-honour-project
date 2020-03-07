# Mean square error (MSE) Comparison of the proposed, shrinkage of
# Schafer & Strimmer (2005) and
# maximum likelihood methods sum of absolute errors in estimated eigenvalues.
# Here Ture covariance matrix and target is AR(1)

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

p <- n.var <- 10
n <- c(10, 20,30,40,50,60,70,90,120,200)
t <- 0.5
sigma <-  t ^ outer(1:p, 1:p, function(aa, bb) abs(aa - bb)) # true covariance matrix AR(1)
gamma <- rep(NA, length(n))
MSE <- matrix(NA, 3, length(n)) # Null matrix for the MSE of all three competing procedures to be store in it and use for further analysis

res <- replicate(1000,{
  for(i in 1:length(n)) {
    x <- mvrnorm(n=n[i], mu= rep(0,p), Sigma=sigma)
    # Compute target matrix, i.e, AR(1)
    t.hat <- t.ar1(x)
    TAR.AR1 <- t.hat ^ outer(1:p, 1:p, function(aa, bb) abs(aa - bb))
    S <- cor(x)
    
    # Compute gamma for the proposed estimator
    k <- seq(0,p-1)
    gamma[i] <- 1/(1+sum(abs(S - TAR.AR1))/(p+sum(2*k*(t.hat)^(p-k))))
    # Compute proposed estimator of the covariance matrix
    sigma.gamma = gamma[i]*S + (1-gamma[i])*TAR.AR1
    # Compute shrinkage estimator of the covariance matrix
    sigma.shrink <- cor.shrink(x,verbose = FALSE)
    
    # Compute the MSE of all three estimators
    MSE_OF_MLE <- sum((S - sigma)^2)
    MSE_OF_SIGMA.GAMMA <- sum((sigma.gamma - sigma)^2)
    MSE_OF_SIGMA.SHRINK <- sum((sigma.shrink - sigma)^2)
    
    # Store the resulting MSE in the null matrix
    MSE[,i] <- c(MSE_OF_MLE,MSE_OF_SIGMA.GAMMA,MSE_OF_SIGMA.SHRINK)
  }
  MSE
})
# Average the MSE repeated 1000 times
ave_MSE <- apply(res, c(1,2), mean)

y.min <- min(ave_MSE)
y.max <- max(ave_MSE)

pdf(file="Ch04_Fig02(a).pdf")
par(mar=c(6,6,2,1), mgp = c(4,1,0))
plot(ave_MSE[1,], type="b",lwd=3,pch=1,ylim = c(y.min,y.max),
     cex.lab=2, ylab = "MSE", xlab = "",las=2, xaxt='n',
     cex.axis=2, col=20)
axis(1, at=c(1,2,3,4,5,6,7,8,9,10),
     labels = c(10, 20,30,40,50,60,70,90,120,200),
     cex.axis=1.3)
mtext("n",at=6, line = -31, cex= 2)
points(ave_MSE[2,], type="b", lwd=3,pch=1 ,col= "red")
points(ave_MSE[3,], type="b", lwd=3,pch=1 ,col= "green")
mtext("(a)", side = 3, line = 0.5, cex = 2)

legend("topright", legend = c("MLA", "Shrinkage", "MLE"),
       cex = 1.5,box.lwd = 1,box.lty = 1,
       text.font = 3,col = c("red", "green", 20),
       lwd=3, lty = 1, pch = 1)
dev.off()