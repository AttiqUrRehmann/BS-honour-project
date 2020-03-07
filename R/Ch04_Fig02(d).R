# Mean square error (MSE) Comparison of the proposed, shrinkage of
# Schafer & Strimmer (2005) and
# maximum likelihood methods sum of absolute errors in estimated eigenvalues.
# Here Ture covariance matrix Identity and target is exchangeable.

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

p <- n.var <- 10
n <- c(10, 20,30,40,50,60,70,90,120,200)
t <- 0.5
sigma <-  diag(p) # Identity as a true covariance matrix
gamma <- rep(NA, length(n))
MSE <- matrix(NA, 3, length(n))

res <- replicate(1000,{
  for(i in 1:length(n)) {
    x <- mvrnorm(n=n[i], mu= rep(0,p), Sigma=sigma)
    # Compute target matrix, i.e, exchangeable
    t.hat <- t.exch(x)
    TAR.EXCH <- matrix(t.hat, p, p)
    diag(TAR.EXCH) = 1
    S <- cor(x)
    
    # Calculate gamma values for the proposed estimator
    gamma[i] <- 1/(1+sum(abs(S - TAR.EXCH))/(p + (p*(p-1)*t.hat)))
    # Calculate proposed and shrinkage estimator
    sigma.gamma <- gamma[i]*S + (1-gamma[i])*TAR.EXCH
    sigma.shrink <- cor.shrink(x,verbose = FALSE)
    
    # Calculate MSE of all three competing procedures
    MSE_OF_MLE <- sum((S - sigma)^2)
    MSE_OF_SIGMA.GAMMA <- sum((sigma.gamma - sigma)^2)
    MSE_OF_SIGMA.SHRINK <- sum((sigma.shrink - sigma)^2)
    
    MSE[,i] <- c(MSE_OF_MLE,MSE_OF_SIGMA.GAMMA,MSE_OF_SIGMA.SHRINK)
  }
  MSE
})
# Average the resulting MSE over 1000 simulations
ave_MSE <- apply(res, c(1,2), mean)
y.min <- min(ave_MSE)
y.max <- max(ave_MSE)

pdf(file="Ch04_Fig02(d).pdf")
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
mtext("(d)", side = 3, line = 0.5, cex = 2)

legend("topright", legend = c("MLA", "Shrinkage", "MLE"),
       cex = 1.5,box.lwd = 1,box.lty = 1,
       text.font = 3,col = c("red", "green", 20),
       lwd=3, lty = 1, pch = 1)
dev.off()