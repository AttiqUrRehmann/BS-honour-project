# Eigenvalues Comparison of the proposed, shrinkage of
# Schafer & Strimmer (2005) and maximum likelihood methods.
# Here Ture covariance matrix is Identity and target is AR(1).

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

p <- n.var <- 50
n <- 30
t <- 0.5

sigma <-  diag(p) # Identity as a true covariance matrix
e.tre <- eigen(sigma)$values
eigen.matrix <- matrix(NA,3,n.var)

res <- replicate(1000,{
  x <- mvrnorm(n=n, mu= rep(0,p), Sigma=sigma)
  # Calculate target matrix, i.e, AR(1)
  t.hat <- t.ar1(x)
  TAR.AR1 <- t.hat ^ outer(1:p, 1:p, function(aa, bb) abs(aa - bb))
  S <- cor(x)
  
  # Calculate gamma values for the proposed method
  k <- seq(0,p-1)
  gamma <- 1/(1+sum(abs(S - TAR.AR1))/(p+sum(2*k*(t.hat)^(p-k))))
  # Calculate the proposed estimator of the true covariance matrix
  sigma.gamma = gamma*S + (1-gamma)*TAR.AR1
  
  # Calculate the eigenvalues of all three proposed estimators
  MLA.eigen_values <- eigen(sigma.gamma)$values
  eigen_values.shrink <- eigen(cor.shrink(x,verbose = FALSE))$values
  MLE.eigen_values <- eigen(S)$values
  
  # Store the eigenvalues in the null matrix
  eigen.matrix[1,] <- MLE.eigen_values
  eigen.matrix[2,] <- MLA.eigen_values
  eigen.matrix[3,] <- eigen_values.shrink
  eigen.matrix
})
# Average the resulting eigenvalue repeated 1000 times
ave_eigen <- apply(res, c(1,2), mean)
yl <- max(rbind(ave_eigen, e.tre))

pdf(file="Ch04_Fig03(b).pdf")
par(mar=c(7,6,2,1), mgp = c(4,1,0))
plot(e.tre, type="b",lwd=3,pch=16,ylim = c(0,yl), xlab="Order",
     ylab = "Eigenvalues", cex.axis=2, yaxt="n", cex.lab=2)
axis(2, las=2, cex.axis=2)
points(ave_eigen[1,], type="b", lwd=3,pch=16 ,col= "20")
points(ave_eigen[2,], type="b", lwd=3,pch=16 ,col= "red")
points(ave_eigen[3,], type="b", lwd=3,pch=16 ,col= "green")
mtext("(b)", side = 3, line = 0.5, cex = 2)

legend("topright", legend = c("True", "MLA", "Shrinkage", "MLE"),
       cex = 1.5,box.lwd = 1,box.lty = 1,
       text.font = 3,col = c("black","red", "green", 20),
       lty = 1, pch = 16)
dev.off()