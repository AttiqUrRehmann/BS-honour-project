# Eigenvalues Comparison of the proposed, shrinkage of
# Schafer & Strimmer (2005) and maximum likelihood methods.
# Here Ture covariance matrix and target is Exchangeable.


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

p <- n.var <- 50
n <- 30
t <- 0.3
sigma <-  matrix(t, p, p) # Exchangeable covariance structure as a true covariance matrix
diag(sigma) = 1
e.tre <- eigen(sigma)$values

eigen.matrix <- matrix(NA,3,n.var)
x <- mvrnorm(n=n, mu= rep(0,p), Sigma=sigma)

res <- replicate(1000,{
  # Calculate the target matrix, i.e, exchangeable
  t.hat <- t.exch(x)
  TAR.EXCH <- matrix(t.hat, p, p)
  diag(TAR.EXCH) = 1
  S <- cor(x)
  
  # Calculate the gamma values for the proposed method
  gamma <- 1/(1+sum(abs(S - TAR.EXCH))/(p + (p*(p-1)*t.hat)))
  # Calculate the proposed estimator 
  sigma.gamma = gamma*S + (1-gamma)*TAR.EXCH
  
  # calculate the eigenvalues of all three competing procedures
  MLA.eigen_values <- eigen(sigma.gamma)$values
  eigen_values.shrink <- eigen(cor.shrink(x,verbose = FALSE))$values
  MLE.eigen_values <- eigen(S)$values
  
  # Store the resulting eigenvalues in the null matrix 
  eigen.matrix[1,] <- MLE.eigen_values
  eigen.matrix[2,] <- MLA.eigen_values
  eigen.matrix[3,] <- eigen_values.shrink
  eigen.matrix
})
# Average the resulting eigenvalues repeated 1000 times
ave_eigen <- apply(res, c(1,2), mean)
yl <- max(rbind(ave_eigen, e.tre))

pdf(file="Ch04_Fig03(c).pdf")
par(mar=c(7,6,2,1), mgp = c(4,1,0))
plot(e.tre, type="b",lwd=3,pch=16,ylim = c(0,yl), xlab="Order",
     ylab = "Eigenvalues", cex.axis=2, yaxt="n", cex.lab=2)
axis(2, las=2, cex.axis=2)
points(ave_eigen[1,], type="b", lwd=3,pch=16 ,col= "20")
points(ave_eigen[2,], type="b", lwd=3,pch=16 ,col= "red")
points(ave_eigen[3,], type="b", lwd=3,pch=16 ,col= "green")
mtext("(c)", side = 3, line = 0.5, cex = 2)

legend("topright", legend = c("True", "MLA", "Shrinkage", "MLE"),
       cex = 1.5,box.lwd = 1,box.lty = 1,
       text.font = 3,col = c("black","red", "green", 20),
       lty = 1, pch = 16)
dev.off()