# Eigenvalues Comparison of the proposed, shrinkage of
# Schafer & Strimmer (2005) and maximum likelihood methods.
# Here Ture covariance matrix is random and target is Identity.


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

n.var <- 50
n <- 30
prop.non.zero <- 0.30
const <- 0.01
sigma <- un.str.cov(n.var, prop.non.zero, const) # random covariance structure as a true covariance matrix

e.tre <- eigen(sigma)$values
eigen.matrix <- matrix(NA,3,n.var)

res <- replicate(1000,{
  x <- mvrnorm(n=n, mu= rep(0,n.var), Sigma=sigma)
  S <- cor(x)
  
  # Calculate gamma values for the proposed method
  gamma <- 1/(1+sum(abs(S - diag(n.var)))/(n.var))
  # Calculate teh proposed estimator of the true covariance matrix
  sigma.gamma = gamma*S + (1-gamma)*diag(n.var)
  
  # Calculate the eigenvalues of all three competing estimators
  MLA.eigen_values <- eigen(sigma.gamma)$values
  eigen_values.shrink <- eigen(cor.shrink(x,verbose = FALSE))$values
  MLE.eigen_values <- eigen(S)$values
  
  eigen.matrix[1,] <- MLE.eigen_values
  eigen.matrix[2,] <- MLA.eigen_values
  eigen.matrix[3,] <- eigen_values.shrink
  eigen.matrix
})
# Average the resulting eigenvalues simulated 1000 times
ave_eigen <- apply(res, c(1,2), mean)
yl <- max(rbind(ave_eigen, e.tre))

pdf(file="Ch04_Fig03(e).pdf")
par(mar=c(7,6,2,1), mgp = c(4,1,0))
plot(e.tre, type="b",lwd=3,pch=16,ylim = c(0,yl), xlab="Order",
     ylab = "Eigenvalues", cex.axis=2, yaxt="n", cex.lab=2)
axis(2, las=2, cex.axis=2)
points(ave_eigen[1,], type="b", lwd=3,pch=16 ,col= "20")
points(ave_eigen[2,], type="b", lwd=3,pch=16 ,col= "red")
points(ave_eigen[3,], type="b", lwd=3,pch=16 ,col= "green")
mtext("(e)", side = 3, line = 0.5, cex = 2)

legend("topright", legend = c("True", "MLA", "Shrinkage", "MLE"),
       cex = 1.5,box.lwd = 1,box.lty = 1,
       text.font = 3,col = c("black","red", "green", 20),
       lty = 1, pch = 16)
dev.off()