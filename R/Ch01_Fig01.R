# To see what happens if we increase number of variables 

library(MASS) 

n <- c(25,50,100,1000)
p <- 50


mat <- matrix(NA, nrow=length(n), ncol=p)

# True covariance and its eigenvalues
sigma <- diag(1, nrow = p, ncol = p)
e.tre <- eigen(sigma)$values

# Calculating sample covariance matrix eigenvalues for different samples "n"
for(i in 1:length(n)){
  # Repeat the results 1000 times. We will repeat the results for all Figures in R codes. 
  res <- replicate(1000, {
    # Generating "n * p" X matrix from multivariate normal distribution 
    x <- mvrnorm(n=n[i], mu=rep(0,p), Sigma=sigma)
    # Sample covarince matrix and its eigenvalues
    S <- cov(x)
    S.eigenvalues <- eigen(S)$values
  })
  # Averaging eigenvalues repeated 1000 times
  mat[i,] <- apply(res,1,mean)
}
# y limit for plot
yl <- max(rbind(mat,e.tre))

pdf(file = "Ch01_Fig01.pdf") 
# plotting true and sample eigenvalues
plot(e.tre,type="b",lwd=3,pch=16, ylim = c(0,yl), xlab="Order", ylab = "Eigenvalues")
points(mat[1,], type="b", lwd=3,pch=16 ,col="red")
points(mat[2,], type="b", lwd=3,pch=16, col="blue")
points(mat[3,], type="b", lwd=3, pch=16,col="5")
points(mat[4,], type="b", lwd=3, pch=16,col="green")

legend("topright",inset = 0.03, legend = c("True","n= 25","n= 50","n= 100","n= 1000"),cex = 1,box.lwd = 2,box.lty = 2,text.font = 4, bg="gray",col = c("black","red","blue","5","green"), lty = 1,  pch = 16)

dev.off()