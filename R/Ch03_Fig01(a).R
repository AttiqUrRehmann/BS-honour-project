# Distribution of γ values for diﬀerent samples of size n = {10,30,300}
# from a multivariate normal distribution with p = {10,30,100} for identity matrix as a true covariance
# matrix and as target

library(MASS)

# Various "n" and "p"
n <- c(10, 30, 300)
p <- c(10, 30, 100)

# Null matrix (gamma matrix) needed for gamma values inside for loops to be used for futher analysis
gamma <- matrix(NA, length(n), length(p))
res<- replicate(1000,{
  for(i in 1:length(p)){
    for(j in 1:length(n)){
      # Identity matrix as a true covariance matrix
      sigma <- diag(p[i])
      x <- mvrnorm(n=n[j], mu= rep(0,p[i]), Sigma=sigma)
      
      # Calculate gamma values and store them in null matrix
      gamma[i,j] <- 1/(1+sum(abs(cor(x) - diag(p[i])))/(p[i]))
    }
  }
  gamma
})

# Make the resulting gamma matrix as a vector repeated 1000 times 
vec <- as.vector(res)
# Boxes positions in the boxplot
arr <- as.vector(array(c(1,4,7,2,5,8,3,6,9),dim = c(3,3,1000)))

pdf(file="Ch03_Fig01(a).pdf")
# Margin
par(mar=c(7,5,1.5,1))
boxplot(vec~arr, outline=FALSE, ylab= expression(gamma), ylim=c(0,1), las=1, xaxt="n",at=c(1,3,5,8,10,12,15,17,19),cex.axis=1.5, cex.lab=1.5)
axis(1, at=c(1,3,5,8,10,12,15,17,19), labels = rep(c(10, 30, 300),3),cex.axis=1.5)
mtext("(a)",at=0.5, line = 0.2, cex= 1.5)
mtext("n =",at=-0.4, line = -28.6, cex= 1.5)
mtext("--",at= 1:5, line = -30, cex= 1.5)
mtext("--",at= 8:12, line = -30, cex= 1.5)
mtext("--",at= 15:19, line = -30, cex= 1.5)
mtext("p =",at= -0.4, line = -31, cex= 1.5)
mtext("10",at= 3, line = -31, cex= 1.5)
mtext("30",at= 10, line = -31, cex= 1.5)
mtext("100",at= 17, line = -31, cex= 1.5)

dev.off()