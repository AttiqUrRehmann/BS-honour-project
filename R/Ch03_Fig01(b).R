# Distribution of γ values for diﬀerent samples of size n = {10,30,300}
# from a multivariate normal distribution with p = {10,30,100} simulated
# 1000 times for AR(1) structure as a true covariance matrix and as a target


library(MASS) 

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

n <- c(10, 30, 300)
p <- c(10, 30, 100)
t <- 0.5 # True t value for AR(1) covariance structure as a true covariance matrix

gamma <- matrix(NA, 3,3)
res <- replicate(1000,{
  
  for(i in 1:length(p)){
    for(j in 1:length(n)){
      sigma <-  t ^ outer(1:p[i], 1:p[i], function(aa, bb) abs(aa - bb))
      x <- mvrnorm(n=n[j], mu= rep(0,p[i]), Sigma=sigma)
      # Estimate t using t.ar1() function
      t.hat <- t.ar1(x)
      # Estimate true covariance matrix 
      AR1.hat <- t.hat ^ outer(1:p[i], 1:p[i], function(aa, bb) abs(aa - bb)) 
      # Calculate gamma values see section 3.3
      k <- seq(0,p[i]-1)
      gamma[i,j] <- 1/(1+sum(abs(cor(x) - AR1.hat))/(p[i]+sum(2*k*(t.hat)^(p[i]-k))))
    }
  }
  gamma
})
vec <- as.vector(res)
arr <- as.vector(array(c(1,4,7,2,5,8,3,6,9),dim = c(3,3,1000)))

pdf(file="Ch03_Fig01(b).pdf")
par(mar=c(7,5,1.5,1))
boxplot(vec~arr, outline=FALSE, ylab= expression(gamma),
        ylim=c(0,1), las=1, xaxt="n", at=c(1,3,5,8,10,12,15,17,19),
        cex.axis=1.5, cex.lab=1.5)
axis(1, at=c(1,3,5,8,10,12,15,17,19), labels = rep(c(10, 30, 300),3),
     cex.axis=1.5)
mtext("(b)",at=0.5, line = 0.2, cex= 1.5)
mtext("n =",at=-0.4, line = -28.6, cex= 1.5)
mtext("--",at= 1:5, line = -30, cex= 1.5)
mtext("--",at= 8:12, line = -30, cex= 1.5)
mtext("--",at= 15:19, line = -30, cex= 1.5)
mtext("p =",at= -0.4, line = -31, cex= 1.5)
mtext("10",at= 3, line = -31, cex= 1.5)
mtext("30",at= 10, line = -31, cex= 1.5)
mtext("100",at= 17, line = -31, cex= 1.5)
dev.off()