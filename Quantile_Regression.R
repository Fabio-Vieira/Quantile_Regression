###############################################################################
#############################Quantile Regression###############################

#This a code to compare quantile regressin to normal linear regression in the
#case the relationship between x a y violates the constant variance assumption

###Simulation
N <- 100
p <- 0.5
tau <- 2/(p * (1 - p)); theta <- (1 - 2*p)/(p * (1 - p))
z <- rexp(N,1)
x <- cbind(rep(1,N), rnorm(N,0,1))
beta <- c(3,1)
y <- x %*% beta + theta*z + sqrt(tau*z)*rnorm(N)

hist(y)

##############################################################################
library(MASS)
library(GIGrvg)

updateBetaQuant <- function(B0, X, Y, Z, p){
    tau <- 2/(p * (1 - p)); theta <- (1 - 2*p)/(p * (1 - p))
    Sig <- diag(tau*Z)
    
    Inv.B <- solve((t(X)%*%Sig%*%X) + B0)
    Beta.T <- Inv.B %*% (t(X)%*%Sig)%*%(Y - theta * Z)
    return(mvrnorm(1,Beta.T,Inv.B))
}

updateZ <- function(Y, X, Beta, p){
    tau <- 2/(p * (1 - p)); theta <- (1 - 2*p)/(p * (1 - p))
  
    lambda <- 0.5
    delta <- sqrt(((Y - X%*%Beta)^2)/tau)
    gamma <- sqrt(2 + (theta^2/(tau)))
  
  v <- NULL
  for(i in 1:length(Y)){v <- c(v,rgig(1,lambda,delta[i],gamma))}
  return(v)
}

##############################################################################
Niter <- 10000
N <- 100
X <- x
Y <- y
p <- 0.5
Beta.out <- array(NA, dim = c(Niter,dim(X)[2]))
Z.out <- array(NA, dim = c(Niter,N))


#####Chain starting points
Beta.out[1,] <- beta
Z.out[1,] <- z

######Prior
B0 <- diag(100,2,2)

for(i in 2:Niter){
    Beta.out[i,] <- updateBeta(B0, X, Y, Z.out[i-1,], p)
    Z.out[i,] <- updateZ(Y, X, Beta.out[i,], p)
    print(i)
}

#####Checking convergence
Nburn <- 1000

plot(Beta.out[-c(1:Nburn),1],type='l')
abline(h=beta[1],col='red')
mean(Beta.out[-c(1:Nburn),1])
hist(Beta.out[-c(1:Nburn),1])
abline(v=beta[1],col='red')
abline(v=quantile(Beta.out[-c(1:Nburn),1],prob=0.025),lty=2,col='blue')
abline(v=quantile(Beta.out[-c(1:Nburn),1],prob=0.975),lty=2,col='blue')

plot(Beta.out[-c(1:Nburn),2],type='l')
abline(h=beta[2],col='red')
mean(Beta.out[-c(1:Nburn),2])
hist(Beta.out[-c(1:Nburn),2])
abline(v=beta[2],col='red')
abline(v=quantile(Beta.out[-c(1:Nburn),2],prob=0.025),lty=2,col='blue')
abline(v=quantile(Beta.out[-c(1:Nburn),2],prob=0.975),lty=2,col='blue')

plot(z,type='l')
lines(colMeans(Z.out[-c(1:Nburn),]),type='l',col='red')

###############################################################################
#Auto-MPG data set from UCI repository
data <- as.data.frame(read.csv("auto-mpgdata.csv", sep = ";"))

#Verifying the relationship between acceleration and the fuel consumption in 
#miles per gallon.
plot(data[,6],data[,1]) #Variance is clearly increasing over time, which
                        # violates a basic assumption of linear regression.

##############################################################################
#Fitting a normal linear regression

updateBeta <- function(Sig2Be,X,Y,Sig2){ 
  #This function provides a Gibbs sampler for the regression parameters, 
  # and it can adapt to simple and multiple linear regression, all that is
  #need is to adjust the dimension of the design matrix
  n <- dim(X)[2]
  B0 <- diag(Sig2Be,n,n)
  Sigma <- Sig2 * diag(1,ncol=n,nrow=n)
  sig2.post <- solve((t(X)%*%X)%*%Sigma + B0)
  m.post <- (t((t(X)%*%Y))%*%Sigma) %*% sig2.post
  
  return(mvrnorm(1,m.post,sig2.post))
}

updateSig2 <- function(a,b,N,Y,X,Beta){
  #This function provides a Gibbs sampler for the contant variance in a normal
  #linear regression model
  a.post <- N/2 + a
  b.post <- b + (t(Y - X%*%Beta)%*%(Y - X%*%Beta))*0.5
  
  return(1/rgamma(1,a.post,b.post))
}

##############################################################################
####Variables
X <- cbind(rep(1,length(data[,1])), data[,6])
Y <- data[,1]
N <- length(Y)

##MCMC
Niter <- 30000
Beta.out <- array(NA, dim = c(Niter,dim(X)[2]))
Sig2.out <- array(NA,dim=Niter)

Beta.out[1,] <- rep(1,1)
Sig2.out[1] <- 1

for(i in 2:Niter){
  Beta.out[i,] <- updateBeta(0.001,X,Y,Sig2.out[i-1])
  Sig2.out[i] <- updateSig2(0.01,0.01,N,Y,X,Beta.out[i,])
  print(i)
}

#####Briefly checking convergence of parameters 
Nburn <- 10000 #burn-in
plot(Beta.out[-c(1:Nburn),1],type='l')
mean(Beta.out[-c(1:Nburn),1])
hist(Beta.out[-c(1:Nburn),1])

plot(Beta.out[-c(1:Nburn),2],type='l')
mean(Beta.out[-c(1:Nburn),2])
hist(Beta.out[-c(1:Nburn),2])

plot(Sig2.out[-c(1:Nburn)],type='l')
mean(Sig2.out[-c(1:Nburn)])
hist(Sig2.out[-c(1:Nburn)])

#Ploting the mean line on the scatter plot
plot(data[,6],data[,1])
means <- colMeans(Beta.out[-c(1:Nburn),])
line <- means[1] + means[2] * data[,6]
lines(data[,6],line,col='red')

##################################################################################
#####Calculating the Fitted Values

Beta <- Beta.out[-c(1:Nburn),]
Sig2 <- Sig2.out[-c(1:Nburn)]

####Fitting the values
Fitted <- array(NA, dim = c(nrow(Beta),N))

for(i in 1:nrow(Fitted)){
  Fitted[i,] <- X %*% Beta[i,] + sqrt(Sig2[i])*rnorm(N)
}
medians <- apply(Fitted,2,median)

####Calculating the residuals 
Residuals <- Y - medians
plot(Residuals)
####Testing for normality
shapiro.test(Residuals) #Residuals do not have a normal distribution

##############################################################################
#Let's see what happens when we fit a quantile regression model

Niter <- 30000
X <- cbind(rep(1,length(data[,1])), data[,6])
Y <- data[,1]
N <- length(Y)
p <- c(0.25, 0.5, 0.75, 0.9)
Beta.out1 <- array(NA, dim = c(Niter,dim(X)[2]))
Beta.out2 <- array(NA, dim = c(Niter,dim(X)[2]))
Beta.out3 <- array(NA, dim = c(Niter,dim(X)[2]))
#Beta.out4 <- array(NA, dim = c(Niter,dim(X)[2]))

Z.out1 <- array(NA, dim = c(Niter,N))
Z.out2 <- array(NA, dim = c(Niter,N))
Z.out3 <- array(NA, dim = c(Niter,N))
Z.out4 <- array(NA, dim = c(Niter,N))


#####Chain starting points
Beta.out1[1,] <- rep(1,dim(X)[2])
Z.out1[1,] <- rep(1,N)
Beta.out2[1,] <- rep(1,dim(X)[2])
Z.out2[1,] <- rep(1,N)
Beta.out3[1,] <- rep(1,dim(X)[2])
Z.out3[1,] <- rep(1,N)
Beta.out4[1,] <- rep(1,dim(X)[2])
Z.out4[1,] <- rep(1,N)

######Prior
B0 <- diag(100,2,2)

for(i in 2:Niter){
  Beta.out1[i,] <- updateBetaQuant(B0, X, Y, Z.out1[i-1,], p[1])
  Z.out1[i,] <- updateZ(Y, X, Beta.out1[i,], p[1])
  
  Beta.out2[i,] <- updateBetaQuant(B0, X, Y, Z.out2[i-1,], p[2])
  Z.out2[i,] <- updateZ(Y, X, Beta.out2[i,], p[2])

  Beta.out3[i,] <- updateBetaQuant(B0, X, Y, Z.out3[i-1,], p[3])
  Z.out3[i,] <- updateZ(Y, X, Beta.out3[i,], p[3])  
  
  Beta.out4[i,] <- updateBetaQuant(B0, X, Y, Z.out4[i-1,], p[4])
  Z.out4[i,] <- updateZ(Y, X, Beta.out4[i,], p[4])  
  
  print(i)
}

#Checking convergence

plot(Beta.out1[-c(1:Nburn),1],type='l')
mean(Beta.out1[-c(1:Nburn),1])
hist(Beta.out1[-c(1:Nburn),1])

plot(Beta.out1[-c(1:Nburn),2],type='l')
mean(Beta.out1[-c(1:Nburn),2])
hist(Beta.out1[-c(1:Nburn),2])

plot(Beta.out2[-c(1:Nburn),1],type='l')
mean(Beta.out2[-c(1:Nburn),1])
hist(Beta.out2[-c(1:Nburn),1])

plot(Beta.out2[-c(1:Nburn),2],type='l')
mean(Beta.out2[-c(1:Nburn),2])
hist(Beta.out2[-c(1:Nburn),2])

plot(Beta.out3[-c(1:Nburn),1],type='l')
mean(Beta.out3[-c(1:Nburn),1])
hist(Beta.out3[-c(1:Nburn),1])

plot(Beta.out3[-c(1:Nburn),2],type='l')
mean(Beta.out3[-c(1:Nburn),2])
hist(Beta.out3[-c(1:Nburn),2])

#Ploting quantile lines on the scatter plot
plot(data[,6],data[,1],ylim=c(9,50))
means1 <- colMeans(Beta.out1[-c(1:Nburn),])
means2 <- colMeans(Beta.out2[-c(1:Nburn),])
means3 <- colMeans(Beta.out3[-c(1:Nburn),])
means4 <- colMeans(Beta.out4[-c(1:Nburn),])
line1 <- means1[1] + means1[2] * data[,6]
line2 <- means2[1] + means2[2] * data[,6]
line3 <- means3[1] + means3[2] * data[,6]
line4 <- means4[1] + means4[2] * data[,6]

lines(data[,6],line1,col='red')
lines(data[,6],line2,col='red')
lines(data[,6],line3,col='red')
lines(data[,6],line4,col='red')

#############################################################################
#This plot shows the slopes of the mean model compared to the quantiles
plot(1, type="n", xlab="", ylab="", xlim=c(1,4), ylim=c(1,1.5),
     main = "Comparing Slopes")
abline(h=mean(Beta[,2]),col='red')
means <- c(means1[2],means2[2],means3[2],means4[2])
lines(1:4,means,type='o')


