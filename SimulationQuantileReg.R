Mu <- 5.00
N <- 100
epsilon <- rnorm(N,0,1)
Y_t <- Mu + epsilon

hist(Y_t)
plot(Y_t,type='p')

postBeta <- function(Beta,Y,N,p){
	post <- N*log(p) + N*log(1-p) - sum((abs(Y-Beta)+((2*p-1)*(Y-Beta)))/2)
	return(post)
}

updateBeta <- function(Beta,Y,N,p,Sm,Sigma,t,BetaMean){
	accept <- NULL
	count <- 0
	c <- .8
	ggamma <- 1/(t^c)
	
	proposed <- Beta + sqrt(Sm*Sigma)*rnorm(1,0,1)
	probab <- min(1,exp(postBeta(proposed,Y,N,p)-postBeta(Beta,Y,N,p)))

	if(runif(1)<probab){
		accept <- proposed
		count <- 1
	} else {
		accept <- Beta
	}

	lsm <- log(Sm)+ggamma*(probab - 0.234)
	Sm <- exp(lsm)
	Sigma <- Sigma + ggamma*(((Beta - BetaMean)^2)-Sigma)
	BetaMean <- BetaMean + ggamma*(Beta-BetaMean)
	return(list(accept,count,Sm,Sigma,BetaMean))
}

##########################################################################
Niter <- 20000
N <- 100
Beta.out1 <- array(NA, dim = Niter)
Count1 <- array(0, dim = Niter)
Sm1 <- array(NA, dim = Niter)
Sigma1 <- array(NA, dim = Niter)
BetaMean1 <- array(NA, dim = Niter)

Beta.out1[1] <- .5
Sm1[1] <- 2.4^2
Sigma1[1] <- 2
BetaMean1[1] <- 3

Beta.out2 <- array(NA, dim = Niter)
Count2 <- array(0, dim = Niter)
Sm2 <- array(NA, dim = Niter)
Sigma2 <- array(NA, dim = Niter)
BetaMean2 <- array(NA, dim = Niter)

Beta.out2[1] <- .5
Sm2[1] <- 2.4^2
Sigma2[1] <- 2
BetaMean2[1] <- 3

Beta.out3 <- array(NA, dim = Niter)
Count3 <- array(0, dim = Niter)
Sm3 <- array(NA, dim = Niter)
Sigma3 <- array(NA, dim = Niter)
BetaMean3 <- array(NA, dim = Niter)

Beta.out3[1] <- .5
Sm3[1] <- 2.4^2
Sigma3[1] <- 2
BetaMean3[1] <- 3

Beta.out4 <- array(NA, dim = Niter)
Count4 <- array(0, dim = Niter)
Sm4 <- array(NA, dim = Niter)
Sigma4 <- array(NA, dim = Niter)
BetaMean4 <- array(NA, dim = Niter)

Beta.out4[1] <- .5
Sm4[1] <- 2.4^2
Sigma4[1] <- 2
BetaMean4[1] <- 3

t1 <- Sys.time()
for(i in 2:Niter){
	Beta1 <- updateBeta(Beta.out1[i-1],Y_t,N,.05,Sm1[i-1],Sigma1[i-1],i,BetaMean1[i-1])
	Beta.out1[i] <- Beta1[[1]]
	Count1[i] <- Beta1[[2]]
	Sm1[i] <- Beta1[[3]]
	Sigma1[i] <- Beta1[[4]]
	BetaMean1[i] <- Beta1[[5]]

	Beta2 <- updateBeta(Beta.out2[i-1],Y_t,N,.25,Sm2[i-1],Sigma2[i-1],i,BetaMean2[i-1])
	Beta.out2[i] <- Beta2[[1]]
	Count2[i] <- Beta2[[2]]
	Sm2[i] <- Beta2[[3]]
	Sigma2[i] <- Beta2[[4]]
	BetaMean2[i] <- Beta2[[5]]

	Beta3 <- updateBeta(Beta.out3[i-1],Y_t,N,.75,Sm3[i-1],Sigma3[i-1],i,BetaMean3[i-1])
	Beta.out3[i] <- Beta3[[1]]
	Count3[i] <- Beta3[[2]]
	Sm3[i] <- Beta3[[3]]
	Sigma3[i] <- Beta3[[4]]
	BetaMean3[i] <- Beta3[[5]]

	Beta4 <- updateBeta(Beta.out4[i-1],Y_t,N,.95,Sm4[i-1],Sigma4[i-1],i,BetaMean4[i-1])
	Beta.out4[i] <- Beta4[[1]]
	Count4[i] <- Beta4[[2]]
	Sm4[i] <- Beta4[[3]]
	Sigma4[i] <- Beta4[[4]]
	BetaMean4[i] <- Beta4[[5]]
	
	print(i)
}
t2 <- Sys.time()
(t2-t1)

plot(Beta.out1[1001:20000],type='l')
plot(Beta.out2[1001:20000],type='l')
plot(Beta.out3[1001:20000],type='l')
plot(Beta.out4[1001:20000],type='l')

hist(Beta.out1[1001:20000])
hist(Beta.out2[1001:20000])
hist(Beta.out3[1001:20000])
hist(Beta.out4[1001:20000])

Beta1 <- mean(Beta.out1[1001:20000])
Beta2 <- mean(Beta.out2[1001:20000])
Beta3 <- mean(Beta.out3[1001:20000])
Beta4 <- mean(Beta.out4[1001:20000])

plot(Y_t,type='p')
abline(h=Beta1,col='red')
abline(h=Beta2,col='red')
abline(h=Beta3,col='red')
abline(h=Beta4,col='red')


meansVec <- c(Beta1,Beta2,Beta3,Beta4)
quants <- quantile(Y_t,probs=c(.05,.25,.75,.95))
cbind(meansVec,quants)
