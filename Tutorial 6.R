######################################################
#beta-binomial distribution
######################################################
a <- 2
b <- 4
n <- 16
num_iter <- 5000
x <- rep(NA, num_iter)
theta <- rep(NA, num_iter)

#initial values
x[1] = 0
theta[1] =0.5

for(i in 2:num_iter){
	x[i] <- rbinom(1, n, theta[i-1])
	theta[i] <- rbeta(1, a + x[i], b + n - x[i])
}

#trace plots
par(mfrow = c(2,1))
plot(x, type = "l")
plot(theta, type = "l")

B <- 50
hist(x[(B+1):num_iter], breaks = seq(-0.5, 16.5, by = 1), freq = FALSE)

######################################################
#beta-binomial-poisson distribution
######################################################
rm(list = ls())
a <- 2
b <- 4
lambda <- 15

num_iter <- 5000
x <- rep(NA, num_iter)
theta <- rep(NA, num_iter)
n <- rep(NA, num_iter)

#initial values
x[1] = 0
theta[1] = 0.5
n[1] = 10

for(i in 2:num_iter){
	x[i] <- rbinom(1, n[i-1], theta[i-1])
	theta[i] <- rbeta(1, a + x[i], b + n[i-1] - x[i])
	n[i] <- x[i] + rpois(1, lambda*(1-theta[i]))
}

#trace plots
par(mfrow = c(3,1))
plot(x, type = "l")
plot(theta, type = "l")
plot(n, type = "l")


B <- 50
hist(x[(B+1):num_iter], breaks = seq(-0.5, 25.5, by = 1), freq = FALSE)


######################################################
#Bivariate normal (Gibbs sampler)
######################################################
rm(list = ls())
mu1 <- 3
mu2 <- 2
sigma1 <- 2
sigma2 <- 3
sigma12 <- 5

num_iter <- 5000
x1 <- rep(NA, num_iter)
x2 <- rep(NA, num_iter)
#initial values
x1[1] <- 0
x2[1] <- 0

for(i in 2:num_iter){
	x1[i] <- rnorm(1, mean = mu1 + sigma12 / sigma2^2 *(x2[i-1] - mu2), sd = 
					sqrt(sigma1^2 - sigma12^2 / sigma2^2))
	x2[i] <- rnorm(1, mean = mu2 + sigma12 / sigma1^2 *(x1[i] - mu1), sd = 
					sqrt(sigma2^2 - sigma12^2 / sigma1^2))
}

par(mfrow = c(2,1))
plot(x1, type = "l")
plot(x2, type = "l")

B <- 500
plot(x1[(B+1):num_iter], x2[(B+1):num_iter], cex = 0.6)
hist(x1[(B+1):num_iter], breaks = 20, freq = FALSE)
hist(x2[(B+1):num_iter], breaks = 20, freq = FALSE)

######################################################
#Bivariate normal (MH algorithm)
######################################################
rm(list = ls())
mu1 <- 3
mu2 <- 2
sigma1 <- 2
sigma2 <- 3
sigma12 <- 5
sigma_matr <- matrix(c(sigma1^2, sigma12, sigma12, sigma2^2), 2, 2)
mu_vec <- c(mu1, mu2)

num_iter <- 5000
x1 <- rep(NA, num_iter)
x2 <- rep(NA, num_iter)
#initial values
x1[1] <- 0
x2[1] <- 0

target_f <- function(x1, x2){
			x_vec <- c(x1, x2)
			temp <- det(sigma_matr)^(-1/2)*exp(-1/2* t(x_vec - mu_vec) %*% 									solve(sigma_matr) %*% (x_vec - mu_vec) )
			return(temp)	
}

accept <- 0

for(i in 2:num_iter){
	#proposal 
	y1 <- x1[i-1] + rnorm(1, sd = 1)
	y2 <- x2[i-1] + rnorm(1, sd = 1)
	
	#acceptance prob
	r <- min(target_f(y1,y2) / target_f(x1[i-1],x2[i-1]), 1)
	if(runif(1) < r){
		x1[i] <- y1
		x2[i] <- y2
		accept <- accept + 1
	}else{
		x1[i] <- x1[i-1]
		x2[i] <- x2[i-1]
	}
	
}

accept / num_iter

par(mfrow = c(2,1))
plot(x1, type = "l")
plot(x2, type = "l")

B <- 500
plot(x1[(B+1):num_iter], x2[(B+1):num_iter], cex = 0.6)
hist(x1[(B+1):num_iter], breaks = 20, freq = FALSE)
hist(x2[(B+1):num_iter], breaks = 20, freq = FALSE)




