##############################################################
#Example 1: Hybrid Gibbs Sampler 
##############################################################
rm(list = ls())
N <- 10000 #iteration number 
lambda_t <- rep(NA, N)   # the parameter we want to infer 
Y_t <- matrix(NA, N, 13) # 13 unobserved variable

#initial values
lambda_t[1] <- 1.5
Y_t[1, ] <- rep(5, 13)


#hybrid Gibbs sampler
accept_num <- 0
for(i in 2:N){
	#sample lambda
	lambda_t[i] <- rgamma(1, shape = 313 + sum(Y_t[i-1, ]), rate = 360)
	
	#sample 13 unobserved variables (MH step)
	for(j in 1:13){
		y_proposal <- Y_t[i-1,j] + sample(c(-1,0,1),1) 
									  # get a proposal, the proposal distribution
									  # depends on the last iteration
		r <- min( (lambda_t[i] )^(y_proposal - Y_t[i-1, j]) * 
					factorial(Y_t[i-1,j]) / factorial(y_proposal), 1)
		if((runif(1) < r) & (y_proposal >= 4)){
			Y_t[i, j] <- y_proposal #accept proposal
			accept_num <- accept_num + 1
		}else{
			Y_t[i, j] <- Y_t[i-1, j]#reject proposal
		}
		
	}
}

accept_num / (13*N) #acceptance rate is good
	
par(mfrow= c(2,1))

plot(lambda_t, type = "l")

plot(Y_t[, 1], type = "l")

#first 2000 as burn-in
B <- 2000
#the estimate for lambda
#posterior mean
mean(lambda_t[(B+1):N])

#estimates for nonobserved variables
#posterior modes 
for(j in 1:13){
	temp <- table(Y_t[(B+1):N, j])
	ind <- which.max(temp)
	print(names(temp)[ind])
}

##############################################################
#Example 2: Dirichlet Distribution
##############################################################
rm(list = ls())
N <- 10000 #iteration number
p_t <- matrix(NA, N, 3)

#function used to sample Dirichlet distributed r.v.
rDirichlet <- function(alpha_vec){
	num <- length(alpha_vec)
	temp <- NULL
	for(i in 1:num){
		temp <- c(temp, rgamma(1, shape = alpha_vec[i], rate = 1))
	}
	return(temp/sum(temp))
} 

##################################
#simulate data
p_truth <- c(1/6,1/2, 1/3)
temp <- sample(1:3, 1000, prob = p_truth, replace = TRUE)
Y1 <- sum(temp == 1)
Y2 <- sum(temp == 2)
Y3 <- sum(temp == 3)
##################################



p_t[1, ] <- c(1/3,1/3,1/3)

for(i in 2:N){
	#sample p_t
	p_t[i, ] <- rDirichlet(c(1/3 + Y1, 1/3 + Y2, 1/3 + Y3))
}


par(mfrow = c(3,1))

plot(p_t[ ,1], type = "l")
abline(h = p_truth[1], col = "red")
plot(p_t[ ,2], type = "l")
abline(h = p_truth[2], col = "red")
plot(p_t[ ,3], type = "l")
abline(h = p_truth[1], col = "red")


	
B <- 2000
mean(p_t[(B+1):N,1])
mean(p_t[(B+1):N,2])
