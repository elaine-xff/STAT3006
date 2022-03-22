################################################
# Example 1: bootstrap 
################################################
rm(list = ls())
weight_data <- c(57, 60, 52, 49, 56, 46, 51, 63, 49, 57, 59, 54, 56, 59, 57, 52,61, 59, 53, 59,51,51, 56, 58,46,53)
N <- 26
hist(weight_data, breaks = 15) #investigate the distribution's property.
T1 <- mean(weight_data)
T2 <- var(weight_data)

B <- 1000 #bootstrap number
T1_star <- rep(NA, B)
T2_star <- rep(NA, B)
for(b in 1:B){
	data_resampled <- sample(weight_data, N, replace = TRUE)
	T1_star[b] <- mean(data_resampled)
	T2_star[b] <- var(data_resampled)
}

#the variance of T1
var(T1_star)

#the bias of T1
mean(T1_star) - T1

#the variance of T2
var(T2_star)

#the bias of T2
mean(T2_star) - T2

################################################
# Example 2: bootstrap check its performance
################################################
rm(list = ls())
set.seed(12345)
N <- 30
original_data <- rnorm(N)
hist(original_data, breaks = 10)

T1 <- mean(original_data)
T2 <- var(original_data)

B <- 1000 #bootstrap number
T1_star <- rep(NA, B)
T2_star <- rep(NA, B)
for(b in 1:B){
	data_resampled <- sample(original_data, N, replace = TRUE)
	T1_star[b] <- mean(data_resampled)
	T2_star[b] <- var(data_resampled)
}

#the variance of T1
var(T1_star)
#truth
1/N

#the bias of T1
mean(T1_star) - T1
#truth
0

#the variance of T2
var(T2_star)
#truth
2 / (N - 1)

#the bias of T2
mean(T2_star) - T2
#truth
0

################################################
# Example 3: permutation test
################################################
rm(list = ls())
set.seed(12345)
tA <- c(258,171, 193, 199, 230, 243 ,248, 250, 267 ,316 ,327, 329)
tB <- c(141, 148, 169, 181, 143, 213, 229, 234, 257, 220, 271, 239)
Tobs <- mean(tA) - mean(tB)
N_A <- length(tA)
N_B <- length(tB)
N <- length(tA) + length(tB)
t_AB <- c(tA, tB)
B <- 1000
T_star <- rep(NA, B)
for(b in 1:B){
	ind <- sample(1:N, N)
	shuffled <- t_AB[ind]
	T_star[b] <- mean(shuffled[1:N_A]) - mean(shuffled[(N_A+1):N])	
}

1 / B * sum(abs(T_star) >= abs(Tobs))


