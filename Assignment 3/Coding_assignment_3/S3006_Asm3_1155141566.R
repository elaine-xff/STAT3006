## Question 1: Hybrid Gibbs Sampler to estimate Poisson Distribution \lambda

rm(list = ls())

N = 2000 # iteration number
lambda_t = rep(NA, N) # the parameter to be inferred
Y_t = matrix(NA, N, 78) # 78 unobserved variables (y_i >= 5)

# Initialization
lambda_t[1] = 1.5
Y_t[1, ] = rep(5, 78)

# Hybrid Gibbs Sampler
accept_num = 0
set.seed(3006)
for (i in 2:N) {
  # sample lambda
  lambda_t[i] = rgamma(1, shape = 422 + sum(Y_t[i-1, ]), rate = 500)
  
  # sample 78 unobserved variables (MH-Step)
  for (j in 1:78) {
    # the proposal distribution of y's
    y_proposal = Y_t[i-1, j] + sample(c(-1, 0, 1), 1)
    # set the accept-reject ratio
    r = min( (lambda_t[i])^(y_proposal - Y_t[i-1, j]) * factorial(Y_t[i-1, j]) / factorial(y_proposal), 1)
    if((runif(1) < r) & (y_proposal >= 5)){
      Y_t[i, j] = y_proposal # accept
      accept_num = accept_num + 1
    }
    else{
      Y_t[i, j] = Y_t[i-1, j] # reject
    }
  }
  
}

accept_ratio = accept_num / (78*N) 

B = 1000 # first 1000 iterations are burn-in
# the estimation for lambda is the posterior mean
estimated_lambda = mean(lambda_t[(B+1):N])


## Question 2: Gibbs Sampler for clustering

rm(list = ls())

# Read Data
Q2Data = read.delim("/Users/elainexfff_/Documents/STAT3006/Assignment 3/Coding_assignment_3/Assg3_Q2.txt", 
                        header = TRUE, sep = " ")

#function used to sample Dirichlet distributed r.v.
rDirichlet <- function(alpha_vec){
  num <- length(alpha_vec)
  temp <- NULL
  for(i in 1:num){
    temp <- c(temp, rgamma(1, shape = alpha_vec[i], rate = 1))
  }
  return(temp/sum(temp))
} 


N = 5000 # iteration number
pi_t = array(NA, dim = c(1000, 3, N) ) # to store pi_1, pi_2, and pi_3 for each iteration
theta_t = array(NA, dim = c(3, 3, N) ) # to store theta for each iteration
z_t = matrix(NA, N, 1000) # to store Z for each iteration
alpha_t = matrix(NA, N, 3) # alpha values for Dirichlet Distribution

set.seed(3006)

for (i in 1:N) {
  if(i==1){
    # initialization
    alpha_t[1,] = rep(2, 1*3)
    for (j in 1:1000) {
      pi_t[j, , 1] = rDirichlet(alpha_t[1,])
    }
    
    theta_t[ , , 1] = matrix(1/9, 3, 3)
    
    binom_k_1 = pi_t[, 1, 1] * dbinom(Q2Data[,1], 10 * 1, theta_t[1,1,1]) * dbinom(Q2Data[,2], 10 * 2, theta_t[2,1,1]) * dbinom(Q2Data[,3], 10 * 3, theta_t[3,1,1])
    binom_k_2 = pi_t[, 2, 1] * dbinom(Q2Data[,1], 10 * 1, theta_t[1,2,1]) * dbinom(Q2Data[,2], 10 * 2, theta_t[2,2,1]) * dbinom(Q2Data[,3], 10 * 3, theta_t[3,2,1])
    binom_k_3 = pi_t[, 3, 1] * dbinom(Q2Data[,1], 10 * 1, theta_t[1,3,1]) * dbinom(Q2Data[,2], 10 * 2, theta_t[2,3,1]) * dbinom(Q2Data[,3], 10 * 3, theta_t[3,3,1])
    
    prob_1 = binom_k_1 / (binom_k_1 + binom_k_2 + binom_k_3)
    prob_2 = binom_k_2 / (binom_k_1 + binom_k_2 + binom_k_3)
    prob_3 = binom_k_3 / (binom_k_1 + binom_k_2 + binom_k_3)
    for (j in 1:1000) {
      z_t[1, j] = sample(1:3, 1, prob = c(prob_1[j], prob_2[j], prob_3[j]), replace = TRUE)
    }
    
  }
  else{
    alpha_t[i, ] = alpha_t[i-1, ] + c(sum(z_t[i-1, ] == 1), sum(z_t[i-1, ] == 2), sum(z_t[i-1, ] == 3))
    for (j in 1:1000) {
      pi_t[j, , i] = rDirichlet(alpha_t[i, ])
    }
    
    for (j in 1:3) {
      for (k in 1:3) {
        sum_x = 0
        for (z in 1:1000) {
          if(z_t[i-1, z] == k){
            sum_x = sum_x + Q2Data[z,j]
          }
        }
        theta_t[j, k, i] = rbeta(1, 1+sum_x, 1+10*j*sum(z_t[i-1, ] == k) - sum_x)
      }
    }
    
    binom_k_1 = pi_t[, 1, i] * dbinom(Q2Data[,1], 10 * 1, theta_t[1,1,i]) * dbinom(Q2Data[,2], 10 * 2, theta_t[2,1,i]) * dbinom(Q2Data[,3], 10 * 3, theta_t[3,1,i])
    binom_k_2 = pi_t[, 2, i] * dbinom(Q2Data[,1], 10 * 1, theta_t[1,2,i]) * dbinom(Q2Data[,2], 10 * 2, theta_t[2,2,i]) * dbinom(Q2Data[,3], 10 * 3, theta_t[3,2,i])
    binom_k_3 = pi_t[, 3, i] * dbinom(Q2Data[,1], 10 * 1, theta_t[1,3,i]) * dbinom(Q2Data[,2], 10 * 2, theta_t[2,3,i]) * dbinom(Q2Data[,3], 10 * 3, theta_t[3,3,i])
    
    prob_1 = binom_k_1 / (binom_k_1 + binom_k_2 + binom_k_3)
    prob_2 = binom_k_2 / (binom_k_1 + binom_k_2 + binom_k_3)
    prob_3 = binom_k_3 / (binom_k_1 + binom_k_2 + binom_k_3)
    for (j in 1:1000) {
      z_t[i, j] = sample(1:3, 1, prob = c(prob_1[j], prob_2[j], prob_3[j]), replace = TRUE)
    }
    
  }

}

B = 3000 # burn-in period

collected_pi_t = pi_t[,,(B+1):N]
#estimated_pi = rowMeans(collected_pi_t, dims = 2)
estimated_pi = apply(collected_pi_t, c(1,2), mean)

estimated_theta = matrix(NA, 3, 3)
for (i in 1:3) {
  for (j in 1:3){
    estimated_theta[i, j] = mean(theta_t[i, j, (B+1):N])
  }
}

estimated_z = matrix(NA, 1000, 1)
for(j in 1:1000){
  temp <- table(z_t[(B+1):N, j])
  ind <- which.max(temp)
  estimated_z[j] = names(temp)[ind]
}


# plot
dev.off()

plot(theta_t[1,1,], ylim = c(0.75,0.85))

plot(z_t[, 1])

plot(pi_t[,,1])



## Question 3: Hybrid Gibbs Sampler

rm(list = ls())
set.seed(3006)

N = 10000 #iteration number
p_1_t = matrix(NA, N, 4) # prob. for stage 1
p_2_t = matrix(NA, N, 4) # prob. for stage 2
y_t = matrix(NA, N, 4) # col 1 for y11, 2 for y12, 3 for y22, 4 for y24

#function used to sample Dirichlet distributed r.v.
rDirichlet <- function(alpha_vec){
  num <- length(alpha_vec)
  temp <- NULL
  for(i in 1:num){
    temp <- c(temp, rgamma(1, shape = alpha_vec[i], rate = 1))
  }
  return(temp/sum(temp))
} 

# initialization
p_1_t[1, ] = rDirichlet(c(2,2,2,2))
p_2_t[1, ] = rDirichlet(c(2,2,2,2))
y_t[1, ] = c(100-22-31-20, 20, 20, 100-28-26-20)

# Hybrid Gibbs Sampler
accept_num = 0
for (i in 2:N) {
  # sample p's
  p_1_t[i, ] = rDirichlet(p_1_t[i-1, ] + c(y_t[i-1, 1], y_t[i-1, 2], 22, 31))
  p_2_t[i, ] = rDirichlet(p_2_t[i-1, ] + c(28, y_t[i-1, 3], 26, y_t[i-1, 4]))
  
  # sample unobserved y's
  # sample yi2 (i = 1,2)
  for (j in 2:3) {
    # get a proposal distribution depends on its last iteration
    if (y_t[i-1, j] == 15){
      y_proposal = y_t[i-1, j] + 1
    }
    if (y_t[i-1, j] == 32){
      y_proposal = y_t[i-1, j] - 1
    }
    if ( (y_t[i-1, j]>15) & (y_t[i-1, j]<32) ){
      y_proposal = y_t[i-1, j] + sample(c(-1, 1), 1)
    }
    
    if( j==2 ){
      p_2 = p_1_t[i, 2]
    }else{
      p_2 = p_2_t[i, 2]
    }
    
    r = min ( p_2^(y_proposal)/factorial(y_proposal)/p_2^(y_t[i-1, j])/factorial(y_t[i-1, j]), 1 )
    
    if(runif(1) < r){
      y_t[i, j] = y_proposal # accept proposal
      accept_num = accept_num + 1
    }else{
      y_t[i, j] = y_t[i-1, j] # reject proposal
    }
    
  }
  
  y_t[i, 1] = 100 - 22 - 31 - y_t[i, 2]
  y_t[i, 4] = 100 - 28 - 26 - y_t[i, 3]
  
}

accept_ratio = accept_num / (2*N) # to check if the acceptance rate is good enough


B = 5000 # burn-in period
estimated_p_1 = rep(NA, 4)
estimated_p_2 = rep(NA, 4)
for (j in 1:4) {
  estimated_p_1[j] = mean(p_1_t[(B+1):N, j])
  estimated_p_2[j] = mean(p_2_t[(B+1):N, j])
}

