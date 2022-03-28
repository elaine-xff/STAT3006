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













