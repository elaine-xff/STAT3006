## Question 1: Inverse method for Poisson Distribution

# Samples from Poisson Distribution
m = 5000 # sample size
lam = 5 # lambda in Poisson
x = 0:20
prob_cumsum = ppois(x, lambda=lam) # cumulative sum of probability masses
u_vec = runif(m) # uniform samples
x_vec = NULL # Poisson samples

for (i in 1:m) {
  x_vec = c(x_vec, min(which(prob_cumsum>=u_vec[i]))-1)
}

hist(x_vec, breaks = seq(-0.5, 20.5, by = 0.5), freq = FALSE) # histogram


## Question 2: Accept-Reject method for truncated Gamma Distribution

# Samples from truncated Gamma Distribution
m = 5000 # sample size
M = 5^0.5 * exp(-5) / ((1-pgamma(5, shape = 1/2, rate = 1)) * gamma(1/2)) # constant M
cnt = 0 # acceptance counter
ratio = function(y){y^(-0.5) / 5^0.5}
x_vec = NULL
for (i in 1:m){
  y = rexp(1, 1) + 5
  u = runif(1)
  if(u <= ratio(y)){
    x_vec = c(x_vec, y)
    cnt = cnt + 1
  }
}
hist(x_vec, breaks=10, freq=FALSE)
acceptance_prob = cnt / m

