## Question 1: Inverse method for Poisson Distribution

# Samples from Poisson Distribution
m = 5000 # sample number
lam = 5 # lambda in Poisson
x = 0:20
prob_cumsum = ppois(x, lambda=lam) # cumulative sum of probability masses
u_vec = runif(m) # uniform samples
x_vec = NULL # Poisson samples

for (i in 1:m) {
  x_vec = c(x_vec, min(which(prob_cumsum>=u_vec[i]))-1)
}

hist(x_vec, breaks = seq(-0.5, 20.5, by = 1), freq = FALSE) # histogram

