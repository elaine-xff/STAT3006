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
pgam = (1-pgamma(5, shape = 1/2, rate = 1)) * gamma(1/2)
M = 5^0.5 * exp(-5) / pgam # constant M
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
#hist(x_vec, breaks=10, freq=FALSE)
acceptance_prob = cnt / m


## Question 3: Integral Estimation

# 1. Use 5000 samples from Q2
sum_cos = sum(cos(x_vec))
integral = pgam/length(x_vec) * sum_cos

# 2. Use importance sampling
m = 5000
x = rexp(m, 1) + 5 # generate samples from g(x)
integral = exp(-5)* sum(cos(x) * x^(-0.5)) /m


## Question 4: Stratified Sampling

# Read Salary Data
SalaryData = read.delim("/Users/elainexfff_/Documents/STAT3006/Assignment 2/Coding_assignment_2/salary_data.txt", 
                        header = TRUE, sep = " ")

# (1) Randomly draw 100 samples from SalaryData and compute each sub-population's s.d.
rdm_salary = SalaryData[sample(nrow(SalaryData), size = 100),]
x1 = subset(rdm_salary, Age_Indicator == 1, select = c("Salary"))
sd_x1 = sd(as.numeric(unlist(x1)))
x2 = subset(rdm_salary, Age_Indicator == 2, select = c("Salary"))
sd_x2 = sd(as.numeric(unlist(x2)))
x3 = subset(rdm_salary, Age_Indicator == 3, select = c("Salary"))
sd_x3 = sd(as.numeric(unlist(x3)))

# (2) haha 











