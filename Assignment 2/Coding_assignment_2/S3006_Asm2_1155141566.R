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

hist(x_vec, breaks = seq(-0.5, 15.5, by = 0.5), freq = TRUE) # histogram


## Question 2: Accept-Reject method for truncated Gamma Distribution

# Samples from truncated Gamma Distribution
m = 5000 # sample size
pgam = (1-pgamma(5, shape = 1/2, rate = 1)) * gamma(1/2)
M = 5^(-0.5) * exp(-5) / pgam # constant M
cnt = 0 # acceptance counter
ratio = function(y){y^(-0.5) / 5^(-0.5)}
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
set.seed(3006)
rdm_salary = SalaryData[sample(nrow(SalaryData), size = 100),]
sd=rep(NA,3)
for(i in 1:3){
  sd[i]=sd(rdm_salary[rdm_salary$Age_Indicator==i,1])
}

# (2) stratified sampling -- sample number for each strata

# the proportion mu
mu1 = 1500 / 11000
mu2 = 4500 / 11000
mu3 = 5000 / 11000
mu = c(mu1, mu2, mu3)

n1 = 1000 * mu1 * sd[1] / sum(mu * sd)
n2 = 1000 * mu2 * sd[2] / sum(mu * sd)
n3 = 1000 * mu3 * sd[3] / sum(mu * sd)
(n=c(n1, n2, n3))

# (3) approximate the mean salary
# the sample number for each subpopulation
n_1 = 68
n_2 = 307
n_3 = 625
salary_1 = subset(SalaryData, Age_Indicator == 1)
salary_2 = subset(SalaryData, Age_Indicator == 2)
salary_3 = subset(SalaryData, Age_Indicator == 3)

set.seed(3006)
rdm_salary_1 = salary_1[sample(nrow(salary_1), size = n_1),]
rdm_salary_2 = salary_2[sample(nrow(salary_2), size = n_2),]
rdm_salary_3 = salary_3[sample(nrow(salary_3), size = n_3),]
stratified_sample = rbind(rdm_salary_1, rdm_salary_2, rdm_salary_3)
# get mean of the salary
mean(stratified_sample$Salary)







