### Question 1: BISECTION METHOD

# Define the function f(x)
f_func <- function(x) {
  return(x^3 + 6 * x^2 + pi * x -12)
}

# Define the bisection method function
# Argument: endpoints, error tolerance, and stopping criterion -- max iteration 
bisec <- function(lower, upper, tolerance, max_iteration){
  turn = 0
  while (turn < max_iteration){
    midpoint = (lower + upper)/2
    # Stop the loop if midpoint is found or the error is less than the tolerance
    if ((f_func(midpoint)==0) | ((upper - lower)/2 < tolerance))
      return(midpoint)
    turn = turn + 1
    # Bisection part
    if (f_func(midpoint) * f_func(lower) > 0)
      lower = midpoint
    else
      upper = midpoint
  }
}

# Determine the number of zeros of function f
a_0 = (-12 - sqrt(12 * (12 - pi)))/6
b_0 = (-12 + sqrt(12 * (12 - pi)))/6
f_func(a_0)
f_func(b_0)

# Initialize variables to find 3 zeros
a_1 = -5
b_1 = -4
a_2 = -3
b_2 = -2
a_3 = 1
b_3 = 2
tolerance = 0.00001
max_iteration = 1000

# Bisection method for the first, second, third zero respectively
bisec(a_1, b_1, tolerance, max_iteration)
bisec(a_2, b_2, tolerance, max_iteration)
bisec(a_3, b_3, tolerance, max_iteration)



### Question 2: POISSON REGRESSION--NEWTON'S METHOD

# Read Poisson Regression Data
PoisRegData = read.delim("/Users/elainexfff_/Documents/STAT3006/Assignment 1/Coding_assignment_1/PoisRegData.txt", 
                         header = TRUE, sep = " ")

# Newton's Iteration for Poisson Regression
# Argument: Initial guess alpha_0, beta_0, gamma_0 and error tolerance
tolerance = 0.00001
alpha_0 = 1
beta_0 = 1
gamma_0 = 1

# Create the 'Goal' vector (to be optimized)
partial_mat <- function(alpha, beta, gamma, reg_data){
  x = reg_data['x']
  y = reg_data['y']
  partial_alpha = sapply(y - exp(alpha + beta * x + gamma * x^2), sum)
  partial_beta = sapply(x*y - x * exp(alpha + beta * x + gamma * x^2), sum)
  partial_gamma = sapply(x^2 * y - x^2 * exp(alpha + beta * x + gamma * x^2), sum)
  vec_f = matrix(c(partial_alpha, partial_beta, partial_gamma), 
                 nrow=3, ncol=1, byrow = FALSE)
  return(vec_f)
}

#  Calculate the Jacobian matrix of the 'Goal' vector
diff_partial_mat <- function(alpha, beta, gamma, reg_data){
  x = reg_data['x']
  y = reg_data['y']
  pow_0 = sapply(- exp(alpha + beta * x + gamma * x^2), sum)
  pow_1 = sapply(- x *  exp(alpha + beta * x + gamma * x^2), sum)
  pow_2 = sapply(- x^2 * exp(alpha + beta * x + gamma * x^2), sum)
  pow_3 = sapply(- x^3 * exp(alpha + beta * x + gamma * x^2), sum)
  pow_4 = sapply(- x^4 * exp(alpha + beta * x + gamma * x^2), sum)
  jacobian_f = matrix(c(pow_0, pow_1, pow_2, pow_1, pow_2, pow_3, pow_2, pow_3, pow_4), 
                      nrow=3, ncol=3, byrow = TRUE)
  return(jacobian_f)
}


# Define the Poisson Regression Function
poisreg <- function(alpha_0, beta_0, gamma_0, tolerance, reg_data){
  # Construct vector to restore former theta
  theta_0 <- matrix(c(alpha_0, beta_0, gamma_0), nrow=3, ncol=1, byrow = FALSE)
  vec_f = partial_mat(alpha_0, beta_0, gamma_0, reg_data)
  # Check if the error is less than the tolerance
  if (abs(vec_f[1])<tolerance & abs(vec_f[2])<tolerance & abs(vec_f[3])<tolerance){
    return(theta_0)
  }
  # If not, go to the next iteration
  else{
    jacobian_f = diff_partial_mat(alpha_0, beta_0, gamma_0, reg_data)
    delta_x = solve(jacobian_f, -vec_f)
    alpha = delta_x[1] + alpha_0
    beta = delta_x[2] + beta_0
    gamma = delta_x[3] + gamma_0
    poisreg(alpha, beta, gamma, tolerance, reg_data)
  }
}

poisreg(alpha_0, beta_0, gamma_0, tolerance, PoisRegData)



### Question 3: LOGISTIC REGRESSION--NEWTON'S METHOD

# Read Poisson Regression Data
LogitRegData = read.delim("/Users/elainexfff_/Documents/STAT3006/Assignment 1/Coding_assignment_1/LogitRegData.txt", 
                         header = TRUE, sep = " ")

# Newton's Iteration for Logistic Regression
# Argument: Initial guess alpha_0, beta_0 and error tolerance
tolerance = 0.00001
alpha_0 = 1
beta_0 = 1

# Create the 'Goal' vector (to be optimized)
partial_mat <- function(alpha, beta, reg_data){
  x = reg_data['x']
  y = reg_data['y']
  partial_alpha = sapply(y - exp(alpha + beta * x)/(1 + exp(alpha + beta * x)), sum)
  partial_beta = sapply(x * y - x * exp(alpha + beta * x)/(1 + exp(alpha + beta * x)), sum)
  vec_f = matrix(c(partial_alpha, partial_beta), 
                 nrow=2, ncol=1, byrow = FALSE)
  return(vec_f)
}

#  Calculate the Jacobian matrix of the 'Goal' vector
diff_partial_mat <- function(alpha, beta, reg_data){
  x = reg_data['x']
  y = reg_data['y']
  pow_0 = sapply(- exp(alpha + beta * x)/(1 + exp(alpha + beta * x))^2, sum)
  pow_1 = sapply(- x *  exp(alpha + beta * x)/(1 + exp(alpha + beta * x))^2, sum)
  pow_2 = sapply(- x^2 * exp(alpha + beta * x)/(1 + exp(alpha + beta * x))^2, sum)
  jacobian_f = matrix(c(pow_0, pow_1, pow_1, pow_2), 
                      nrow=2, ncol=2, byrow = TRUE)
  return(jacobian_f)
}

# Define the Logistic Regression Function
logitreg <- function(alpha_0, beta_0, tolerance, reg_data){
  # Construct vector to restore former theta
  theta_0 <- matrix(c(alpha_0, beta_0), nrow=2, ncol=1, byrow = FALSE)
  vec_f = partial_mat(alpha_0, beta_0, reg_data)
  # Check if the error is less than the tolerance
  if ((abs(vec_f[1])<tolerance) && (abs(vec_f[2])<tolerance)){
    return(theta_0)
  }
  # If not, go to the next iteration
  else{
    jacobian_f = diff_partial_mat(alpha_0, beta_0, reg_data)
    delta_x = solve(jacobian_f, -vec_f)
    alpha = delta_x[1] + alpha_0
    beta = delta_x[2] + beta_0
    logitreg(alpha, beta, tolerance, reg_data)
  }
}

logitreg(alpha_0, beta_0, tolerance, LogitRegData)



### Question 4: EM - ALGORITHM

# Read Salary Data
SalaryData = read.delim("/Users/elainexfff_/Documents/STAT3006/Assignment 1/Coding_assignment_1/SalaryData.txt", 
                          header = TRUE, sep = " ")
train_data = SalaryData['x']

# Initial guess of parameters pi, mu, and sigma
pi1_0 = 0.5
pi2_0 = 0.25
mu1_0 = 3000
mu2_0 = 8000
mu3_0 = 30000
sigma1_0 = 300
sigma2_0 = 1500
sigma3_0 = 8000

# Stopping criterion
tolerance = 0.0001

# E-step
z_estimation <- function(pi_1, pi_2, mu_1, mu_2, mu_3, sigma_1, sigma_2, sigma_3, salary_data){
  y = salary_data
  
  p_1 = (pi_1/(sqrt(2*pi) * sigma_1)) * exp(-(y - mu_1)^2/(2 * sigma_1^2))
  p_2 = (pi_2/(sqrt(2*pi) * sigma_2)) * exp(-(y - mu_2)^2/(2 * sigma_2^2)) 
  p_3 = ((1 - pi_1 - pi_2)/(sqrt(2*pi) * sigma_3)) * exp(-(y - mu_3)^2/(2 * sigma_3^2))
  
  # Q-function
  estimated_z1 = p_1/(p_1 + p_2 + p_3)
  estimated_z2 = p_2/(p_1 + p_2 + p_3)
  estimated_z3 = p_3/(p_1 + p_2 + p_3)
  
  return(c(estimated_z1, estimated_z2, estimated_z3))
}

# Compute the observed- data likelihood function
obs_data_likelihood <- function(pi_1, pi_2, mu_1, mu_2, mu_3, sigma_1, sigma_2, sigma_3, y){
  p_1 = (pi_1/(sqrt(2*pi) * sigma_1)) * exp(-(y - mu_1)^2/(2 * sigma_1^2))
  p_2 = (pi_2/(sqrt(2*pi) * sigma_2)) * exp(-(y - mu_2)^2/(2 * sigma_2^2)) 
  p_3 = ((1 - pi_1 - pi_2)/(sqrt(2*pi) * sigma_3)) * exp(-(y - mu_3)^2/(2 * sigma_3^2))
  
  p = sum(p_1 + p_2 + p_3)
  
  return(p)
}

# M-step
maximization <- function(pi_1, pi_2, mu_1, mu_2, mu_3, sigma_1, sigma_2, sigma_3, y, tolerance){
  
  n = length(y) # data size
  z = z_estimation(pi_1, pi_2, mu_1, mu_2, mu_3, sigma_1, sigma_2, sigma_3, y)
  z1 = unlist(z[1])
  z2 = unlist(z[2])
  z3 = unlist(z[3])

  # Update the parameters
  new_pi_1 = sum(z1) / sum(z1 + z2 + z3)
  new_pi_2 = sum(z2) / sum(z1 + z2 + z3)
  new_mu_1 = sum(z1 * y) / sum(z1)
  new_mu_2 = sum(z2 * y) / sum(z2)
  new_mu_3 = sum(z3 * y) / sum(z3)
  new_sigma_1 = sqrt(sum(z1  * (y - new_mu_1)^2) / sum(z1))
  new_sigma_2 = sqrt(sum(z2  * (y - new_mu_2)^2) / sum(z2))
  new_sigma_3 = sqrt(sum(z3  * (y - new_mu_3)^2) / sum(z3))

  # check if it should stop by comparing the difference of observed-data likelihood function of two adjacent iterations
  obs_data = obs_data_likelihood(pi_1, pi_2, mu_1, mu_2, mu_3, sigma_1, sigma_2, sigma_3, y)
  new_obs_data = obs_data_likelihood(new_pi_1, new_pi_2, new_mu_1, new_mu_2, new_mu_3, new_sigma_1, new_sigma_2, new_sigma_3, y)
  if (abs(obs_data - new_obs_data) < tolerance){
    # list out the first 50 classification of individuals
    df <- data.frame (low_income_1 = z1[1:50],
                      middle_income_2 = z2[1:50],
                      high_income_3 = z3[1:50],
                      class = 0
    )
    df['class'] = apply(df,1,function(x) which(x==max(x)))
    print(df)
    
    return(c(new_pi_1, new_pi_2, new_mu_1, new_mu_2, new_mu_3, new_sigma_1, new_sigma_2, new_sigma_3))
  }
  
  maximization(new_pi_1, new_pi_2, new_mu_1, new_mu_2, new_mu_3, new_sigma_1, new_sigma_2, new_sigma_3, y, tolerance)
}

maximization(pi1_0, pi2_0, mu1_0, mu2_0, mu3_0, sigma1_0, sigma2_0, sigma3_0, train_data, tolerance)






