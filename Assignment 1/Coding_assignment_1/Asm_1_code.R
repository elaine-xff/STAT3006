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
alpha_0 = 10
beta_0 = 44
gamma_0 = 50

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

library(dplyr)

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





