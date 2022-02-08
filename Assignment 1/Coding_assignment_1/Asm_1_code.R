### Question 1: BISECTION METHOD

# Define the function f(x)
f_func <- function(x) {
  return(x^3 + 6 * x^2 + pi * x -12)
}

# Define the bisection method function
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









