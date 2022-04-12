## Question 1: Parallel computing for EM alogorithm

rm(list=ls())
#install.packages("parallel")
#install.packages("foreach")
#install.packages("doParallel")
library(parallel)
library(iterators)
library(foreach)
library(doParallel)

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


### original code version
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
  new_sigma_1 = sqrt(sum(z1  * (y - mu_1)^2) / sum(z1))
  new_sigma_2 = sqrt(sum(z2  * (y - mu_2)^2) / sum(z2))
  new_sigma_3 = sqrt(sum(z3  * (y - mu_3)^2) / sum(z3))
  
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
    
    return(c(new_pi_1, new_pi_2, new_mu_1, new_mu_2, new_mu_3, new_sigma_1, new_sigma_2, new_sigma_3))
  }
  
  maximization(new_pi_1, new_pi_2, new_mu_1, new_mu_2, new_mu_3, new_sigma_1, new_sigma_2, new_sigma_3, y, tolerance)
}

system.time(maximization(pi1_0, pi2_0, mu1_0, mu2_0, mu3_0, sigma1_0, sigma2_0, sigma3_0, train_data, tolerance))


### parallel computing version
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
  new_sigma_1 = sqrt(sum(z1  * (y - mu_1)^2) / sum(z1))
  new_sigma_2 = sqrt(sum(z2  * (y - mu_2)^2) / sum(z2))
  new_sigma_3 = sqrt(sum(z3  * (y - mu_3)^2) / sum(z3))
  
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
    
    return(c(new_pi_1, new_pi_2, new_mu_1, new_mu_2, new_mu_3, new_sigma_1, new_sigma_2, new_sigma_3))
  }
  
  maximization(new_pi_1, new_pi_2, new_mu_1, new_mu_2, new_mu_3, new_sigma_1, new_sigma_2, new_sigma_3, y, tolerance)
}

num_core = detectCores() # 8
cl = makeCluster(num_core, type = "FORK")



stopCluster(cl)


## Question 2: Database Access from R
rm(list=ls())
#install.packages("RMySQL")
#install.packages("DBI")
library("RMySQL")

drv=dbDriver("MySQL")
con=dbConnect(drv,user="student", password="HappyStudy2022",
              dbname="Library", port=3306,
              host="rds-mysql-statclass.czyn7pdbk60s.us-west-2.rds.amazonaws.com")

dbListTables(con)
dbReadTable(con,"Student", "Book", "Record")

# (a) list the 'Book' table
dbGetQuery(con,"SELECT * FROM Book;")

# (b) find the ids and the entry years of students 
#     who have borrowed book from natural science
dbGetQuery(con, "SELECT Student.StudentID, EntryYear FROM Student, Record, Book
    where Book.Classification = 'Natural Science'
    And Record.BookNumber = Book.BookNumber
    And Record.StudentID = Student.StudentID;")

# (c) find the ids and majors of the students who have 
#     borrowed book 008 and occupied it for more than 30 days
dbGetQuery(con, "SELECT Student.StudentID, Major FROM Student, Record
           where Record.BookNumber = '8'
           And Record.StudentID = Student.StudentID
           And TIMESTAMPDIFF(day, BorrowingTime, ReturnTime)>30;")
        
dbDisconnect(con)
dbUnloadDriver(drv)


## Question 3: Parse HTML

rm(list=ls())
#install.packages("XML")
#install.packages("httr")
#install.packages("RCurl")
library(XML)
library(httr)
library(RCurl)

# (a) find all the companies and ticker symbols
url_complist = "https://www.slickcharts.com/nasdaq100"
doc_complist = htmlTreeParse(rawToChar(GET(url_complist)$content), useInternalNodes = TRUE)

#comp_info = xpathSApply(doc_complist, "//table/tbody/tr/td/a")
comp_info = readHTMLTable(doc_complist, which = 1)
comp = comp_info[, c(2:3)]

# (b) retrieve Market Cap, Price to Book Value, and Dividend Yield from Y-Charts
url_marklist = "https://ycharts.com/companies/FB"
doc_marklist =  htmlTreeParse(rawToChar(GET(url_marklist)$content), useInternalNodes = TRUE)













