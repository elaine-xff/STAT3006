## Question 1: Parallel computing for EM alogorithm
rm(list=ls())
library(parallel)

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

group <- function(pi_1, pi_2, mu_1, mu_2, mu_3, sigma_1, sigma_2, sigma_3, salary_data, i){
  y=salary_data[i,1]

  p_1 = (pi_1/(sqrt(2*pi) * sigma_1)) * exp(-(y - mu_1)^2/(2 * sigma_1^2))
  p_2 = (pi_2/(sqrt(2*pi) * sigma_2)) * exp(-(y - mu_2)^2/(2 * sigma_2^2))
  p_3 = ((1 - pi_1 - pi_2)/(sqrt(2*pi) * sigma_3)) * exp(-(y - mu_3)^2/(2 * sigma_3^2))
  p = c(p_1, p_2, p_3)
  
  return(c(p, y*p, p_1*(y-mu_1)^2, p_2*(y-mu_2)^2, p_3*(y-mu_3)^2))
}

update_parameter <-function(summed_col, n, i){
  return(c(summed_col[i]/n, summed_col[i+3]/summed_col[i], sqrt(summed_col[i+6] / summed_col[i])))
}

obs_data_likelihood <- function(pi_1, pi_2, mu_1, mu_2, mu_3, sigma_1, sigma_2, sigma_3, y){
  p_1 = (pi_1/(sqrt(2*pi) * sigma_1)) * exp(-(y - mu_1)^2/(2 * sigma_1^2))
  p_2 = (pi_2/(sqrt(2*pi) * sigma_2)) * exp(-(y - mu_2)^2/(2 * sigma_2^2)) 
  p_3 = ((1 - pi_1 - pi_2)/(sqrt(2*pi) * sigma_3)) * exp(-(y - mu_3)^2/(2 * sigma_3^2))
  
  p = sum(p_1 + p_2 + p_3)
  return(p)
}
# original version
maximization <- function(pi_1, pi_2, mu_1, mu_2, mu_3, sigma_1, sigma_2, sigma_3, salary_data, tolerance){
  
  n = nrow(salary_data) # data size
  
  # E-step and pre-computation of intermediate parameters
  group_l <- function(x){group(pi_1, pi_2, mu_1, mu_2, mu_3, sigma_1, sigma_2, sigma_3, salary_data, x)}
  group_missing = lapply(X = 1:n, FUN = group_l)
  grouped_data = do.call(rbind, group_missing)

  denominator = c(grouped_data[,1]+grouped_data[,2]+grouped_data[,3])
  group_data = grouped_data / denominator

  
  # M-step
  summed_col = colSums(group_data)
  
  update_parameter_l <- function(x){update_parameter(summed_col, n, x)}
  updated_out = lapply(X = 1:3, FUN = update_parameter_l)
  param_out = do.call(rbind, updated_out)

  new_pi_1 = param_out[1,1]
  new_mu_1 = param_out[1,2]
  new_sigma_1 = param_out[1,3]
  
  new_pi_2 = param_out[2,1]
  new_mu_2 = param_out[2,2]
  new_sigma_2 = param_out[2,3]
  
  new_mu_3 = param_out[3,2]
  new_sigma_3 = param_out[3,3]
  # check if it should stop by comparing the difference of observed-data likelihood function of two adjacent iterations
  obs_data = sum(denominator)
  new_obs_data = obs_data_likelihood(new_pi_1, new_pi_2, new_mu_1, new_mu_2, new_mu_3, new_sigma_1, new_sigma_2, new_sigma_3, salary_data)
  if (abs(obs_data - new_obs_data) < tolerance){
    return(c(new_pi_1, new_pi_2, new_mu_1, new_mu_2, new_mu_3, new_sigma_1, new_sigma_2, new_sigma_3))
  }

  maximization(new_pi_1, new_pi_2, new_mu_1, new_mu_2, new_mu_3, new_sigma_1, new_sigma_2, new_sigma_3, salary_data, tolerance)
}

# parallel computing version
maximization_l <- function(pi_1, pi_2, mu_1, mu_2, mu_3, sigma_1, sigma_2, sigma_3, salary_data, tolerance){
  
  n = nrow(salary_data) # data size
  
  # E-step and pre-computation of intermediate parameters
  group_l <- function(x){group(pi_1, pi_2, mu_1, mu_2, mu_3, sigma_1, sigma_2, sigma_3, salary_data, x)}
  group_missing = mclapply(X = 1:n, FUN = group_l, mc.cores = num_core/4)
  grouped_data = do.call(rbind, group_missing)

  denominator = c(grouped_data[,1]+grouped_data[,2]+grouped_data[,3])
  group_data = grouped_data / denominator

  # M-step
  summed_col = colSums(group_data)
  update_parameter_l <- function(x){update_parameter(summed_col, n, x)}
  updated_out = lapply(X = 1:3, FUN = update_parameter_l)
#  updated_out = mclapply(X = 1:3, FUN = update_parameter_l, mc.cores = )
  param_out = do.call(rbind, updated_out)
  
  new_pi_1 = param_out[1,1]
  new_mu_1 = param_out[1,2]
  new_sigma_1 = param_out[1,3]
  
  new_pi_2 = param_out[2,1]
  new_mu_2 = param_out[2,2]
  new_sigma_2 = param_out[2,3]
  
  new_mu_3 = param_out[3,2]
  new_sigma_3 = param_out[3,3]
  
  # check if it should stop by comparing the difference of observed-data likelihood function of two adjacent iterations
  obs_data = sum(denominator)
  new_obs_data = obs_data_likelihood(new_pi_1, new_pi_2, new_mu_1, new_mu_2, new_mu_3, new_sigma_1, new_sigma_2, new_sigma_3, salary_data)
  if (abs(obs_data - new_obs_data) < tolerance){
    return(c(new_pi_1, new_pi_2, new_mu_1, new_mu_2, new_mu_3, new_sigma_1, new_sigma_2, new_sigma_3))
  }
  
  maximization_l(new_pi_1, new_pi_2, new_mu_1, new_mu_2, new_mu_3, new_sigma_1, new_sigma_2, new_sigma_3, salary_data, tolerance)
}

# original version
system.time(maximization(pi1_0, pi2_0, mu1_0, mu2_0, mu3_0, sigma1_0, sigma2_0, sigma3_0, train_data, tolerance))

# parallel version
num_core = detectCores()
cl = makeCluster(num_core/2, type = "FORK")
system.time(maximization_l(pi1_0, pi2_0, mu1_0, mu2_0, mu3_0, sigma1_0, sigma2_0, sigma3_0, train_data, tolerance))
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
#url_marklist = "https://ycharts.com/companies/FB"
#destfile = "/Users/elainexfff_/Documents/STAT3006/Assignment 4/Coding_assignment_4/FB.html"
#download.file(url_marklist, destfile)
#doc_marklist = htmlTreeParse(destfile, useInternalNodes = TRUE)

url_marketcap = "https://ycharts.com/companies/FB/market_cap"
destfile_mc = "/Users/elainexfff_/Documents/STAT3006/Assignment 4/Coding_assignment_4/FB_market_cap.html"
download.file(url_marketcap, destfile_mc)
doc_marketcap = htmlTreeParse(destfile_mc, useInternalNodes = TRUE)
marketcap_table = readHTMLTable(doc_marketcap, which = 3)
colnames(marketcap_table) = c("Company", "MarketCap")
temp1 = merge(comp, marketcap_table, by=c("Company"), all.x = TRUE)
  
url_pricetbv = "https://ycharts.com/companies/FB/price_to_book_value"
destfile_pricetbv = "/Users/elainexfff_/Documents/STAT3006/Assignment 4/Coding_assignment_4/FB_price_to_book_value.html"
download.file(url_pricetbv, destfile_pricetbv)
doc_pricetbv = htmlTreeParse(destfile_pricetbv, useInternalNodes = TRUE)
pricetbv_table = readHTMLTable(doc_pricetbv, which = 3)
colnames(pricetbv_table) = c("Company", "PriceToBookValue")
temp2 = merge(temp1, pricetbv_table, by=c("Company"), all.x = TRUE)

url_dy = "https://ycharts.com/companies/FB/dividend_yield"
destfile_dy = "/Users/elainexfff_/Documents/STAT3006/Assignment 4/Coding_assignment_4/FB_dividend_yield.html"
download.file(url_dy, destfile_dy)
doc_dy = htmlTreeParse(destfile_dy, useInternalNodes = TRUE)
dy_table = readHTMLTable(doc_dy, which = 3)
colnames(dy_table) = c("Company", "DividendYield")
whole_table = merge(temp2, dy_table, by=c("Company"), all.x = TRUE)

#install.packages("gridExtra")
library("gridExtra")
pdf("Q3b_data.pdf")       # Export PDF
grid.table(whole_table)
dev.off()












