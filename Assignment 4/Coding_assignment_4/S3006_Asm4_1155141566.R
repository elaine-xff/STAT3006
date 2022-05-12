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

gp1 <- function(pi_1, mu_1, sigma_1, y){
  p_1 = (pi_1/(sqrt(2*pi) * sigma_1)) * exp(-(y - mu_1)^2/(2 * sigma_1^2))
  p_1_y = p_1 * y
  p_1_y2 = p_1 * (y-mu_1)^2
  
  return(c(p_1, p_1_y, p_1_y2))
}

gp2 <- function(pi_2, mu_2, sigma_2, y){
  p_2 = (pi_2/(sqrt(2*pi) * sigma_2)) * exp(-(y - mu_2)^2/(2 * sigma_2^2))
  p_2_y = p_2 * y
  p_2_y2 = p_2 * (y-mu_2)^2
  
  return(c(p_2, p_2_y, p_2_y2))
}

gp3 <- function(pi_1, pi_2, mu_3, sigma_3, y){
  p_3 = ((1 - pi_1 - pi_2)/(sqrt(2*pi) * sigma_3)) * exp(-(y - mu_3)^2/(2 * sigma_3^2))
  p_3_y = p_3 * y
  p_3_y2 = p_3 * (y-mu_3)^2
  
  return(c(p_3, p_3_y, p_3_y2))
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
  
  n = length(salary_data) # data size
  
  # E-step and pre-computation of intermediate parameters
  gp1_l <- function(x) {gp1(pi_1, mu_1, sigma_1, x)}
  gp1_missing = lapply(X = salary_data, FUN = gp1_l)
  gp1_data = do.call(rbind, gp1_missing)
  
  gp2_l <- function(x) {gp2(pi_2, mu_2, sigma_2, x)}
  gp2_missing = lapply(X = salary_data, FUN = gp2_l)
  gp2_data = do.call(rbind, gp2_missing)
  
  gp3_l <- function(x) {gp3(pi_1, pi_2, mu_3, sigma_3, x)}
  gp3_missing = lapply(X = salary_data, FUN = gp3_l)
  gp3_data = do.call(rbind, gp3_missing)
  
  denominator = gp1_data[,1] + gp2_data[,1] + gp3_data[,1]
  
  group_1 = gp1_data/denominator
  group_2 = gp2_data/denominator
  group_3 = gp3_data/denominator
  
  # M-step
  new_pi_1 = sum(group_1[,1]) / n
  new_mu_1 = sum(group_1[,2]) / sum(group_1[,1])
  new_sigma_1 = sqrt(sum(group_1[,3]) / sum(group_1[,1]))
  
  new_pi_2 = sum(group_2[,1]) / n
  new_mu_2 = sum(group_2[,2]) / sum(group_2[,1])
  new_sigma_2 = sqrt(sum(group_2[,3]) / sum(group_2[,1]))
  
  new_mu_3 = sum(group_3[,2]) / sum(group_3[,1])
  new_sigma_3 = sqrt(sum(group_3[,3]) / sum(group_3[,1]))
  
  # check if it should stop by comparing the difference of observed-data likelihood function of two adjacent iterations
  y=salary_data
  obs_data = obs_data_likelihood(pi_1, pi_2, mu_1, mu_2, mu_3, sigma_1, sigma_2, sigma_3, y)
  new_obs_data = obs_data_likelihood(new_pi_1, new_pi_2, new_mu_1, new_mu_2, new_mu_3, new_sigma_1, new_sigma_2, new_sigma_3, y)
  if (abs(obs_data - new_obs_data) < tolerance){
    return(c(new_pi_1, new_pi_2, new_mu_1, new_mu_2, new_mu_3, new_sigma_1, new_sigma_2, new_sigma_3))
  }
  
  maximization(new_pi_1, new_pi_2, new_mu_1, new_mu_2, new_mu_3, new_sigma_1, new_sigma_2, new_sigma_3, y, tolerance)
}
# parallel computing version
maximization_l <- function(pi_1, pi_2, mu_1, mu_2, mu_3, sigma_1, sigma_2, sigma_3, salary_data, tolerance){
  
  n = length(salary_data) # data size
  
  # E-step and pre-computation of intermediate parameters
  gp1_l <- function(x) {gp1(pi_1, mu_1, sigma_1, x)}
  gp1_missing = mclapply(X = salary_data, FUN = gp1_l, mc.cores = num_core)
  gp1_data = do.call(rbind, gp1_missing)
  
  gp2_l <- function(x) {gp2(pi_2, mu_2, sigma_2, x)}
  gp2_missing = mclapply(X = salary_data, FUN = gp2_l, mc.cores = num_core)
  gp2_data = do.call(rbind, gp2_missing)
  
  gp3_l <- function(x) {gp3(pi_1, pi_2, mu_3, sigma_3, x)}
  gp3_missing = mclapply(X = salary_data, FUN = gp3_l, mc.cores = num_core)
  gp3_data = do.call(rbind, gp3_missing)
  
  denominator = gp1_data[,1] + gp2_data[,1] + gp3_data[,1]
  
  group_1 = gp1_data/denominator
  group_2 = gp2_data/denominator
  group_3 = gp3_data/denominator
  
  # M-step
  new_pi_1 = sum(group_1[,1]) / n
  new_mu_1 = sum(group_1[,2]) / sum(group_1[,1])
  new_sigma_1 = sqrt(sum(group_1[,3]) / sum(group_1[,1]))
  
  new_pi_2 = sum(group_2[,1]) / n
  new_mu_2 = sum(group_2[,2]) / sum(group_2[,1])
  new_sigma_2 = sqrt(sum(group_2[,3]) / sum(group_2[,1]))
  
  new_mu_3 = sum(group_3[,2]) / sum(group_3[,1])
  new_sigma_3 = sqrt(sum(group_3[,3]) / sum(group_3[,1]))
  
  # check if it should stop by comparing the difference of observed-data likelihood function of two adjacent iterations
  y=salary_data
  obs_data = obs_data_likelihood(pi_1, pi_2, mu_1, mu_2, mu_3, sigma_1, sigma_2, sigma_3, y)
  new_obs_data = obs_data_likelihood(new_pi_1, new_pi_2, new_mu_1, new_mu_2, new_mu_3, new_sigma_1, new_sigma_2, new_sigma_3, y)
  if (abs(obs_data - new_obs_data) < tolerance){
    return(c(new_pi_1, new_pi_2, new_mu_1, new_mu_2, new_mu_3, new_sigma_1, new_sigma_2, new_sigma_3))
  }
  
  maximization_l(new_pi_1, new_pi_2, new_mu_1, new_mu_2, new_mu_3, new_sigma_1, new_sigma_2, new_sigma_3, y, tolerance)
}
# parallel version
num_core = detectCores() # 8
cl = makeCluster(num_core, type = "FORK")
system.time(maximization_l(pi1_0, pi2_0, mu1_0, mu2_0, mu3_0, sigma1_0, sigma2_0, sigma3_0, train_data, tolerance))
stopCluster(cl)  

# original version
system.time(maximization(pi1_0, pi2_0, mu1_0, mu2_0, mu3_0, sigma1_0, sigma2_0, sigma3_0, train_data, tolerance))







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












