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
  group_missing = parLapply(cl, 1:n, group_l)
  grouped_data = do.call(rbind, group_missing)

  denominator = c(grouped_data[,1]+grouped_data[,2]+grouped_data[,3])
  group_data = grouped_data / denominator

  # M-step
  summed_col = colSums(group_data)
  update_parameter_l <- function(x){update_parameter(summed_col, n, x)}
  updated_out = parLapply(cl, 1:3, update_parameter_l)
#  updated_out = lapply(X = 1:3, FUN = update_parameter_l)
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
cl = makeCluster(num_core, type = "FORK")
system.time(maximization_l(pi1_0, pi2_0, mu1_0, mu2_0, mu3_0, sigma1_0, sigma2_0, sigma3_0, train_data, tolerance))
stopCluster(cl)



## Question 2: Database Access from R
rm(list=ls())
#install.packages("RMySQL")
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
#install.packages("DBI")


# (a) find all the companies and ticker symbols
url_complist = "https://www.slickcharts.com/nasdaq100"
doc_complist = htmlTreeParse(rawToChar(GET(url_complist)$content), useInternalNodes = TRUE)

#comp_info = xpathSApply(doc_complist, "//table/tbody/tr/td/a")
comp_info = readHTMLTable(doc_complist, which = 1)
comp = comp_info[, c(2:3)]

# (b) retrieve Market Cap, Price to Book Value, and Dividend Yield from Y-Charts
comp$url = paste0("https://ycharts.com/companies/", comp$Symbol)

empty_cols = c("MarketCap", "PriceToBookValue", "DividendYield")
comp[, empty_cols] = NA

n = nrow(comp)
root = "/Users/elainexfff_/Documents/STAT3006/Assignment 4/Coding_assignment_4/"
file_address = paste0(root, comp[, 'Symbol'], ".html")
for (i in 1:n) {
  url_testlist = comp[i, 'url']
  destfile = file_address[i]
  # download the .html files to local disk if there's no target file
  if(!file.exists(destfile)){
    download.file(url_testlist, destfile)
  }
  doc_file = htmlTreeParse(destfile, useInternalNodes = TRUE)
  
  market_temp = xpathSApply(doc_file, "//div[1]/table/tbody[1]/tr[1]/td[2]/text()", xmlValue)[1]
  if (!is.na(market_temp)){
    market_cap_temp = gsub("\n", "", market_temp) # delete '\n' in the string
    market_cap = gsub(" ","",market_cap_temp) # delete the white-space
    comp[i, 'MarketCap'] = market_cap
  }
  
  price_temp = xpathSApply(doc_file, "//div[1]/table/tbody[2]/tr[4]/td[2]/text()", xmlValue)
  if (!is.na(price_temp)){
    price_value_temp = gsub("\n", "", price_temp) # delete '\n' in the string
    price_to_book_value = gsub(" ", "", price_value_temp) # delete the white-space
    comp[i, 'PriceToBookValue'] = price_to_book_value
  }
  
  dividend_temp = xpathSApply(doc_file, "//div[2]/table/tbody[1]/tr[1]/td[2]/text()", xmlValue)[1]
  if (!is.na(dividend_temp)){
    dividned_yield_temp = gsub("\n", "", dividend_temp) # delete '\n' in the string
    dividend_yield = gsub(" ","",dividned_yield_temp) # delete the white-space
    comp[i, 'DividendYield'] = dividend_yield
  }
}

col_to_export = c("Company", "Symbol", "MarketCap", "PriceToBookValue", "DividendYield")
data_to_export = comp[, col_to_export]

# export the table out
library(gridExtra)
png("Q3b_data.png", height = 50*nrow(data_to_export), width = 200*ncol(data_to_export))
grid.table(data_to_export)
dev.off()

# (c) list the yop 3 companies with the highest Market Cap
mc_cols = c("Company", "Symbol", "MarketCap")
mc_df = comp[, mc_cols]
mc_df[, 'converted_mv'] = NA
unit_1 = 'B'
unit_2 = 'T'
for (i in 1:n) {
  target = mc_df[i, 'MarketCap']
  if (grepl(unit_2, target, fixed = TRUE)){
    mc_df[i, 'converted_mv'] = 1000 * as.numeric(gsub("T", "", target))
  }
  if (grepl(unit_1, target, fixed = TRUE)){
    mc_df[i, 'converted_mv'] = as.numeric(gsub("B", "", target))
  }
}
highest_three = head(mc_df[, 'converted_mv'], 3)
row_ind = which(mc_df$converted_mv == highest_three)
print(mc_df[row_ind,])






