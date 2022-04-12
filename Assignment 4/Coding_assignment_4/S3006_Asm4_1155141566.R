## Question 1:
rm(list=ls())
install.packages("parallel")
library(parallel)

detectCores()

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


## Question 3:
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














