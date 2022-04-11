#install.packages("RMySQL")
#install.packages("DBI")
library("RMySQL")

drv=dbDriver("MySQL")
con=dbConnect(drv,user="student", password="HappyStudy2022",
	dbname="Bank", port=3306,
	host="rds-mysql-statclass.czyn7pdbk60s.us-west-2.rds.amazonaws.com")
	
dbListTables(con)
dbReadTable(con,"Account")
dbGetQuery(con,"SELECT * FROM Customer;")

dbGetQuery(con, "SELECT CID, SUM(Balance) AS Total
	FROM Registration, Account
    where Account.AcctNo=Registration.AcctNo
    group by CID;")
	
dbDisconnect(con)
dbUnloadDriver(drv)
