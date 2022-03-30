###################################
#lapply function
###################################
lapply(1:3, function(x){ c(x, x^2, x^3)})

lapply(1:3/3, round, digits = 3)

sapply(1:5, function(x) x^2)

sapply(1:5, function(x) c(x+1, 4*x))

###################################
#parallel package
###################################
library(parallel)

detectCores()

cl <- makeCluster(8)

parLapply(cl, 2:4, function(i) i^2)
parSapply(cl, 1:4, function(i) 2*i)
parSapply(cl, 1:4, function(i) c(100*i, 5*i))
stopCluster(cl) #Never forget to close the cluster if you don't use these cores again.

#elapsed time for a sequential function
system.time(lapply(1:10^7, function(x) x^2))


#elapased time for a parallelized function
cl <- makeCluster(8)
system.time(parLapply(cl, 1:10^7, function(x) x^2))
stopCluster(cl)



###################################
#variable scope
###################################
cl <- makeCluster(4)
base <- 2

parLapply(cl, 2:4, function(exponent){base^exponent})
#Error in checkForRemoteErrors(val) : 
#  3 nodes produced errors; first error: object 'base' not found
stopCluster(cl)


cl <- makeCluster(4)
base <- 2
clusterExport(cl, "base")
parLapply(cl, 2:4, function(exponent){base^exponent})
stopCluster(cl)

#what are the outputs of the following code?
cl <- makeCluster(4)
base <- 2
clusterExport(cl, "base")
base <- 100
parLapply(cl, 2:4, function(exponent){base^exponent})
stopCluster(cl)


###################################
#foreacah package
###################################
#install.packages("foreach")
#install.packages("doParallel")

library(foreach)
library(doParallel)

cl <- makeCluster(4)
registerDoParallel(cl)
base <- 2

foreach(exponent = 2:4, .combine = c) %dopar%{
		base^exponent
	}
	
foreach(exponent = 2:4, .combine = rbind) %dopar%{
		base^exponent
	}
stopCluster(cl)

#
cl <- makeCluster(4)
registerDoParallel(cl)
base <- 2

test <- function (){
	foreach(exponent = 2:4, .combine = c) %dopar%{
		base^exponent
	}
}
test()
#Error in { : task 1 failed - "object 'base' not found"
	
#you can use ".export" in "foreach" function
test <- function (){
	foreach(exponent = 2:4, .combine = c, .export = "base") %dopar%{
		base^exponent
	}
}
test()
stopCluster(cl)

# create output files
cl <- makeCluster(4)
registerDoParallel(cl)
foreach(x = list(1, 2, "a")) %dopar%{
  cat(dput(x), file = paste0("file_", x, ".txt"))
}

stopCluster(cl)

