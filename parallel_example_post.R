library(parallel)

sum_split<-function(x,Total,core)
{
	start=Total/core*(x-1)+1
	end=Total/core*(x)
	result<-0
	for(i in start:end)
	{
		result<-result+i
	}
	result
}

sum_up<-function(Total)
{
	for(i in 1:Total)
	{
		result<-result+i
	}
	result
}

# parallel computing
sum_parallel <- function(Total, numcore=detectCores()){
    require(parallel)
    # clusterApply() for Windows
    if (Sys.info()[1] == "Windows"){
	
	cl <- makeCluster(numcore)
		#clusterExport(cl, c("sum_split","numcore","Total"))
		clusterExport(cl, "sum_split")
        runtime <- system.time({
            sumed_up <- clusterApply(cl,1:numcore, function(y) sum_split(y,Total,numcore))
        })[3]
        stopCluster(cl) # Don't forget to do this--I frequently do
 
    # mclapply() for everybody else
    } else {
		sum_split_l<-function(x) sum_split(x,Total,numcore)
		runtime <- system.time({
            sumed_up<-mclapply(X=1:numcore, FUN=sum_split_l,mc.cores=numcore)
        })[3]
    }
    return(list(sumed_up=sum(unlist(sumed_up)), runtime=runtime,numcore=numcore))
}

run_result<-sapply(c(1,2,3,4,8), function(y) sum_parallel(1200000,y))
run_result<-rbind(unlist(run_result[1,]),unlist(run_result[2,]),unlist(run_result[3,]))
plot(run_result[2,]~run_result[3,],xlab="no. of cores", ylab="time")
plot(run_result[2,1]/run_result[2,]~run_result[3,],
	xlim=c(0,8),ylim=c(0,8), xlab="no. of cores", ylab="speed")
abline(0,1)
