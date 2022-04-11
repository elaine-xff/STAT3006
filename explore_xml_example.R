#install.packages("XML")
library(XML)
setwd("E:/CUHK/teaching/STAT3006")

#############################################################
#a plant catalog example
#place plant_catalog.xml in your working directory
#############################################################
#read in the xml file
doc = xmlTreeParse('plant_catalog.xml')
root = xmlRoot(doc)
class(root)
#take a look the the sample elements
head(root)
#the elements of the first layer
names(root)
#the total number of first layer elements
table(names(root))

#look into the first element
oneplant=root[[1]]
#the value of the "COMMON" tag for the first element
xmlValue(oneplant[['COMMON']])
#obtain the values of the "COMMON" tag for all the elements
commons=xmlSApply(root, function(x) xmlValue(x[['COMMON']]))
#construct a function to obtain a given tag "var" of the element
getvar=function(x,var) xmlValue(x[[var]])
#collect all the attributes of each element into a table
res = lapply(names(root[[1]]),
	function(var)xmlSApply(root,getvar,var))
plants=data.frame(res)

#############################################################
#a nested structure example
#location of counties within each state
#place counties.gml.xml in your working directory
#############################################################
rm(list=ls())
#read in the xml file
doc = xmlTreeParse('counties.gml.xml')
root=xmlRoot(doc)
#take a look the the sample elements
head(root)
#total number of states
table(names(root))
#retrieve the names for all the first layer elements
statenames = xmlSApply(root,function(x)xmlValue(x[['name']]))
head(statenames)
#retrieve the  name-value pairs of the "name" attribute 
#(i.e. abbreviation for states in the example) 
#for all the first layer elements
stateabbs = xmlSApply(root,function(x)xmlAttrs(x[['name']]))
head(stateabbs)

#look into the first element
onestate=root[[1]]
#count the number of second layer elements for 
#the first element of the first layer
table(names(onestate))
#take out the elements under the "county" tags
counties = xmlElementsByTagName(onestate,'county')
#xmlElementsByTagName returns a list
class(counties)
length(counties)

#construct a function to collect the second layer elements for
#a given first layer element
SingleState = function(state){
counties = xmlElementsByTagName(state,'county')
countynames = sapply(counties,function(x)xmlValue(x[['name']]))
coords = lapply(counties,function(x)x[['location']][['coord']])
x = as.numeric(sapply(coords,function(x)xmlValue(x[['X']])))
y = as.numeric(sapply(coords,function(x)xmlValue(x[['Y']])))
data.frame(county=countynames,x=x,y=y)
}

#survey all the first layer elements (i.e. states in the example)
res = xmlApply(root,SingleState)
names(res) = xmlSApply(root,function(x)xmlValue(x[['name']]))

