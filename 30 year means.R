# Extracting CRU Climate data: CRU TS v4.01
# Complete guide available at: http://www.benjaminbell.co.uk
# Load packages
library(raster)
library(ncdf4)
# Load the CRU TS datasets into R 
tmp <- brick("Documents/cru_ts4.02.1901.2017.tmp.dat.nc", varname="tmp") # Mean monthly temperature

# Import sample site information
samples <- read.csv("Desktop/samples.csv", header=TRUE, row.names="site", sep=",")

# Extract climate data from the RasterBrick as a data.frame

tmp.sites <- data.frame(extract(tmp, samples, ncol=2)) # Mean monthly temperature
# Add sample site names to the data.frame
row.names(tmp.sites) <- row.names(samples)
# Change column names
years <- 1901:2017
month <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
names(tmp.sites) <- paste(rep(years, each=12), rep(month, times= 117), sep="_")
# Save the extracted climate data to a .csv file

write.csv(tmp.sites, file="Temperature Data.csv")

#read csv back into r after closing out of it 
temp <- read.csv("Temperature Data.csv", header=TRUE, row.names=1, sep=",", check.names=FALSE) # Mean monthly temperature)

years <- 1950:2017

# Calculate mean annual temperature values
temp.year.mean <- as.data.frame(sapply(years, function(x) rowMeans(temp[, grep(x, names(temp))])))
names(temp.year.mean) <- years # Rename columns in the new data frame object

#30 year mean for 1950-1980 - will probs be useful for anom maps 
#See the length of the colnames or how long the range is
which(colnames(temp)=="1950_Jan")
which(colnames(temp)=="1979_Dec")
temps1950.1980 <- temp[589:948]
# Calculate mean monthly data
temp1950.1980.mean <- as.data.frame(sapply(month, function(x) rowMeans(temps1950.1980[, grep(x, names(temps1950.1980))])))
names(temp1950.1980.mean) <- month # Rename columns in the new data frame object

#Actual 30yr mean 1950-80 ; gives us each monthly avg
years50.80 <- 1950:1979
mean.1950.80 <- as.data.frame(sapply(years50.80, function(x) rowMeans(temps1950.1980[, grep(x, names(temp))])))
names(mean.1950.80) <- years50.80

