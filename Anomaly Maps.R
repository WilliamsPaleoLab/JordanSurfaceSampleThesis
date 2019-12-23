#Load packages
library(maps)
library(raster)
library(ncdf4)

# Downloads modern pollen and climate data
tmp <- brick("Documents/cru_ts4.02.1901.2017.tmp.dat.nc", varname="tmp") # Mean monthly temperature

sites <- read.csv("Documents/pssSitesClimate.csv", header=TRUE, row.names="site", sep=",")
sites <- sites[,1:2]
sites.tmp <- data.frame(extract(tmp, sites, ncol=2)) # Mean monthly temperature
# Add sample site names to the data.frame
row.names(sites.tmp) <- row.names(sites)
# Change column names
years <- 1951:2015
month <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
names(sites.tmp) <- paste(rep(years, each=12), rep(month, times= 117), sep="_")

# Save the extracted climate data to a .csv file
write.csv(sites.tmp, file="PSS Sites Temp Data.csv")

#read csv back into r after closing out of it 
sites.temperatures <- read.csv("PSS Sites Temp Data.csv", header=TRUE, row.names=1, sep=",", check.names=FALSE) # Mean monthly temperature)

sites.temps <- sites.temperatures[,601:1380] #1951Jan-2015Dec
sites <- read.csv("Documents/pssSitesClimate.csv", header=TRUE, row.names="site", sep=",")
sites <- sites[,1:2]

#next 2 lines same as code below basically
# temp51.80 <- sites.temps[,1:360] #Gets mean monthly temp for Jan1951-Dec1980
# temps86.15 <- sites.temps[,421:780] #Gets mean monthly temp for Jan1986-Dec2015

#From mean monthly temp need to create annual means and then from
#annual means, create 30 year means (or 60 depending on what we need)

#30 year means for January temps
tjan51.80 <- sites.temps[,c(1,13,25,37,49,61,73,85,97,109,121,133,145,157,
                            169,181,193,205,217,229,241,253,265,277,289,
                            301,313,325,337,349)]
tjan51.80$Jan51.80mean <- rowMeans(tjan51.80[,1:30]) #adds means column
tjan86.15 <- sites.temps[,c(421,433,445,457,469,481,
                            493,505,517,529,541,553,565,577,589,601,613,625,
                            637,649,661,673,685,697,709,721,733,745,757,
                            769)]
tjan86.15$Jan86.15mean <- rowMeans(tjan86.15[,1:30]) #adds mean column

#30 year means for July temps
tjul51.80 <- sites.temps[,c(7,19,31,43,55,67,79,91,103,115,127,139,151,163,175,
                            187,199,211,223,235,247,259,271,283,295,307,319,
                            331,343,355)]
tjul51.80$Jul51.80mean <- rowMeans(tjul51.80[,1:30]) #adds mean column

tjul86.15 <- sites.temps[,c(427,439,451,463,475,487,499,511,523,535,
                            547,559,571,583,595,607,619,631,643,655,
                            667,679,691,703,715,727,739,751,763,775)]
tjul86.15$Jul86.15mean <- rowMeans(tjul86.15[,1:30]) #adds mean column

#Create temp avg for all 65 years
tavg <- rowMeans(sites.temps[,1:780])
#Create temp avg for 1951-1980 and 1986-2015
tavg51.80 <- rowMeans(sites.temps[,1:360])
tavg.86.15 <- rowMeans(sites.temps[,421:780])

all.temp.data <- data.frame(sites[,0:2], tjan51.80$Jan51.80mean,
                            tjan86.15$Jan86.15mean, tjul51.80$Jul51.80mean,
                            tjul86.15$Jul86.15mean, tavg, tavg51.80, tavg.86.15)

anom.temps <- data.frame(site = rownames(all.temp.data),
                         temp.anom = all.temp.data[,9] - all.temp.data[,8],
                         lat = all.temp.data$lat,
                         long = all.temp.data$lon)

jan.anom <- data.frame(site = rownames(all.temp.data),
                       januaryanom = all.temp.data[,4] - all.temp.data[,3],
                       lat = all.temp.data$lat,
                       long = all.temp.data$lon)

jul.anom <- data.frame(site = rownames(all.temp.data), 
                       julyanom = all.temp.data[,6] - all.temp.data[,5],
                       lat = all.temp.data$lat,
                       long = all.temp.data$lon)
#Anomaly Maps
library(ggplot2)
worldmap <- map_data("world")
USCanada <- c("United States", "Canada")
map_USCan <- map_data("world", region = USCanada)

gg <- ggplot(worldmap, aes(x = long, y = lat)) +
               geom_polygon(aes(group = group), color = "black", alpha = 0) + 
  geom_point(data = anom.temps, aes(x = long, y = lat, color = temp.anom)) +
  scale_color_gradient2(low="blue", mid="white", high="red") +
  coord_cartesian(xlim = c(-40, -175), ylim = c(20, 85)) #specifies a region
gg
                      
#Jan86.15.mean - Jan51.80.mean
gg.january <- ggplot(worldmap, aes(x = long, y = lat)) +
  geom_polygon(aes(group = group), color = "black", alpha = 0) +
  geom_point(data = jan.anom, aes(x = long, y = lat, color = januaryanom)) +
  scale_color_gradient2(low="blue", mid="white", high="red") +
  coord_cartesian(xlim = c(-40, -175), ylim = c(20, 85)) 
gg.january

#July86.15.mean - July51.80.mean
gg.july <- ggplot(worldmap, aes(x = long, y = lat)) +
  geom_polygon(aes(group = group), color = "black", alpha = 0) +
  geom_point(data = jul.anom, aes(x = long, y = lat, color = julyanom)) +
  scale_color_gradient2(low="blue", mid="white", high="red") + 
  coord_cartesian(xlim = c(-40, -175), ylim = c(20, 85))
gg.july


#now need to change colors/make them inverse so that it shows trends of getting warmer
#look into scale_color_gradient2 or scale_color_gradientn for changing the temp scales
#then add the line into each gg plot 

