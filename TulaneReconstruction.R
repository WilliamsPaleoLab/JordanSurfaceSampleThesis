#Reconstruction for Lake Tulane

library(neotoma)
library(analogue)
library(rioja)
library(vegan)
library(maps)

# Gets historic datasets
tulane_dataset <- get_dataset(19620) #Lake Tulane Dataset ID

# Stores data associated with each Lake Tulane dataset
tulane_data <- get_download(tulane_dataset)

# Assigns the first dataset, which is pollen, to tulane_pollen
tulane_pollen <- tulane_data[[1]]

# Downloads modern pollen and climate data
library(raster)

#read sites temp data.csv back into R
sites.temperatures <- read.csv("PSS Sites Temp Data.csv", header=TRUE, row.names=1, sep=",", check.names=FALSE) # Mean monthly temperature)

sites.temps <- sites.temperatures[,601:1380] #1951Jan-2015Dec
#not sure if these next 2 lines are needed
# temp51.80 <- sites.temps[,1:360] #Gets mean monthly temp for Jan1951-Dec1980
# temps86.15 <- sites.temps[,421:780] #Gets mean monthly temp for Jan1986-Dec2015

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

#Creates data frame with necessary info
#If get an error might have to load sites csv in before all.temp.data
all.temp.data <- data.frame(tjan51.80$Jan51.80mean,
                            tjan86.15$Jan86.15mean, tjul51.80$Jul51.80mean,
                            tjul86.15$Jul86.15mean, tavg, tavg51.80, tavg.86.15) #sites[,0:2]

# Assigns climate variables to vectors
tempavg51.80 <- all.temp.data$tavg51.80
tempavg86.15 <- all.temp.data$tavg.86.15

#Read surface sample data in 
surfacesamples <- read.csv("Documents/actual modern pollen data.csv")[,1:4]
#Download dataset info i.e. taxa counts whatever 
ss.downloaded <- neotoma::get_download(surfacesamples[,1])
clean.surfacesamples <- compile_taxa(ss.downloaded, list.name = "WS64")

#create matrix for modern pollen data to pull from later
mod.pollen <- matrix(ncol = 65, #65 columns is all taxa from ws64 + 'other'
                     nrow = length(clean.surfacesamples), #makes it length 1971
                     data = 0) #makes everything 0 so can put data in later
mod.clim <- matrix(ncol = 5, #climate matrix
                   nrow = length(clean.surfacesamples),
                   data = 0)
mod.pollen <- as.data.frame(mod.pollen)
mod.clim <- as.data.frame(mod.clim)
colnames(mod.clim) <- c("site.name", "lat", "long", "51.80", "86.15") #what we want row names of mod clim to be
alt.table <- read.csv("Downloads/poll_temps/alt_table.csv", stringsAsFactors = FALSE)
colnames(mod.pollen) <- sort(c(unique(alt.table$WS64), "Other")) 

#make a loop to basically import data into the empty matrices
i=1
#val <- holds spot for data to go
#assign val to mod pollen
#when want to add new temp data, add it above in matrix, then add another
#row of data below 
for(i in 1:length(clean.surfacesamples)){
  val <- merge(as.data.frame(clean.surfacesamples[[i]]$counts), mod.pollen[1,], all.x = TRUE)
  mod.pollen[i,] <- val
  mod.clim[i,1] <- clean.surfacesamples[[i]]$dataset$site.data$site.name
  mod.clim[i,2] <- clean.surfacesamples[[i]]$dataset$site.data$lat
  mod.clim[i,3] <- clean.surfacesamples[[i]]$dataset$site.data$long
  mod.clim[i,4] <- tempavg51.80[i]
  mod.clim[i,5] <- tempavg86.15[i]
}
#bind data frames together so everything is in one place to grab data 
mod.data <- cbind(mod.clim, mod.pollen)
mod.data[is.na(mod.data)] <- 0 #set na values to 0
mod.data <- mod.data[,-grep("Other", names(mod.data))] #removes "other" from df

# Compiles tulane pollen data with Williams & Shuman 2008 taxa list
tulane_pollen_clean <- compile_taxa(tulane_pollen, list.name = "WS64")

# Cleans counts & gives us percentages for historic site
tulane_counts <- tulane_pollen_clean$counts
tulane_counts_clean <- replace(tulane_counts, is.na(tulane_counts), 0)
tulane_counts_final <- tulane_counts_clean[,c(-26)]
tulane_counts_final <- (tulane_counts_final/rowSums(tulane_counts_final)) * 100

# Calculates pollen percentages for modern pollen sites 
#Not sure if need this --> Pollen_percentages <- 100*(Pollen_counts/rowSums(Pollen_counts))
mod.data[,6:69] <- (mod.data[,6:69]/rowSums(mod.data[,6:69])) * 100

#First build model then predict temps

#model for MAT 51.80
mat.51.80 <- rioja::MAT(y = mod.data[,6:42], x = mod.data$`51.80`, dist.method = "sq.chord", k = 5, lean = FALSE)
cross.51.80 <- rioja::crossval(mat.51.80, cv.method = "boot", nboot = 100)
tulane.temp.51.80 <- predict(cross.51.80, newdata = tulane_counts_final)
MAT.ts.51.80 <- plot(tulane_pollen[2]$sample.meta$age, tulane.temp.51.80$fit[,1], type = "o", main = 'MAT 1951-1980', xlab = "Ages", ylab = "Predicted Temps")

#Model for MAT 86.15
mat.86.15 <- rioja::MAT(y = mod.data[,6:42], x = mod.data$`86.15`, dist.method = "sq.chord", k = 5, lean = FALSE)
cross.86.15 <- rioja::crossval(mat.86.15, cv.method = "boot", nboot = 100)
tulane.temp.86.15 <- predict(cross.86.15, newdata = tulane_counts_final)
MAT.ts.86.15 <- plot(tulane_pollen[2]$sample.meta$age, tulane.temp.86.15$fit[,1], type = "o", main = 'MAT 1986-2015', xlab = "Ages", ylab = "Predicted Temps (˚Celsius)")

#Model for WA 51.80
wa.51.80 <- rioja::WA(y = mod.data[,6:42], x = mod.data$`51.80`, lean = FALSE) # has to be 6:42 otherwise get error about species data have 0 abundances for following columns
cross.wa.51.80 <- rioja::crossval(wa.51.80, cv.method = "boot", nboot = 100)
wa.tulane.temp.51.80 <- predict(cross.wa.51.80, newdata = tulane_counts_final)
WA.ts.51.80 <- plot(tulane_pollen[2]$sample.meta$age, wa.tulane.temp.51.80$fit[,1], type = "o", main = 'WA 1951-1980', xlab = "Ages", ylab = "Predicted Temps")

#Model for WA 86.15
wa.86.15 <- rioja::WA(y = mod.data[,6:42], x = mod.data$`86.15`, lean = FALSE)
cross.wa.86.15 <- rioja::crossval(wa.86.15, cv.methods = "boot", nboot = 100)
wa.tulane.temp.86.15 <- predict(cross.wa.86.15, newdata= tulane_counts_final)
WA.ts.86.15 <- plot(tulane_pollen[2]$sample.meta$age, wa.tulane.temp.86.15$fit[,1], type = "o", main = 'WA 1986-2015', xlab = "Ages", ylab = "Predicted Temps")

#Model for WAPLS 51.80
wapls.51.80 <- rioja::WAPLS(y = mod.data[,6:42], x = mod.data$`51.80`, lean = FALSE)
cross.wapls.51.80 <- rioja::crossval(wapls.51.80, cv.methods = "boot", nboot = 100)
wapls.tulane.temp51.80 <- predict(cross.wapls.51.80, newdata = tulane_counts_final)
WAPLS.ts.51.80 <- plot(tulane_pollen[2]$sample.meta$age, wapls.tulane.temp51.80$fit[,1], type = "o", main = 'WAPLS 1951-1980', xlab = 'Ages', ylab = 'Predicted Temps')

#Model for WAPLS 86.15
wapls.86.15 <- rioja::WAPLS(y=mod.data[,6:42], x = mod.data$`86.15`, lean = FALSE)
cross.wapls.86.15 <- rioja::crossval(wapls.86.15, cv.methods = "boot", nboot = 100)
wapls.tulane.temp86.15 <- predict(cross.wapls.86.15, newdata = tulane_counts_final)
WAPLS.ts.86.15 <- plot(tulane_pollen[2]$sample.meta$age, wapls.tulane.temp86.15$fit[,1], type = "o", main = 'WAPLS 1986-2015', xlab = 'Ages', ylab = 'Predicted Temps')


# #use square root pollen percentages for each model(TF)???????
# #keep Pollen_percentages where it is 
# # Runs MAT on modern pollen percentages and 1951-80 30-year mean
# MAT.tempavg51.80 <- rioja::MAT(Pollen_percentages, tempavg51.80, lean = FALSE)
# # Plots MAT results
# plot(MAT.tempavg51.80)
# # cross-validates the MAT using the leave-one-out technique
# cv.tempavg51.80.mat.model <- rioja::crossval(MAT.tempavg51.80, cv.method='lgo', verbose=FALSE)
# plot(cv.tempavg51.80.mat.model)

# # Runs MAT on modern pollen percentages and 1986-2015 30-year mean
# MAT.tempavg86.15 <- rioja::MAT(Pollen_percentages, tempavg86.15, lean = FALSE)
# # Plots MAT results
# plot(MAT.tempavg86.15)
# # cross-validates the MAT using the leave-one-out technique
# cv.tempavg86.15.mat.model <- rioja::crossval(MAT.tempavg86.15, cv.method='lgo', verbose=FALSE)
# plot(cv.tempavg86.15.mat.model)


#data frame to compare model temps
mat.temps.86.15 <- tulane.temp.86.15$fit[,1]
mat.temps.51.80 <- tulane.temp.51.80$fit[,1]
wapls.temps.51.80 <- wapls.tulane.temp51.80$fit[,1]
wapls.temps.86.15 <- wapls.tulane.temp86.15$fit[,1]
wa.temps.51.80 <- wa.tulane.temp.51.80$fit[,1]
wa.temps.86.15 <- wa.tulane.temp.86.15$fit[,1]
tulane.times <- tulane_pollen[2]$sample.meta$age
tulanetimetemps <- matrix(c(tulane.times, mat.temps.51.80 ,mat.temps.86.15, wapls.temps.51.80, wapls.temps.86.15, wa.temps.51.80, wa.temps.86.15), nrow = 191, ncol = 7)
tulane.data <- as.data.frame(tulanetimetemps)
tulane.data <- setNames(tulane.data, c("Time", "MAT5180", "MAT 86-15", "WAPLS 51-80", "WAPLS 86-15", "WA5180", "WA 86-15"))

#Crystal Scatter Plots 

library(ggplot2)
#MAT vs WA 1951-1980 Temps
ggplot(tulane.data, aes(x=mat.temps.51.80, y=wa.temps.51.80)) + geom_point() + 
  geom_smooth(method = lm) + 
  labs(title="MAT vs WA 1951-1980",
       x="MAT Temps", y = "WA Temps")

#MAT vs WAPLS 1951-1980 Temps
ggplot(tulane.data, aes(x=mat.temps.51.80, y=wapls.temps.51.80)) + geom_point() + 
  geom_smooth(method = lm) + 
  labs(title="MAT vs WAPLS 1951-1980",
       x="MAT Temps", y = "WAPLS Temps")

#WA vs WAPLS 1951-1980 Temps
ggplot(tulane.data, aes(x=wapls.temps.51.80, y=wa.temps.51.80)) + geom_point() + 
  geom_smooth(method = lm) + 
  labs(title="WAPLS vs WA 1951-1980",
       x="WAPLS Temps", y = "WA Temps")

#MAT vs WA 1986-2015 Temps
ggplot(tulane.data, aes(x=mat.temps.86.15, y=wa.temps.86.15)) + geom_point() + 
  geom_smooth(method = lm) + 
  labs(title="MAT vs WA 1986-2015",
       x="MAT Temps", y = "WA Temps")

#MAT vs WAPLS 1986-2015 Temps
ggplot(tulane.data, aes(x=mat.temps.86.15, y=wapls.temps.86.15)) + geom_point() + 
  geom_smooth(method = lm) + 
  labs(title="MAT vs WAPLS 1986-2015",
       x="MAT Temps", y = "WAPLS Temps")

#WA vs WAPLS 1986-2015
ggplot(tulane.data, aes(x=wapls.temps.86.15, y=wa.temps.86.15)) + geom_point() + 
  geom_smooth(method = lm) + 
  labs(title="WAPLS vs WA 1986-2015",
       x="WAPLS Temps", y = "WA Temps")

#WA 1951-1980 vs WA 1986-2015 
ggplot(tulane.data, aes(x=wa.temps.51.80, y=wa.temps.86.15)) + geom_point() + 
  geom_smooth(method = lm) + 
  labs(title="WA 51-80 vs WA 86-15",
       x="WA 51.80 Temps", y = "WA 86.15 Temps")

#WAPLS 1951-1980 vs WAPLS 1986-2015
ggplot(tulane.data, aes(x=wapls.temps.51.80, y=wapls.temps.86.15)) + geom_point() + 
  geom_smooth(method = lm) + 
  labs(title="WAPLS 51-80 vs WAPLS 86-15",
       x="WAPLS 51.80 Temps", y = "WAPLS 86.15 Temps")

#MAT 1951-1980 vs MAT 1986-2015 
ggplot(tulane.data, aes(x=mat.temps.51.80, y=mat.temps.86.15)) + geom_point() + 
  geom_smooth(method = lm) + 
  labs(title="MAT 51-80 vs MAT 86-15",
       x="MAT 51.80 Temps", y = "MAT 86.15 Temps")


