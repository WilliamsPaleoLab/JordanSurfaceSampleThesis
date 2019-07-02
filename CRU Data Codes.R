library(ncdf4)
library(fields)
library(ncdf4.helpers)

#Open cru data - with code help from kevin
cru_data_all <- ncdf4::nc_open
cru_data_all <- ncdf4::nc_open("Documents/cru_ts4.02.1901.2017.tmp.dat.nc")

#Get variables from open data
#Temp is only variable in list but can get lat/long as well

#Getting longitude for north america 
cru_long <- ncvar_get(cru_data_all, 'lon', start = 400, count = 111)
long_df <- data.frame(cru_long)
# NorthAm_lon <- long_df[100:239, 1]
# NorthAm_lon <- data.frame(NorthAm_lon)

#Getting latitude for north america
cru_lat <- ncvar_get(cru_data_all, 'lat', start = 11, count = 77)
lat_df <- data.frame(cru_lat)
# NorthAm_lat <- lat_df[210:349,1]
# NorthAm_lat <- data.frame(NorthAm_lat)

#Getting time from file
time <- ncvar_get(cru_data_all, 'time')

#Getting temp for north america ** 
cru_temp <- ncvar_get(cru_data_all, varid = 'tmp') #, start = NA, count = NA, verbose = FALSE)
extracted_crutemp <- subset(cru_temp, select = c(50:100, 50:100, 50:100))
cru_temp_df <- data.frame(extracted_crutemp)
temp
nc_close(cru_data_all)

temp_df[400:600,]
dim(cru_temp)

#Whole world map
fields::image.plot(cru_temp[,,1])
#Plot for North America only
fields::image.plot(cru_temp[77,111,1])

# # For analysis below, a year starts with december of the previous year to satisfy:
# # winter = djf, sprg = mam, summer = jja, fall = son
# season_list = c('DJF','MAM','JJA', 'SON') 
# season_num = 4 # for annual, season_num = 1; for monthly, season_num = 12
# season =  matrix(c(12,1,2,3,4,5,6,7,8,9,10,11), c(12/season_num, season_num))
# 
# # variable lists
# cru_var_list <- c('tmp') # variable names, we could add 'prcp' if we needed the precip data. 
# cru_start <- c(1, 1, 720) # start at December 1960
# cru_count <- c(-1, -1, 360) # include all indexes of lat and lon, and then 30 years of data (12*30=360).
# 
# # create a climatology for each variable
# create_climatology <- function() {
#   f_cru <- ncdf4::nc_open("Documents/cru_ts4.02.1901.2017.tmp.dat.nc")
#   cru_clim <- ncvar_get(f_cru, cru_var_name, start = cru_start, count = cru_count)
#   cru_lon <- ncvar_get(f_cru, 'lon')
#   cru_lat <- ncvar_get(f_cru, 'lat')
#   dim(cru_clim) <- c(dim(cru_clim)[1], dim(cru_clim)[2], 12, dim(cru_clim)[3] / 12) #reshape our data to lon:lat:months:years
#   for (i in 1:season_num) {
#     months <- season[, i]
#     
#     # pick out the corresponding months. Subset cru data to only DJF, MAM etc.
#     clim_seasonal <- apply(cru_clim[, , months, ], c(1, 2, 4), mean) # this takes the average across months, so we have 30 annual means
#     
#     # calculate the 1960-1990 climatologies. (Now take mean of those 30 years for given season) 
#     # clim_seasonal_mean <- apply(clim_seasonal, c(1,2), mean) # mean
#     # clim_seasonal_sd <- apply(clim_seasonal, c(1,2), sd) # standard deviation
#     # 
#     # assign(paste0('CRU_1960_1990_', cru_var_name, '_season_', i, '_mean'), clim_seasonal_mean) # seasonal variables for debiasing
#     # assign(paste0('CRU_1960_1990_', cru_var_name, '_season_', i, '_sd'), clim_seasonal_sd) # seasonal variables for debiasing
#   } # end for seaason_num
#   nc_close(f_cru)
# }
# 
# create_climatology()
# 
# #Separate lat and lon for North America from f_cru
# for(i in length(f_cru$dim$lon$vals))
