
#i think works when don't have gpid in there, so would just need to filter out to get us and canda somehow?
#sites_w_coldat <- neotoma::get_dataset(datasettype = "pollen surface sample",ageold = 0, ageyoung = -60)
#also trying to get coll date; below lets me see age young and old when i click on it
#dwndld <- get_download(sites_w_coldat, verbose = TRUE)
#######################

# Creates sorted collection sites list
sorted_collection_sites <- sites_w_pubs

sorted_sites_NA <- list()

#Getting Collection date 
CollectionData <- neotoma::get_table('CollectionUnits')

# David - loop below; loop works to add collection date to list of sites_w_pubs
id_index <- grep("id", sorted_collection_sites) # ID_index connects collection data table to sites_w_pubs list
for(i in id_index){
  site_id1 <- sorted_collection_sites[[i]]$meta$id
  site_collection_date <- CollectionData[site_id1 == CollectionData$SiteID, 9] #Adds collect date to list, gets info from collection date in table, column 9 
  site_data_na_false <- !is.na(site_collection_date) # If there is col. date, assigns sites 
  final_collect_date <- site_collection_date[site_data_na_false]
  if (is_empty(final_collect_date)) {
    final_collect_date = NA
  } #If - checks to see if site has actual collection date; if it doesn't then assigned value NA
  if (length(final_collect_date) >= 2){
    # If site has more than one collection date, sorts them from oldest to youngest 
    # Keeps oldest as the value 
    sort(final_collect_date, decreasing = FALSE)
    final_collect_date <- final_collect_date[[1]]
  }
  sorted_collection_sites[[i]]$meta$collection_date <- final_collect_date
}

#if (length(currentsite) >= 2) {
    #for (D in length(currentsite))
     #collection_dates <- c(sites_w_pubs, collection_dates) 
  #}
  #sort(collection_dates, decreasing = FALSE)
  #for(D in length(currentsite)){
   #collection_dates <- list(currentsite[[D]])
    #if (currentsite[[D]]$meta$collection_date == collection_date){
     #collection_date[[1]]$meta$collection_date <- final_collect_date
      #final_collect_date <- c(sites_w_pubs, final_collect_date)
    #}
 #}
#}

