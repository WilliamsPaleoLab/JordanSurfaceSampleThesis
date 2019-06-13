#if use simons code, look in ssamp: submission data and also dataset.meta has some 
#collection data (if there is any)

#to sort sites with collection dates and sites without, but doesn't work; returns error of argument length 0
#for (i in 1:3032){
  #currentsite <- ssamp[[i]]
  #if(is.na(currentsite[[1]]$meta.ageold)) {
    #sites_no_coldat <- c(sites_no_coldat, currentsite)
  #}
  #if(!is.na(currentsite[[1]]$meta.ageold)) {
   # sites_w_coldat <- c(sites_w_coldat, currentsite)
#  }
#}

#i think works when don't have gpid in there, so would just need to filter out to get us and canda somehow?
#sites_w_coldat <- neotoma::get_dataset(datasettype = "pollen surface sample",ageold = 0, ageyoung = -60)
#also trying to get coll date; below lets me see age young and old when i click on it
#dwndld <- get_download(sites_w_coldat, verbose = TRUE)


#ACTUAL CODE TO USE/THAT WORKS FOLLOWS BELOW
#Getting Collection date 
CollectionData <- neotoma::get_table('CollectionUnits')

#David - loop below; loop works to add collection date to sites w pubs list
id_index <- grep("id", sites_w_pubs)
for(i in id_index){
  site_id1 <- sites_w_pubs[[i]]$meta$id
  site_collection_date <- CollectionData[site_id1 == CollectionData$SiteID, 9]
  site_data_na_false <- !is.na(site_collection_date)
  final_collect_date <- site_collection_date[site_data_na_false]
  if (is_empty(final_collect_date)) {
    final_collect_date = NA
  }
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
  #sort smallest using other sort example in simonfixbadcode file
