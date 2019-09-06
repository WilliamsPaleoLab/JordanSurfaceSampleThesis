library(neotoma)
library(dplyr)
library(purrr)


ssamp <- get_dataset(datasettype = "pollen surface sample",
                     gpid = c("Canada", "United States")) #Gets all pollen surface sample datasets in Canada & US
sspub <- get_publication(ssamp) #Gets publications for the above datasets 

assertthat::assert_that(length(ssamp) == length(sspub), msg = "There are missing publication objects.")

dsids <- (1:length(ssamp)) %>% 
  map(function(x) {
    data.frame(dsid = ssamp[[x]]$dataset.meta$dataset.id,
               map(sspub[[x]], function(y) y$meta) %>% bind_rows()) }) %>% 
  bind_rows()

#Creates list for sites with no publications 
sites_no_pubs <- list()

#Creates list for sites with publications
sites_w_pubs <- list()

#Creates vector for publications
pub_years  <- vector()

#Loop to sort sites from ones w pubs to ones w no pubs
for (i in 1:3032){ 
  currentsite <- sspub[[i]]
  if (is.na(currentsite[[1]]$meta$id)){
    #^^ checks if current site has no publication (by going into metadata), if no pubs, put into list
    sites_no_pubs <- c(sites_no_pubs, currentsite)
  } 
  if(!is.na(currentsite[[1]]$meta$id)) {
    # If site does have publication (!is.na - if it does), follows loop below 
    # Checks if site has only 1 pub, then put into list - sites_w_pubs
    if(length(currentsite) == 1){
      site_id <- as.numeric(names(sspub[i]))
      currentsite[[1]]$meta$dataset.id <- site_id
      sites_w_pubs <- c(sites_w_pubs, currentsite)
    }
    # If site has 2 or more pubs, loop checks which publication is older 
    # By looking at metadata, then, if checks for the oldest pub to keep 
    # And put into sites_w_pubs list so each site has the oldest pub
    if(length(currentsite) >= 2) {
      for(P in length(currentsite)){
        currentpub_year <- currentsite[[P]]$meta$year
        pub_years <- c(pub_years,currentpub_year)
      }
      sort(pub_years, decreasing = FALSE)
      earliest_pub_year <- pub_years[1]
      pub_years <- vector()
      for(P in length(currentsite)){ 
        currentpub <- list(currentsite[[P]])
        if (currentsite[[P]]$meta$year == earliest_pub_year){
          site_id <- as.numeric(names(sspub[i]))
          currentpub[[1]]$meta$dataset.id <- site_id
          sites_w_pubs <- c(sites_w_pubs, currentpub)
        }
      }
    }
  }
}

for(N in 3032) {
  sitename <- ssamp[[N]]
  if (!is.na(sitename[[1]]$'site.data'$'site.name'))
    ss.site.names <- c(ss.site.names, sitename)
 }

#writes ssamp printed to a txt file
options(max.print = 10000)
sink("output.txt")
print(ssamp)
sink()


