library(neotoma)
library(dplyr)
library(purrr)


ssamp <- get_dataset(datasettype = "pollen surface sample",
                     gpid = c("Canada", "United States"))
sspub <- get_publication(ssamp)

assertthat::assert_that(length(ssamp) == length(sspub), msg = "There are missing publication objects.")

dsids <- (1:length(ssamp)) %>% 
  map(function(x) {
    data.frame(dsid = ssamp[[x]]$dataset.meta$dataset.id,
               map(sspub[[x]], function(y) y$meta) %>% bind_rows()) }) %>% 
  bind_rows()

#create list for sites with no publications 
sites_no_pubs <- list()

#create list for sites with publications
sites_w_pubs <- list()

#create list for publications
pub_years  <- vector()

for (i in 1:3032){ 
  currentsite <- sspub[[i]]
  if (is.na(currentsite[[1]]$meta$id)){
    sites_no_pubs <- c(sites_no_pubs, currentsite)
  } 
  if(!is.na(currentsite[[1]]$meta$id)) {
    # Need to sort sites w/ 1 pub from sites w/ 1+ pub
    # add sites to sites_w_pubs that only have one pub
    # make new list e.g sites_multi_pubs
    # add sites w/ 1+ pubs to sites_multi_pubs
    # for sites_multi_pubs, need to keep oldest
    # to keep oldest, go into metadata and look at years
    # so, smaller# pub will be oldest
    # when finished sorting, get 2 lists:
      #1) sites with one pub & year 
      #2) sites w/ multi pubs & oldest pub, and
    if(length(currentsite) == 1){
      site_id <- as.numeric(names(sspub[i]))
      currentsite[[1]]$meta$dataset.id <- site_id
      sites_w_pubs <- c(sites_w_pubs, currentsite)
    }
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



    site_id <- as.numeric(names(sspub[i]))
    currentsite[[1]]$meta$dataset.id <- site_id
    sites_w_pubs <- c(sites_w_pubs, currentsite)

#added the site_id and site_info lines into loop to help connect it to Collection Data - dates