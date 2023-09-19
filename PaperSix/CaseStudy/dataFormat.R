#This script is used to format the data from Pantraps and FIT counts
#for analysis with NIMBLE

# Packages need to run this script
library(tidyverse)
library(lubridate)
library(dplyr)
library(pbapply)

#Reading the data
genus_counts <- read.csv("CaseStudy/FIT_counts_with_effort.csv")
pan_traps <- read.csv("CaseStudy/pan_traps.csv")
species_lookup <- read.csv("CaseStudy/species_lookup.csv")
gr_ref <- read.csv("CaseStudy/gr_ref.csv")
colnames(gr_ref)[1] <- "X1km_square"

#Formating the date into year, month and day
pan_traps$date <- dmy_hms(pan_traps$date)
genus_counts$date <- dmy_hms(genus_counts$date)

###joining the insect group data to the species data
#Rename the column names as taxon aggregated to species
names(pan_traps)[names(pan_traps) == "taxon_aggregated"] <- "species"

#Merge the species pantrap data to the species_lookup and match them to their taxon_group
update_pantraps <- merge(pan_traps, species_lookup, by= "species")

# Format the group taxon data (FSV counts)
update_genus <-genus_counts %>%
  filter(insect_group != "all_insects_total")

pantraps_surveys <- unique(
  update_pantraps[,2:3]) %>%
  mutate(pan_trap = "survey")

FIT_surveys <- unique(
  update_genus[,1:2]) %>%
  mutate(FIT_count = "survey")

all_surveys <- pantraps_surveys %>%
  full_join(FIT_surveys)

write.csv(all_surveys, "CaseStudy/all_surveys.csv",
          row.names = FALSE)


#######################
# attempt at formatting
#######################
update_pantraps_wide <- update_pantraps %>%
  dplyr::select(-c("taxon_group", "n_traps", "insect_group")) %>%
  pivot_wider(id_cols = c("X1km_square", "date"),
              names_from = species, values_from = n_traps_present) %>%
  mutate_at(vars(3:169), ~replace_na(., 0))


update_genus_wide <- update_genus %>%
  pivot_wider(id_cols = c("X1km_square", "date"),
              names_from = insect_group, values_from = count) %>%
  mutate_at(vars(3:12), ~replace_na(., 0))

update_genus_nsurveys_wide <- update_genus %>%
  pivot_wider(id_cols = c("X1km_square", "date"),
              names_from = insect_group, values_from = n_surveys) %>%
  mutate_at(vars(3:12), ~replace_na(., 0))
colnames(update_genus_nsurveys_wide)[3:12] <- c("beetlesNsurvey",
                                                "bumblebeesNsurvey",
                                                "butterflies_mothsNsurvey",
                                                "honeybeesNsurvey",
                                                "hoverfliesNsurvey",
                                                "insect_otherNsurvey",
                                                "insect_smallNsurvey",
                                                "other_fliesNsurvey",
                                                "solitary_beesNsurvey",
                                                "waspsNsurvey")

all_data_wide <- full_join(update_pantraps_wide,
                           update_genus_wide,
                           by = c("X1km_square", "date"))%>%
  full_join(., update_genus_nsurveys_wide,
            by = c("X1km_square", "date"))

#Number of sites
nsites <- length(unique(all_data_wide$X1km_square))

#unique sites for the study
sites <- data.frame(X1km_square = sort(unique(all_data_wide$X1km_square)))

data_gen <- function(group_name){
  if(!is.character(group_name)) stop("Group name must be a character")

  group_data <- all_data_wide %>%
    dplyr::select(c(
      "X1km_square",
      "date",
      group_name, #Must be a character
      paste0(group_name, "Nsurvey"),
      which(names(all_data_wide)%in%species_lookup[species_lookup$insect_group==group_name, "species"] == TRUE)))%>%
    dplyr::mutate(year = lubridate::year(date),
                  month = lubridate::month(date),
                  day = lubridate::day(date))

  species_names = colnames(group_data)[5:(ncol(group_data)-3)]

  years = sort(unique(group_data$year))
  nyears = length(years)

  #number of species
  nspecies = ncol(group_data) - 7

  #number of sampling visits in a year
  nvisits= 4

  species_data <- array(NA, dim = c(nsites, nspecies, nvisits, nyears))
  genus_data <- array(NA, dim = c(nsites, nvisits, nyears))
  genus_data_Nsurvey <- array(1, dim = c(nsites, nvisits, nyears))


  for(year.tag in 1:nyears){
    group_data_new <-  group_data%>%
      dplyr::filter(year == years[year.tag])%>%
      dplyr::arrange(c("date"))%>%
      dplyr::group_by(X1km_square)%>%
      mutate(rank=as.numeric(as.factor((date))))

    #number of visits
    nvisits = length(unique(group_data_new$rank))

    for(visit.tag in 1:nvisits){
      visit_group_data = group_data_new %>%
        filter(rank == visit.tag)%>%
        full_join(sites, by="X1km_square")%>%
        arrange(X1km_square)%>%
        ungroup()%>%
        dplyr::select(-c("X1km_square", "date", "year", "month", "day", "rank"))

      species_data[ , , visit.tag, year.tag] <- visit_group_data[ , 3:ncol(visit_group_data)]%>%
        as.matrix()

      genus_data[ ,visit.tag, year.tag] =  visit_group_data[ , 1]%>%
        as.matrix()

      genus_data_Nsurvey[ ,visit.tag, year.tag] =  visit_group_data[ , 2]%>%
        as.matrix()

    }
  }

  # # ####
  #  speciesData <- rbind(species_data[ , , 1, 1],
  #                      species_data[ , , 2, 1],
  #                       species_data[ , , 3, 1],
  #                       species_data[ , , 4, 1],
  #                       species_data[ , , 1, 2],
  #                       species_data[ , , 2, 2],
  #                       species_data[ , , 3, 2],
  #                       species_data[ , , 4, 2])
  # #
  #  genusData <- c(genus_data[ ,1, 1],
  #                     genus_data[ ,2, 1],
  #                     genus_data[ ,3, 1],
  #                     genus_data[ ,4, 1],
  #                    genus_data[ ,1, 2],
  #                     genus_data[ ,2, 2],
  #                     genus_data[ ,3, 2],
  #                     genus_data[ ,4, 2])
  # #
  #  Nsurveys <- c(genus_data_Nsurvey[ ,1, 1],
  #                genus_data_Nsurvey[ ,2, 1],
  #                genus_data_Nsurvey[ ,3, 1],
  #                genus_data_Nsurvey[ ,4, 1],
  #                genus_data_Nsurvey[ ,1, 2],
  #                genus_data_Nsurvey[ ,2, 2],
  #                genus_data_Nsurvey[ ,3, 2],
  #                genus_data_Nsurvey[ ,4, 2])
  #  Nsurveys[is.na(Nsurveys)] <- 1
  genus_data_Nsurvey[is.na(genus_data_Nsurvey)] <- 1

  Nsurveys = cbind(genus_data_Nsurvey[,,1], genus_data_Nsurvey[,,2])
  genusData = cbind(genus_data[,,1], genus_data[,,2])
  speciesData <- array(NA, dim = c(nsites, nspecies, nvisits*2))
  for(i in 1:8){
    if(i < 5){
      speciesData[,,i] <-  species_data[,,i,1]
    }else{
      speciesData[,,i] <-  species_data[,,(i-4),2]
    }
  }
  # speciesData = cbind(species_data[,,1,1],
  #                     species_data[,,2,1],
  #                     species_data[,,3,1],
  #                     species_data[,,4,1],
  #                     species_data[,,1,2],
  #                     species_data[,,2,2],
  #                     species_data[,,3,2],
  #                     species_data[,,4,2])

latitude = sites%>%
  left_join(gr_ref, by = "X1km_square")%>%
  dplyr::select("northing")

visit <- rep(1:8, each = nsites)

n.replicates = rep(5, 74)

  # used for CV in nimble
  indx <- sample(1:74, 74, replace = FALSE)

  return(list(mat.species=speciesData,
              mat.genus =  genusData,
              Nsurveys = Nsurveys,
              latitude = latitude$northing,
              visit = visit,
              sites = sites,
              species_names = species_names,
              indx = indx))
}

#The years and group name
#years = list(2017, 2018)
group= list("bumblebees",  "hoverflies", "solitary_bees")

# Getting the data
simulations_all <- pblapply(group, function(x){
  #pblapply(years, function(z){
  data <- data_gen(group_name=x)
  # }, cl=4)
}, cl=4)

#save the data
save(simulations_all, file="CaseStudy/data_idm.RData")

