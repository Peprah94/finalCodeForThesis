options(ggplot2.continuous.colour="viridis")
options(ggplot2.continuous.fill = "viridis")

#Functions needed to run this script
source('spde-book-functions.R')

#packages
library(sp)
library(sf)
library(rgeos)
library(INLA)
library(dplyr)
library(raster)
library(pbapply)
library(reshape)
library(tiff)
library(ggplot2)
library(gridExtra)
library(readr)
library(terra)
library(tidyr)
library(stringr)
library(inlabru)

#BNGproj <-  CRS(SRS_string = "EPSG:27700")
BNGproj <- CRS("+proj=robin +datum=WGS84 +units=km")
#BNGproj <- CRS("+proj=utm +zone=33 +datum=WGS84 +units=km")
#proj=utm +zone=33 +datum=WGS84
#crsModel <- "+proj=robin +datum=WGS84 +units=km"
#> Warning in showSRID(SRS_string, format = "PROJ", multiline = "NO", prefer_proj =
#> prefer_proj): Discarded datum OSGB 1936 in Proj4 definition
wkt_BNGproj <- wkt(BNGproj)
cat(wkt_BNGproj)

# Set the number of species and the names of species interested
n.species = nspecies = 2
speciesName = c("Larus argentatus", "Larus canus")

## Initial data manipulation ##
datagbif <- read_delim("~/Documents/GitHub/Sampling-process-in-CS/gullData/allGullsDataset/occurrence.txt",
                       delim = "\t",
                       escape_double = FALSE,
                       trim_ws = TRUE)%>%
  dplyr::filter(countryCode == "NO")

#set the species name
datagbif <- datagbif[which(!is.na(datagbif$decimalLatitude)),]%>%
  mutate(reported = species)

event <- read_delim("dwca-tove_birdsampling/event.txt",
                    delim = "\t", escape_double = FALSE,
                    trim_ws = TRUE)

#Add states
states <- geodata::gadm(country = "NOR", level = 2, path = tempdir())

# Define the data as spatialPointsDataFrame and set the correct projection
datagbifDF <- sp::SpatialPointsDataFrame(cbind(as.numeric(datagbif$decimalLongitude),
                                               as.numeric(datagbif$decimalLatitude)),
                                         datagbif,
                                         proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))%>%
  sf::st_as_sf()

datagbifNew <- terra::vect(datagbifDF)
datagbifNew <- terra::project(datagbifNew,
                              BNGproj)

# set new cordinate referenec system
#newcrs <- CRS("+proj=longlat +datum=WGS84 +units=km")
#BNGproj <-  sp::CRS(SRS_string = "+proj=longlat +datum=WGS84 +units=km")
#crsModel <- "+proj=longlat +datum=WGS84 +units=km"
#BNGproj<-  CRS(SRS_string = "EPSG:27700")
#> Warning in showSRID(SRS_string, format = "PROJ", multiline = "NO", prefer_proj =
#> prefer_proj): Discarded datum OSGB 1936 in Proj4 definition
# wkt_BNGproj <- wkt(BNGproj)
# cat(wkt_BNGproj)
# # Extract administrative regions of Norway
# norwaybu <-raster::getData("GADM",
#                           country="NOR",
#                           level=1)%>%
#  sf::st_as_sf()%>%   terra::vect()
#
# #transform the regions to match the new crs
# norwaybu <- terra::project(norwaybu,
#                          BNGproj)


#extract the polygon for Trondheim
#index <- which(norwaybu$NAME_1 == 'Nordland')
#index1 <- which(states$NAME_1 == "Sør-Trøndelag")
index1 <- which(states$NAME_1 == "Nord-Trøndelag")
#index1 <- which(states$NAME_1 == "Nord-Trøndelag")
#states$NAME_2

#index <- which(norwaybu$NAME_1 == "Sør-Trøndelag")
#index1 <- which(norwaybu$NAME_1 == "Nord-Trøndelag")
# norwaybu <- states[c(index,index1),]%>%
#   sf::st_as_sf()%>%
#   terra::vect()

#norwaybu <- states[c(index,index1),]
norwaybu <- states[c(index1),]

norway <- terra::project(norwaybu,
                         BNGproj)

spNorway <- norway%>%
  as(., "Spatial") #convert from spatVector to spatialPolygons for use in INLA

# Data cleansing: selecting human obs only, formatting the year and
# formatting the data and selecting data collected after 2006
dataGbifTrondheim <- terra::intersect(datagbifNew, norway)
#select only human observation records


# extract count data
countData <- dataGbifTrondheim[dataGbifTrondheim$datasetName%in% c("Seabirds in Norway - Estimated population sizes") ,]
nbic <- dataGbifTrondheim[dataGbifTrondheim$institutionCode%in% c("miljolare") ,]

# Extract presence absence
dataGbifTrondheim <- dataGbifTrondheim[which(dataGbifTrondheim$basisOfRecord=="HUMAN_OBSERVATION"),]
dates <- paste(dataGbifTrondheim$day,dataGbifTrondheim$month,dataGbifTrondheim$year,sep="/")
dates <- as.Date(dates,format = "%d/%m/%Y")
dataGbifTrondheim <- dataGbifTrondheim[which(dates > "2006-01-01"),]
dataGbifTrondheim <- dataGbifTrondheim[!dataGbifTrondheim$institutionCode%in% c("miljolare") ,]





indxDuplicate <- which(!duplicated(cbind(dataGbifTrondheim$decimalLatitude, dataGbifTrondheim$decimalLongitude)))
dataGbifTrondheim <- dataGbifTrondheim[indxDuplicate,]

#get the coordinates
coordsGbifTrondheim <- dataGbifTrondheim[ ,c("decimalLongitude","decimalLatitude")]
  


# Transform coordinates to spatial points dataframe
# spPointsCoordsTrondheim <- SpatialPoints(coordsGbifTrondheim,
#                                          proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))%>%
#   sf::st_as_sf()%>%
#   terra::vect()

spPointsCoordsTrondheim <- coordsGbifTrondheim
#project to the new coordinate system
spPointsCoordsTrondheim <- terra::project(spPointsCoordsTrondheim,
                                          BNGproj)

#ensure the points are within the boundary of Trondheim
spPointsGbif <- terra::intersect(norway,
                                 spPointsCoordsTrondheim)%>%
  as(., "Spatial" ) # INLA only used sp class, so need to comvert sf and terra class to sp


norwaySP <- norway%>%
  as(., "Spatial")%>%
  spTransform(., CRS=BNGproj)

norwaySt <- norwaySP%>%
  st_as_sf()

spPointsGbifSf <- spPointsGbif%>%
  st_as_sf()

spPointsGbifEdit <- gIntersection(norwaySP,spPointsGbif)

#dataGbifTrondheim <- terra::intersect(dataGbifTrondheim, vect(spPointsGbifEdit))
#dataGbifTrondheimEdit <- st_intersection(, norwaySt)
  
  #terra::intersect(norway, 
   #                                       dataGbifTrondheim)
# # Presence-Absence data
# occurrenceNINA <- read_delim("~/Documents/GitHub/Sampling-process-in-CS/gullData/presenceAbsence/occurrence.txt", 
#                              delim = "\t", escape_double = FALSE, 
#                              trim_ws = TRUE)%>%
#   dplyr::filter(countryCode == "NO")%>%
#   dplyr::filter(institutionCode %in% c("NINA"))
# 
# occurrenceNTNU <- read_delim("~/Documents/GitHub/Sampling-process-in-CS/gullData/presenceAbsence/occurrence.txt", 
#                              delim = "\t", escape_double = FALSE, 
#                              trim_ws = TRUE)%>%
#   dplyr::filter(countryCode == "NO")%>%
#   dplyr::filter(institutionCode %in% c("NTNU-VM"))
# 
# occurrencenof <- read_delim("~/Documents/GitHub/Sampling-process-in-CS/gullData/presenceAbsence/occurrence.txt", 
#                             delim = "\t", escape_double = FALSE, 
#                             trim_ws = TRUE)%>%
#   dplyr::filter(countryCode == "NO")%>%
#   dplyr::filter(institutionCode %in% c("nof"))
# Defining the MESH
#https://rpubs.com/jafet089/886687
premesh <- inla.sp2segment(norwaySP)

#BNGproj<-  CRS(SRS_string = "EPSG:27700")

#Calculate maximum edge
max.edge = diff(range(spPointsGbifEdit@coords[,1])/(3*5))
bound.outer <- diff(range(spPointsGbifEdit@coords[,1])/(5))

proj4string = CRS(projection(norwaySP))
mesh <- fm_mesh_2d(boundary = norwaySP,
                     loc = spPointsGbifEdit,
                     max.edge  = c(1,3)*max.edge,
                     offset = c(max.edge,bound.outer),
                     cutoff = 5,
                     #min.angle = 60,
                     crs = BNGproj
)

#df <- pixels(mesh, mask = norwaySP )
#identicalCRS(mesh,boundary)
plot(mesh)

# Format the TOV-E Sampling data


# NOTE: Currebly INLA only uses spatial class:sp.
#All other classes such as spatVec and sf are not supported



##################
# Formatting detection data
###################
# Use time spent as a covariate
#load data
eventDataTOVE <- read_delim("dwca-tove_birdsampling/event.txt",
                            delim = "\t", escape_double = FALSE,
                            trim_ws = TRUE)

occurrence <- read_delim("dwca-tove_birdsampling/occurrence.txt",
                         delim = "\t", escape_double = FALSE,
                         trim_ws = TRUE)

mergeDF <- merge(occurrence,eventDataTOVE[, c("samplingEffort","id" )], by = "id", all.x = TRUE, all.y = TRUE )

toveSamplingGbif <- read_delim("tove_sampling.csv",
                               delim = "\t", escape_double = FALSE,
                               trim_ws = TRUE)%>%
  merge(.,
        mergeDF[, c("catalogNumber", "samplingEffort")],
        by = "catalogNumber",
        all.x = TRUE,
        all.y = TRUE)%>%
  #filter(stateProvince == "Trøndelag")%>%
  filter(species %in% c("Larus argentatus", "Larus canus"))%>%
  arrange(year)%>%
  group_by(locality)%>%
  mutate(visit = as.numeric(as.factor(eventDate)),
         presence = 1)%>%
  ungroup()%>%
  separate(samplingEffort, c("minutes","samplingUnits"), "/")%>%
  mutate(minutes = as.numeric(unlist(str_extract_all(minutes, "[0-9]+"))),
         samplingUnits = as.numeric(unlist(str_extract_all(samplingUnits, "[0-9]+"))))%>%
  mutate(timeSpent = minutes/samplingUnits)#%>%
  #mutate(timeSpent = scale(timeSpent))

# Select the data in the given region
toveSamplingGbifSf <- st_as_sf(toveSamplingGbif, 
                               coords = c("decimalLongitude",
                                          "decimalLatitude"),
                               crs = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))%>%
  st_transform(., crs = BNGproj)
  #as(., "Spatial")%


toveSamplingGbifTrond <- st_intersection(toveSamplingGbifSf, norwaySt)

 #gIntersection(norwaySP,toveSamplingGbifSf)
# Make the data in wider format
toveSamplingGbifWider <- toveSamplingGbifTrond%>%
  as.data.frame()%>%
  dplyr::select(c("visit", "species", "presence", "locality",  "timeSpent", "geometry"))%>%
  group_by(locality, visit, geometry)%>%
  mutate(index = cur_group_id())%>%
  ungroup()%>%
  #select(-c("visit"))%>%
  tidyr::pivot_wider(.,
                     names_from = species,
                     values_from = presence,
                     values_fill = list(presence = 0))%>%
  dplyr::mutate(nTrials = length(unique(eventDataTOVE$year)))

#Organising as detection probabilty
#coordsTOVE <- toveSamplingGbifWider[,c("decimalLongitude","decimalLatitude")]
detectData <- toveSamplingGbifWider[,c("locality","Larus argentatus","Larus canus", "timeSpent", "visit", "index", "geometry", "nTrials")]%>%
  dplyr::group_by(locality, geometry, nTrials)%>%
  dplyr::summarise_at(.vars = names(.)[2:4],
                      .funs = c(sum = "sum"))%>%
  ungroup()

#coordsTOVE <- detectData[,c("decimalLongitude","decimalLatitude")]

#nVisitsForEvents = max(c(max(detectData$`Larus argentatus_sum`), max(detectData$`Larus canus_sum`)))

detectionDF <- occurenceDF <- countDf <- nbicDf <- list()
toveSF <- detectData%>%
  sf::st_as_sf()


#Lets go back and format the count data
countDataSf <- countData

#remember we have count data from NBIC
nbicSf <- nbic%>%
  st_as_sf()%>%
 # st_set_crs(., CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))%>%
  st_transform(., BNGproj)
# %>%
#   as.data.frame()%>%
#   tidyr::pivot_wider(.,
#                      names_from = species,
#                      values_from = individualCount,
#                      values_fill = list(individualCount = 0))

for(i in 1:nspecies){
  #For detection
  # sp::SpatialPointsDataFrame(coordsTOVE,
  #                            data = data.frame(detectData[ ,i+3],
  #                                              rep(nVisitsForEvents, nrow(coordsTOVE)),
  #                                              detectData[ ,6]),
  #                            proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

  
  #occupied sites
  # toveSF1 <- sp::SpatialPointsDataFrame(coordsTOVE,
  #                                       data = data.frame(as.numeric(detectData[ ,i+3]>0),
  #                                                         rep(1, nrow(coordsTOVE)),
  #                                                         detectData[ ,6]),
  #                                       proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))%>%
  #   sf::st_as_sf()
  # 
  # #convert to terra
  # toveSFNew <- terra::vect(toveSF)
  # toveSFNew1 <- terra::vect(toveSF1)
  
  #extract detection for each species
  detectionDF[[i]] <- toveSF%>%
    dplyr::select(c(i+3, 3, 6))
  
  #extract occupancy history
  occurenceDF[[i]] <- detectionDF[[i]]
  if(i == 1){
  occurenceDF[[i]]$`Larus argentatus_sum` <- as.numeric(occurenceDF[[i]]$`Larus argentatus_sum` > 0)
  }else{
   occurenceDF[[i]]$`Larus canus_sum` <- as.numeric(occurenceDF[[i]]$`Larus canus_sum` > 0) 
  }
  
  names(detectionDF[[i]])[1:3] <- c(paste0("detdata",i),
                                    "Ntrials",
                                    paste0("timeSpent",i))
  names(occurenceDF[[i]])[1:3] <- c(paste0("detdata",i),
                                    "Ntrials",
                                    paste0("timeSpent",i))
    
  
  
  countDf[[i]] <- countData%>%
    st_as_sf()%>%
    filter(species %in% speciesName[i])%>%
    dplyr::select(c("individualCount", "geometry"))%>%
    dplyr::filter(individualCount < 1000)
  
  nbicDf[[i]] <- nbicSf%>%
    filter(species %in% speciesName[i])%>%
    dplyr::select(c("individualCount", "year", "decimalLatitude", "geometry"))
         #   terra::project(toveSFNew,
  #                                    "+proj=robin +datum=WGS84 +units=km")%>%
  #   sf::st_as_sf()%>%
  #   as(., "Spatial" )
  # 
  # detectionDFOcc[[i]] <- terra::project(toveSFNew1,
  #                                       "+proj=robin +datum=WGS84 +units=km")%>%
  #   sf::st_as_sf()%>%
  #   as(., "Spatial" )
  # names(detectionDF[[i]]) <- c(paste0("detdata",i),
  #                              "Ntrials",
  #                              paste0("timeSpent",i))
  # names(detectionDFOcc[[i]]) <- c(paste0("detdata",i),
  #                                 "Ntrials",
  #                                 paste0("timeSpent",i))
}
detdata = detectionDF


#Citizen science data
#classifationData <- dataGbifTrondheim[which(dataGbifTrondheim$institutionCode %in% c("iNaturalist", "CLO")),]


#table(classifationData$species, classifationData$reported)
####
# Let's add the ML Prediction scores now
####
#Add the ML data
observations <- read_csv("observations.csv")%>%
  filter(`Prediction 1` %in% c("Larus argentatus",
                               "Larus marinus",
                               "Larus canus", 
                               "Larus fuscus"))

rr <- cbind(observations$`Prediction 1 score`,
            observations$`Prediction 2 score`,
            observations$`Prediction 3 score`,
            observations$`Prediction 4 score`,
            observations$`Prediction 5 score`)%>%
  as.matrix()


data_obs <- data.frame(gbifID = rep(observations$gbifID, 5),
                       scientificName = rep(observations$scientificName, 5),
                       occurenceID = rep(observations$occurrenceID, 5),
                       predicted_species = c(observations$`Prediction 1`,
                                             observations$`Prediction 2`,
                                             observations$`Prediction 3`,
                                             observations$`Prediction 4`,
                                             observations$`Prediction 5`
                       ),
                       predicted_score = c(observations$`Prediction 1 score`,
                                           observations$`Prediction 2 score`,
                                           observations$`Prediction 3 score`,
                                           observations$`Prediction 4 score`,
                                           observations$`Prediction 5 score`))%>%
  tidyr::pivot_wider(., id_cols = c("gbifID", "occurenceID"),
                     names_from = predicted_species,
                     values_from = predicted_score,
                     values_fill = 0.000001
  ) %>%
  dplyr::rowwise() %>% 
  dplyr::mutate(others = sum(c_across(-c(1,2,3,5)), na.rm = T))%>%
  dplyr::ungroup()%>%
  #dplyr::select(1:5,8, 377)%>%
  dplyr::select(1,2,3,5, 320)%>%
  dplyr::arrange(gbifID)

colnames(data_obs) <- c("gbifID", "occurenceID", "herring",
                        "common",
                        "other")

#reorder column names
col_order <- c("gbifID", "occurenceID", "common",
               "herring",
               "other")
data_obs <- data_obs[, col_order]

classifationDataWithWeights <- as.data.frame(dataGbifTrondheim)%>%
  dplyr::left_join(data_obs, by = c("gbifID"))%>%
  mutate_at(vars(common, herring), ~replace_na(., 1))

classificationDataWithML <-sp::SpatialPointsDataFrame(cbind(as.numeric(classifationDataWithWeights$decimalLongitude),
                                                            as.numeric(classifationDataWithWeights$decimalLatitude)),
                                                      classifationDataWithWeights,
                                                      proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))%>%
  sf::st_as_sf()






dataBeforeClassification <- list()
for(i in 1:nspecies){
  message(paste0("Sorting detection dataframe for ", speciesName[i]))
  dataBeforeClassification[[i]] <- classificationDataWithML[which(classificationDataWithML$species== speciesName[i] ),]%>%
    sf::st_as_sf()%>%
    terra::vect(.)%>%
    terra::project(.,
                   "+proj=robin +datum=WGS84 +units=km")%>%
    sf::st_as_sf()%>%
    as(., "Spatial" )
  
  #slot(dataBeforeClassification[[i]], "proj4string") <-  BNGproj
  
  if(i == 1){
    dataBeforeClassification[[i]]$weight <- dataBeforeClassification[[i]]$common
  }else{
    dataBeforeClassification[[i]]$weight <- dataBeforeClassification[[i]]$herring 
  }
}

# classifationData <- as.data.frame(classifationData)
# # Track the misclassified data from the iNaturalist dataset
# classifationData$reported[classifationData$occurrenceID %in% c("https://www.inaturalist.org/observations/15915907",
#                                                                "https://www.inaturalist.org/observations/128864933",
#                                                                "https://www.inaturalist.org/observations/53400091",
#                                                                "https://www.inaturalist.org/observations/43040008",
#                                                                "https://www.inaturalist.org/observations/28570282",
#                                                                "https://www.inaturalist.org/observations/27668962",
#                                                                "https://www.inaturalist.org/observations/16598948",
#                                                                "https://www.inaturalist.org/observations/13974098",
#                                                                "https://www.inaturalist.org/observations/15915907",
#                                                                "https://www.inaturalist.org/observations/120659167",
#                                                                "https://www.inaturalist.org/observations/115519068",
#                                                                "https://www.inaturalist.org/observations/112944668",
#                                                                "https://www.inaturalist.org/observations/110009380",
#                                                                "https://www.inaturalist.org/observations/86631681",
#                                                                "https://www.inaturalist.org/observations/50428943",
#                                                                "https://www.inaturalist.org/observations/26390486",
#                                                                "https://www.inaturalist.org/observations/15792254",
#                                                                "https://www.inaturalist.org/observations/138227844",
#                                                                "https://www.inaturalist.org/observations/80586887",
#                                                                "https://www.inaturalist.org/observations/73305432",
#                                                                "https://www.inaturalist.org/observations/72270329",
#                                                                "https://www.inaturalist.org/observations/48892129",
#                                                                "https://www.inaturalist.org/observations/16926432",
#                                                                "https://www.inaturalist.org/observations/109471390",
#                                                                "https://www.inaturalist.org/observations/37614241",
#                                                                "https://www.inaturalist.org/observations/33940419",
#                                                                "https://www.inaturalist.org/observations/62904456",
#                                                                "https://www.inaturalist.org/observations/111344780",
#                                                                "https://www.inaturalist.org/observations/52498871"
# )] = "Larus argentatus"
# 
# 
# ##############################
# classifationData$reported[classifationData$occurrenceID %in% c("https://www.inaturalist.org/observations/133393346",
#                                                                "https://www.inaturalist.org/observations/7250296",
#                                                                "https://www.inaturalist.org/observations/54139803",
#                                                                "https://www.inaturalist.org/observations/52803570",
#                                                                "https://www.inaturalist.org/observations/46400645",
#                                                                "https://www.inaturalist.org/observations/19071786",
#                                                                "https://www.inaturalist.org/observations/110477951",
#                                                                "https://www.inaturalist.org/observations/32751987"
# )] = "other"
# 
# ########################
# 
# classifationData$reported[classifationData$occurrenceID %in% c("https://www.inaturalist.org/observations/120922696",
#                                                                "https://www.inaturalist.org/observations/92021997",
#                                                                "https://www.inaturalist.org/observations/19031769",
#                                                                "https://www.inaturalist.org/observations/134749540",
#                                                                "https://www.inaturalist.org/observations/101276215",
#                                                                "https://www.inaturalist.org/observations/4071081",
#                                                                "https://www.inaturalist.org/observations/116287400"
#                                                                
# )] = "other"
# 
# ########################
# 
# classifationData$reported[classifationData$occurrenceID %in% c("https://www.inaturalist.org/observations/72681262"
# )] = "Larus canus"
# 
# 
# ########################
# classifationData$reported[classifationData$occurrenceID %in% c("https://www.inaturalist.org/observations/58353519",
#                                                                "https://www.inaturalist.org/observations/54614542",
#                                                                "https://www.inaturalist.org/observations/25931117",
#                                                                "https://www.inaturalist.org/observations/34776968",
#                                                                "https://www.inaturalist.org/observations/29059722",
#                                                                "https://www.inaturalist.org/observations/17899929",
#                                                                "https://www.inaturalist.org/observations/7188731",
#                                                                "https://www.inaturalist.org/observations/67677363",
#                                                                "https://www.inaturalist.org/observations/56822719",
#                                                                "https://www.inaturalist.org/observations/40743121",
#                                                                "https://www.inaturalist.org/observations/128284445",
#                                                                "https://www.inaturalist.org/observations/96565757",
#                                                                "https://www.inaturalist.org/observations/120817745",
#                                                                "https://www.inaturalist.org/observations/49045975",
#                                                                "https://www.inaturalist.org/observations/47869589",
#                                                                "https://www.inaturalist.org/observations/113020448",
#                                                                "https://www.inaturalist.org/observations/121395057",
#                                                                "https://www.inaturalist.org/observations/87531960",
#                                                                "https://www.inaturalist.org/observations/119391382",
#                                                                "https://www.inaturalist.org/observations/88829645",
#                                                                "https://www.inaturalist.org/observations/45060713"
#                                                                
# )] = "other"

# Check the cross-tabulation of reported and verified data
#classifationDataFiltered <- classifationData[which(classifationData$scientificName %in% c("Larus argentatus Pontoppidan, 1763", "Larus canus Linnaeus, 1758")),]



#Get the extent of covariates needed
extentPlots <- extent(c(min(mesh$loc[,1]),
                        max(mesh$loc[,1]),
                        min(mesh$loc[,2]),
                        max(mesh$loc[,2])))

extentMesh <- extent(c(min(mesh$loc[,1])-2,
                       max(mesh$loc[,1])+2,
                       min(mesh$loc[,2])-2,
                       max(mesh$loc[,2])+2))

#extentSp <- spTransform(SpatialPoints(matrix(extentMesh, ncol=2, byrow=F),
#                                      proj4string = CRS(proj4string(spNorway)) ),
#                       CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

extentSt <- SpatialPoints(matrix(extentMesh, ncol=2, byrow=F),
                          proj4string = CRS(proj4string(spNorway)) )%>%
  st_as_sf()

extentSp <- extentSt%>%
  as(., "Spatial")

extentSpPoints <- extent(c(min(extentSp@coords[,1]),
                           max(extentSp@coords[,1]),
                           min(extentSp@coords[,2]),
                           max(extentSp@coords[,2])))

extSpPoints <- terra::ext(c(min(extentSp@coords[,1]),
                            max(extentSp@coords[,1]),
                            min(extentSp@coords[,2]),
                            max(extentSp@coords[,2]))) #created from terra package
extentSpPoly <- vect(extentSp, "polygons")

# Loading distance to road data
distanceTrondheim <- terra::rast("otherInfo/distanceRaster.tif")
#r_up_0.1 <- terra::disagg(distanceTrondheim, fact = c(3,2))
#change the raster dimension of distance
#a1 <- aggregate(altitudeTrondheim, 2)

cropDistanceTrondheim <- terra::crop(distanceTrondheim,
                                     extentSpPoly)
#cropDistanceTrondheim1 <- terra::mask(r_up_0.1,
#                                     norway)

plot(cropDistanceTrondheim,
     col=viridis(64, direction=-1),
     ext = extentPlots)
plot(mesh, add = TRUE)
points(spPointsGbifEdit, col = "white", pch = 19,cex=0.2)


# Loading distance to waterbody
distanceTrondheimWaterBody <- terra::rast("otherInfo/distanceRasterWaterBody.tif")
#r_up_0.1 <- terra::disagg(distanceTrondheim, fact = c(3,2))
#change the raster dimension of distance
#a1 <- aggregate(altitudeTrondheim, 2)

cropDistanceTrondheimWaterBody <- terra::crop(distanceTrondheimWaterBody,
                                              extentSpPoly)
#cropDistanceTrondheim1 <- terra::mask(r_up_0.1,
#                                     norway)

plot(cropDistanceTrondheimWaterBody,
     col=viridis(64, direction=-1),
     ext = extentPlots)
plot(mesh, add = TRUE)
points(spPointsGbif, col = "white", pch = 19,cex=0.2)
# loading elevation data
alt <- raster::getData('worldclim', var='alt', res=2.5)
#alt <- terra::rast("otherinfo/elevation.tif")
alt <- terra::rast(alt)
altitudeTrondheim <- terra::project(alt,
                                    "+proj=robin +datum=WGS84 +units=km",
                                    method = "near")
#altitudeTrondheim <- terra::disagg(altitudeTrondheim, fact = c(10,10))
#change the raster dimension of altitude
#a2 <- aggregate(altitudeTrondheim, 2)
#newNorwayRaster <- terra::rasterize(norway)
cropAltitudeTrondheim <- terra::crop(altitudeTrondheim,
                                     extentSpPoly
                                     )

plot(cropAltitudeTrondheim,
     ext = extentPlots)
plot(mesh, add = TRUE)
points(spPointsGbif, col = "white", pch = 19,cex=0.2)

#cropAltitudeTrondheimFilled <- focal(cropAltitudeTrondheim, w=3, na.policy = "only", na.rm = t)

#rad_hed1 <- mask(alt,norway)
#plot(rad_hed1,col=viridis(256))

# loading annual precipitation data
#bio <- raster::getData('worldclim', var='prec', res=2.5)

for(i in 1:12){
  if(i<10){assign(paste("precMonth",i,sep=""),terra::rast(paste("Precipitation25/prec0",i,".tif",sep="")))}
  else{assign(paste("precMonth",i,sep=""),terra::rast(paste("Precipitation25/prec",i,".tif",sep="")))}
  #assign(paste("precMonth",i,sep=""),terra::crop(get(paste("precMonth",i,sep="")),extSpPoints))
  assign(paste("precMonth",i,sep=""),terra::project(get(paste("precMonth",i,sep="")),"+proj=robin +datum=WGS84 +units=km",method = "near"))
}

#calculating average precipitation
averagePrecipitation  <- mean(precMonth1,
                              precMonth2,
                              precMonth3,
                              precMonth4,
                              precMonth5,
                              precMonth6,
                              precMonth7,
                              precMonth8,
                              precMonth9,
                              precMonth10,
                              precMonth11,
                              precMonth12)

#a1 <- aggregate(averagePrecipitation , 4)

cropPrecipitationTrondheim <- terra::crop(averagePrecipitation,
                                          extentSpPoly)

plot(cropPrecipitationTrondheim,
     ext = extentPlots)
plot(mesh, add = TRUE)
points(spPointsGbif, col = "white", pch = 19,cex=0.1)

#covariates definition
#cov1.spix <- as(raster(averagePrecipitation),"SpatialPixelsDataFrame")
#cov2.spix <- as(raster(distanceTrondheim),"SpatialPixelsDataFrame")
#cropDistanceTrondheim
#proj4string(cov2.spix) <- newcrs
#cov3.spix <- as(raster(altitudeTrondheim),"SpatialPixelsDataFrame")
#covslist <- list(cov1.spix,cov2.spix,cov3.spix)
#covs = covslist
cov1.spix <- raster(cropPrecipitationTrondheim)
cov1.spix <- as(raster(cropPrecipitationTrondheim),"SpatialPixelsDataFrame")
cov2.spix <- as(raster(cropDistanceTrondheim),"SpatialPixelsDataFrame")
proj4string(cov2.spix) <- "+proj=robin +datum=WGS84 +units=km"
cov4.spix <- as(raster(cropDistanceTrondheimWaterBody),"SpatialPixelsDataFrame")
proj4string(cov4.spix) <- "+proj=robin +datum=WGS84 +units=km"
cov3.spix <- as(raster(cropAltitudeTrondheim),"SpatialPixelsDataFrame")
covslist <- list(cov1.spix,cov2.spix,cov3.spix, cov4.spix)
covs = covslist

#choosing range and sd
#https://groups.google.com/g/r-inla-discussion-group/c/7p5n9KCfZU8
#Define SPDES
#max.edge = diff(range(spPointsGbifEdit@coords[,1])/(3*5))
range <- diff(range(mesh$loc[,2])) 
spdeEological <- list()
#for(i in 1: nspecies){
spdeEological[[1]] <- inla.spde2.pcmatern(mesh = mesh,
                                          # PC-prior on range: P(practic.range < 0.05) = 0.01
                                          prior.range = c(50, 0.1),
                                          # PC-prior on sigma: P(sigma > 0.8) = 0.01
                                          prior.sigma = c(5, 0.1))

spdeEological[[2]] <- inla.spde2.pcmatern(mesh = mesh,
                                          # PC-prior on range: P(practic.range < 0.05) = 0.01
                                          prior.range = c(50, 0.5),
                                          # PC-prior on sigma: P(sigma > 0.8) = 0.01
                                          prior.sigma = c(5, 0.1))
#}

#SPDEs for the thinning
spdeSampling <- inla.spde2.pcmatern(mesh = mesh,
                                    # PC-prior on range: P(practic.range < 0max.edge+3) = 0.01
                                    prior.range = c(range, 0.1),
                                    # PC-prior on sigma: P(sigma > 0.8) = 0.01
                                    prior.sigma = c(1, 0.1))
spdeslist <- list(spdes = spdeEological, spde2 = spdeSampling)

# Put all data together in a list
allDataForModels <- list(covariates = covs,
                         mesh = mesh,
                         boundary = spNorway,
                         dataBeforeClassification = dataBeforeClassification,
                         spPointsGbif = spPointsGbif,
                         detdata = detectionDF,
                         nspecies = nspecies,
                         spdeslist = spdeslist,
                         df = df,
                         countDf = countDf,
                         nbicDf = nbicDf,
                         occurenceDF = occurenceDF
                         )

save(allDataForModels, file = "allDataForModelNew.RData")

library(ggplot2)
library(inlabru)
names(cropDistanceTrondheim) <- "layer"
gg1 <- ggplot() +
  gg(mask(raster(cropDistanceTrondheim), allDataForModels$boundary)) +
  gg(allDataForModels$boundary, alpha = 0) +
  scale_fill_viridis_c(direction = -1)+
  #gg(spPointsGbif, shape = "+", col = "red", size = 3) +
  coord_equal()+
  theme_minimal()+
  theme(legend.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  xlab("")+
  ylab("")+
  ggtitle("Distance to road")


names(cropAltitudeTrondheim) <- "layer"
gg2 <- ggplot() +
  gg(mask(raster(cropAltitudeTrondheim), allDataForModels$boundary)) +
  gg(allDataForModels$boundary, alpha = 0) +
  scale_fill_viridis_c(direction = -1)+
  #gg(dataBeforeClassification[[1]], shape = "+", col = "red", size = 3) +
  #gg(dataBeforeClassification[[2]], shape = "+", col = "orange", size = 3) +
  #gg(spPointsGbif, shape = "+", col = "red", size = 3) +
  theme_minimal()+
  theme(legend.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  xlab("")+
  ylab("")+
  ggtitle("Elevation")

names(cropPrecipitationTrondheim) <- "layer"
gg3 <- ggplot() +
  gg(mask(raster(cropPrecipitationTrondheim), allDataForModels$boundary)) +
  gg(allDataForModels$boundary, alpha = 0) +
  scale_fill_viridis_c(direction = -1)+
  #gg(dataBeforeClassification[[1]], shape = "+", col = "red", size = 3) +
 # gg(dataBeforeClassification[[2]], shape = "+", col = "orange", size = 3) +
  #gg(spPointsGbif, shape = "+", col = "red", size = 3) +
  theme_minimal()+
  theme(legend.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  xlab("")+
  ylab("")+
  ggtitle("Annual precipitation ")

detdata <- detectionDF[[1]]%>%
  st_as_sf()%>%
  as(., "Spatial")

gg4 <- ggplot()+
  gg(allDataForModels$boundary, alpha = 0) +
  gg(detdata, shape = "*", col = "red", size = 10) +
  theme_minimal()+
  theme(legend.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  xlab("")+
  ylab("")+
  ggtitle("ToV-E sites")

gg5 <- ggplot()+
  gg(allDataForModels$boundary, alpha = 0) +
  gg(dataBeforeClassification[[1]], shape = "+", col = "red", size = 3) +
  gg(dataBeforeClassification[[2]], shape = "+", col = "orange", size = 3) +
  theme_minimal()+
  theme(legend.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  xlab("")+
  ylab("")+
  ggtitle("Selected species")

gg6 <- ggplot()+
  gg(allDataForModels$boundary, alpha = 0) +
  gg(spPointsGbif, shape = "+", col = "red", size = 3) +
  coord_equal()+
  theme_bw()+
  xlab("")+
  ylab("")+
  ggtitle("Sampled locations")

ggData <- ggpubr::ggarrange(gg5, gg4,
                            nrow = 1, 
                            ncol = 2,
                            labels = c("a)", "b)"))

ggsave("inferenceDataPlotNoSampling.png", plot = ggData)



# Add new covariate
names(cropDistanceTrondheimWaterBody) <- "layer"
gg7 <- ggplot() +
  gg(mask(raster(cropDistanceTrondheimWaterBody),allDataForModels$boundary)) +
  gg(allDataForModels$boundary, alpha = 0) +
  scale_fill_viridis_c(direction = -1)+
  #gg(spPointsGbif, shape = "+", col = "red", size = 3) +
  theme_minimal()+
  theme(legend.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  xlab("")+
  ylab("")+
  ggtitle("Distance to water body")


ggAll <- ggpubr::ggarrange(gg1, gg2, gg3, gg7,
                           nrow = 2, ncol = 2,
                           labels = c("a)", "b)", "c)", "d)"))

ggsave("inferenceCovariatesPlotNoPoints.png", plot = ggAll)


#Extract covariates at mesh points
meshSpPoints <- sp::SpatialPointsDataFrame(mesh$loc[, 1:2],
                                           data = data.frame(mesh$loc[, 3]),
                                           proj4string = CRS("+proj=robin +datum=WGS84 +units=km"))%>%
  sf::st_as_sf()
dist <- raster::extract(raster(cropDistanceTrondheim), meshSpPoints)
water <- raster::extract(cropDistanceTrondheimWaterBody, meshSpPoints)
prec <- raster::extract(cropPrecipitationTrondheim, meshSpPoints)
altt <- raster::extract(cropAltitudeTrondheim, meshSpPoints)
dist
alldata <- cbind(dist, water$layer, prec$layer, altt$layer)

corVals <- cor(na.omit(alldata))
colnames(corVals) <- c("DISTr", "DISTw", "PREC", "ELEV")
rownames(corVals) <- c("DISTr", "DISTw", "PREC", "ELEV")
fig <- ggcorrplot::ggcorrplot(corVals, lab = TRUE, type = "upper")
ggsave("corrPlots.png", plot = fig)


# plot simulation region
simRegion <- ggplot()+
  gg(allDataForModels$boundary, alpha = 0) +
  geom_rect(aes(xmin = 800,
                  xmax = 850,
                  ymin = 6700,
                  ymax = 6750),
            alpha = 0.3,
            color = "green")+
  theme_minimal()+
  theme(legend.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  xlab("")+
  ylab("")

ggsave("results/simulationRegion.png",
       plot = simRegion,
       width = 10,
       height = 5,
       units = "in")
