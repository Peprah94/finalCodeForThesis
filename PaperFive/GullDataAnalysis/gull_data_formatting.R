library(readr)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(gganimate)
library(parallel)
library(pbapply)
library("rnaturalearth")
library("rnaturalearthdata")
library(spatstat)
library(dplyr)
library(sp)
library(gstat)
library(ggplot2)
library(INLA)
library(rgeos)
theme_set(theme_bw())
unix::rlimit_as(100*10^9)
gull <- read_delim("~/Documents/GitHub/misclassification-project/GullDataAnalysis/gull_new_dataset/occurrence.txt",
                         delim = "\t", escape_double = FALSE,
                         trim_ws = TRUE)%>%
  mutate(truth = ifelse(scientificName == "Chroicocephalus ridibundus (Linnaeus, 1766)", "black_headed",
                        ifelse(scientificName == "Ichthyaetus melanocephalus (Temminck, 1820)", "mediterranean",
                               "laughing")))%>%
  mutate(reported = truth)%>%
  dplyr::group_by(recordedBy)%>%
  dplyr::arrange(gbifID)%>%
  dplyr::mutate(count = n(),
                previous_obs = row_number()-1)%>%
  dplyr::ungroup()



gull$reported[gull$occurrenceID %in% c("https://www.inaturalist.org/observations/136303515",
                                        "https://www.inaturalist.org/observations/135261262",
                                       "https://www.inaturalist.org/observations/130976845",
                                       "https://www.inaturalist.org/observations/128587163",
                                       "https://www.inaturalist.org/observations/128130694",
                                       "https://www.inaturalist.org/observations/124897371",
                                       "https://www.inaturalist.org/observations/120806192",
                                       "https://www.inaturalist.org/observations/120806192",
                                       "https://www.inaturalist.org/observations/119634337",
                                       "https://www.inaturalist.org/observations/118490119",
                                       "https://www.inaturalist.org/observations/118260726",
                                       "https://www.inaturalist.org/observations/118260726",
                                       "https://www.inaturalist.org/observations/118260726",
                                       "https://www.inaturalist.org/observations/117321162",
                                       "https://www.inaturalist.org/observations/113131480",
                                       "https://www.inaturalist.org/observations/113103289",
                                       "https://www.inaturalist.org/observations/112535260",
                                       "https://www.inaturalist.org/observations/112533145",
                                       "https://www.inaturalist.org/observations/112530025",
                                      "https://www.inaturalist.org/observations/112529858",
                                      "https://www.inaturalist.org/observations/112508166",
                                      "https://www.inaturalist.org/observations/109098512",
                                      "https://www.inaturalist.org/observations/109029657",
                                      "https://www.inaturalist.org/observations/108976355",
                                      "https://www.inaturalist.org/observations/108584971",
                                      "https://www.inaturalist.org/observations/107953025",
                                      "https://www.inaturalist.org/observations/105424387",
                                      "https://www.inaturalist.org/observations/105080864",
                                      "https://www.inaturalist.org/observations/103878283",
                                      "https://www.inaturalist.org/observations/103589078",
                                      "https://www.inaturalist.org/observations/103576312",
                                      "https://www.inaturalist.org/observations/101728724",
                                      "https://www.inaturalist.org/observations/101603438",
                                      "https://www.inaturalist.org/observations/100494239",
                                      "https://www.inaturalist.org/observations/100140476",
                                      "https://www.inaturalist.org/observations/99805124",
                                      "https://www.inaturalist.org/observations/98963306",
                                      "https://www.inaturalist.org/observations/98103101",
                                      "https://www.inaturalist.org/observations/98036711",
                                      "https://www.inaturalist.org/observations/96965017",
                                      "https://www.inaturalist.org/observations/94958155",
                                      "https://www.inaturalist.org/observations/92843520",
                                      "https://www.inaturalist.org/observations/92717919",
                                      "https://www.inaturalist.org/observations/91553979",
                                      "https://www.inaturalist.org/observations/91147763",
                                      "https://www.inaturalist.org/observations/90978138",
                                      "https://www.inaturalist.org/observations/90768570",
                                      "https://www.inaturalist.org/observations/88601691",
                                      "https://www.inaturalist.org/observations/86175260",
                                      "https://www.inaturalist.org/observations/83625913",
                                      "https://www.inaturalist.org/observations/77747254",
                                      "https://www.inaturalist.org/observations/77741849",
                                      "https://www.inaturalist.org/observations/76683098",
                                      "https://www.inaturalist.org/observations/72773177",
                                      "https://www.inaturalist.org/observations/72082480",
                                      "https://www.inaturalist.org/observations/71703425",
                                      "https://www.inaturalist.org/observations/71147188",
                                      "https://www.inaturalist.org/observations/71131000",
                                      "https://www.inaturalist.org/observations/46374052",
                                      "https://www.inaturalist.org/observations/34781476",
                                      "https://www.inaturalist.org/observations/34773301",
                                      "https://www.inaturalist.org/observations/34080600",
                                      "https://www.inaturalist.org/observations/30882963",
                                      "https://www.inaturalist.org/observations/28121069",
                                      "https://www.inaturalist.org/observations/27043722",
                                      "https://www.inaturalist.org/observations/23615211",
                                      "https://www.inaturalist.org/observations/20983035",
                                      "https://www.inaturalist.org/observations/18289129",
                                      "https://www.inaturalist.org/observations/20518028",
                                      "https://www.inaturalist.org/observations/11402399",
                                      "https://www.inaturalist.org/observations/11402399",
                                      "https://www.inaturalist.org/observations/10779783",
                                      "https://www.inaturalist.org/observations/10505306",
                                      "https://www.inaturalist.org/observations/10505305",
                                      "https://www.inaturalist.org/observations/10224253",
                                      "https://www.inaturalist.org/observations/4583750",
                                      "https://www.inaturalist.org/observations/3841754",
                                      "https://www.inaturalist.org/observations/3841748",
                                      "https://www.inaturalist.org/observations/4583750",
                                      "https://www.inaturalist.org/observations/3841748",
                                      "https://www.inaturalist.org/observations/3841643",
                                      "https://www.inaturalist.org/observations/1635053"
                                      )] = "mediterranean"



##############################
gull$reported[gull$occurrenceID %in% c("https://www.inaturalist.org/observations/136778974",
                                       "https://www.inaturalist.org/observations/135350496",
                                       "https://www.inaturalist.org/observations/133350012",
                                       "https://www.inaturalist.org/observations/131916754",
                                       "https://www.inaturalist.org/observations/130766681",
                                       "https://www.inaturalist.org/observations/107356990",
                                       "https://www.inaturalist.org/observations/107187959",
                                       "https://www.inaturalist.org/observations/106331423",
                                       "https://www.inaturalist.org/observations/106331423",
                                       "https://www.inaturalist.org/observations/96771483",
                                       "https://www.inaturalist.org/observations/96771480",
                                       "https://www.inaturalist.org/observations/95842807",
                                       "https://www.inaturalist.org/observations/92330815",
                                       "https://www.inaturalist.org/observations/91625182",
                                       "https://www.inaturalist.org/observations/85453012",
                                       "https://www.inaturalist.org/observations/84904724",
                                       "https://www.inaturalist.org/observations/84529080",
                                       "https://www.inaturalist.org/observations/83845651",
                                       "https://www.inaturalist.org/observations/80603854",
                                       "https://www.inaturalist.org/observations/79791734",
                                       "https://www.inaturalist.org/observations/72788030",
                                       "https://www.inaturalist.org/observations/71310037",
                                       "https://www.inaturalist.org/observations/67263594",
                                       "https://www.inaturalist.org/observations/66298463",
                                       "https://www.inaturalist.org/observations/61076099",
                                       "https://www.inaturalist.org/observations/39820318",
                                       "https://www.inaturalist.org/observations/39566244",
                                       "https://www.inaturalist.org/observations/30032281",
                                       "https://www.inaturalist.org/observations/29021728",
                                       "https://www.inaturalist.org/observations/27424086",
                                       "https://www.inaturalist.org/observations/13781047",
                                       "https://www.inaturalist.org/observations/11258766",
                                       "https://www.inaturalist.org/observations/8425050",
                                       "https://www.inaturalist.org/observations/6941948"

)] = "black_headed"

########################

gull$reported[gull$occurrenceID %in% c("https://www.inaturalist.org/observations/136270469",
                                       "https://www.inaturalist.org/observations/128355372",
                                       "https://www.inaturalist.org/observations/127914240",
                                       "https://www.inaturalist.org/observations/127204693",
                                       "https://www.inaturalist.org/observations/112126921",
                                       "https://www.inaturalist.org/observations/112126921",
                                       "https://www.inaturalist.org/observations/108325587",
                                       "https://www.inaturalist.org/observations/107581863",
                                       "https://www.inaturalist.org/observations/106448948",
                                       "https://www.inaturalist.org/observations/104658608",
                                       "https://www.inaturalist.org/observations/104127290",
                                       "https://www.inaturalist.org/observations/91241935",
                                       "https://www.inaturalist.org/observations/82668923",
                                       "https://www.inaturalist.org/observations/82513602",
                                       "https://www.inaturalist.org/observations/80709249",
                                       "https://www.inaturalist.org/observations/72591573",
                                       "https://www.inaturalist.org/observations/66817934",
                                       "https://www.inaturalist.org/observations/82513602",
                                       "https://www.inaturalist.org/observations/65088990",
                                      "https://www.inaturalist.org/observations/65088990",
                                      "https://www.inaturalist.org/observations/62953214",
                                      "https://www.inaturalist.org/observations/62730414",
                                      "https://www.inaturalist.org/observations/62730211",
                                      "https://www.inaturalist.org/observations/59647058",
                                      "https://www.inaturalist.org/observations/59643349",
                                      "https://www.inaturalist.org/observations/59636127",
                                      "https://www.inaturalist.org/observations/59633540",
                                      "https://www.inaturalist.org/observations/59629524",
                                      "https://www.inaturalist.org/observations/59627050",
                                      "https://www.inaturalist.org/observations/59625954",
                                      "https://www.inaturalist.org/observations/59613942",
                                      "https://www.inaturalist.org/observations/59613814",
                                      "https://www.inaturalist.org/observations/54776892",
                                      "https://www.inaturalist.org/observations/54776884",
                                      "https://www.inaturalist.org/observations/54776876",
                                      "https://www.inaturalist.org/observations/54776870",
                                      "https://www.inaturalist.org/observations/54776868",
                                      "https://www.inaturalist.org/observations/54776867",
                                      "https://www.inaturalist.org/observations/54288345",
                                      "https://www.inaturalist.org/observations/49204683",
                                      "https://www.inaturalist.org/observations/48248104",
                                      "https://www.inaturalist.org/observations/48248104",
                                      "https://www.inaturalist.org/observations/41289516",
                                      "https://www.inaturalist.org/observations/39937780",
                                      "https://www.inaturalist.org/observations/39745498",
                                      "https://www.inaturalist.org/observations/39351520",
                                      "https://www.inaturalist.org/observations/38220777",
                                      "https://www.inaturalist.org/observations/37954268",
                                      "https://www.inaturalist.org/observations/37954267",
                                      "https://www.inaturalist.org/observations/37954265",
                                      "https://www.inaturalist.org/observations/37740806",
                                      "https://www.inaturalist.org/observations/37740497",
                                      "https://www.inaturalist.org/observations/37222315",
                                      "https://www.inaturalist.org/observations/37221215",
                                      "https://www.inaturalist.org/observations/37189409",
                                      "https://www.inaturalist.org/observations/29901859",
                                      "https://www.inaturalist.org/observations/29390388",
                                      "https://www.inaturalist.org/observations/25115915",
                                      "https://www.inaturalist.org/observations/25001518",
                                      "https://www.inaturalist.org/observations/23941419",
                                      "https://www.inaturalist.org/observations/20944816",
                                      "https://www.inaturalist.org/observations/20802153",
                                      "https://www.inaturalist.org/observations/18576217",
                                      "https://www.inaturalist.org/observations/18515579",
                                      "https://www.inaturalist.org/observations/10696764",
                                      "https://www.inaturalist.org/observations/10150774",
                                      "https://www.inaturalist.org/observations/131464726",
                                      "https://www.inaturalist.org/observations/59628039"
)] = "laughing"

########################

gull$reported[gull$occurrenceID %in% c("https://www.inaturalist.org/observations/116887772",
                                       "https://www.inaturalist.org/observations/116886647",
                                      " https://www.inaturalist.org/observations/92769847",
                                      "https://www.inaturalist.org/observations/98036672",
                                      "https://www.inaturalist.org/observations/108584971",
                                      "https://www.inaturalist.org/observations/101320045",
                                      "https://www.inaturalist.org/observations/136197235",
                                      "https://www.inaturalist.org/observations/108126681",
                                      "https://www.inaturalist.org/observations/106180910",
                                      "https://www.inaturalist.org/observations/74052761",
                                      "https://www.inaturalist.org/observations/25270540",
                                      "https://www.inaturalist.org/observations/7143752",
                                      "https://www.inaturalist.org/observations/7143749",
                                      "https://www.inaturalist.org/observations/16742971",
                                      "https://www.inaturalist.org/observations/89621080",
                                      "https://www.inaturalist.org/observations/127878853",
                                      "https://www.inaturalist.org/observations/4954433",
                                      "https://www.inaturalist.org/observations/137321549",
                                      "https://www.inaturalist.org/observations/115610736",
                                      "https://www.inaturalist.org/observations/105568150",
                                      "https://www.inaturalist.org/observations/105143275",
                                      "https://www.inaturalist.org/observations/67395028",
                                      "https://www.inaturalist.org/observations/53604350",
                                      "https://www.inaturalist.org/observations/128119102",
                                      "https://www.inaturalist.org/observations/69044013",
                                      "https://www.inaturalist.org/observations/133530183",
                                      "https://www.inaturalist.org/observations/130952531",
                                      "https://www.inaturalist.org/observations/130800473",
                                      "https://www.inaturalist.org/observations/118562819",
                                      "https://www.inaturalist.org/observations/111117137",
                                      "https://www.inaturalist.org/observations/14275065",
                                      "https://www.inaturalist.org/observations/106597294",
                                      "https://www.inaturalist.org/observations/106441748",
                                      "https://www.inaturalist.org/observations/105365159",
                                      "https://www.inaturalist.org/observations/102715491",
                                      "https://www.inaturalist.org/observations/92457888",
                                      "https://www.inaturalist.org/observations/68014104",
                                      "https://www.inaturalist.org/observations/64365859",
                                      "https://www.inaturalist.org/observations/102714344",
                                      "https://www.inaturalist.org/observations/98883545",
                                      "https://www.inaturalist.org/observations/62883940",
                                      "https://www.inaturalist.org/observations/34311187",
                                      "https://www.inaturalist.org/observations/22080197",
                                      "https://www.inaturalist.org/observations/19696035",
                                      "https://www.inaturalist.org/observations/17430049",
                                      "https://www.inaturalist.org/observations/17430007",
                                      "https://www.inaturalist.org/observations/9338560",
                                      "https://www.inaturalist.org/observations/17316054",
                                      "https://www.inaturalist.org/observations/109127263",
                                      "https://www.inaturalist.org/observations/99076982",
                                      "https://www.inaturalist.org/observations/98963107",
                                      "https://www.inaturalist.org/observations/65759489",
                                      "https://www.inaturalist.org/observations/64514218",
                                      "https://www.inaturalist.org/observations/63635493",
                                      "https://www.inaturalist.org/observations/37167746",
                                      "https://www.inaturalist.org/observations/36218525",
                                      "https://www.inaturalist.org/observations/31343815",
                                      "https://www.inaturalist.org/observations/17792311",
                                      "https://www.inaturalist.org/observations/17792304",
                                      "https://www.inaturalist.org/observations/135261233",
                                      "https://www.inaturalist.org/observations/129832977",
                                      "https://www.inaturalist.org/observations/116250189",
                                      "https://www.inaturalist.org/observations/96575434",
                                      "https://www.inaturalist.org/observations/64515485",
                                      "https://www.inaturalist.org/observations/55779920",
                                      "https://www.inaturalist.org/observations/35403570",
                                      "https://www.inaturalist.org/observations/18125038",
                                      "https://www.inaturalist.org/observations/136955391",
                                      "https://www.inaturalist.org/observations/136033928",
                                      "https://www.inaturalist.org/observations/110397866",
                                      "https://www.inaturalist.org/observations/107774052",
                                      "https://www.inaturalist.org/observations/104969359",
                                      "https://www.inaturalist.org/observations/103921910",
                                      "https://www.inaturalist.org/observations/101085639",
                                      "https://www.inaturalist.org/observations/99912126",
                                      "https://www.inaturalist.org/observations/136937892",
                                      "https://www.inaturalist.org/observations/136036671",
                                      "https://www.inaturalist.org/observations/133586629",
                                      "https://www.inaturalist.org/observations/105909404",
                                      "https://www.inaturalist.org/observations/105838050",
                                      "https://www.inaturalist.org/observations/103825266",
                                      "https://www.inaturalist.org/observations/105838050",
                                      "https://www.inaturalist.org/observations/103825266",
                                      "https://www.inaturalist.org/observations/102460130",
                                      "https://www.inaturalist.org/observations/102157012",
                                      "https://www.inaturalist.org/observations/99804728",
                                      "https://www.inaturalist.org/observations/97923686",
                                      "https://www.inaturalist.org/observations/97539399",
                                      "https://www.inaturalist.org/observations/69956672",
                                      "https://www.inaturalist.org/observations/66091678",
                                      "https://www.inaturalist.org/observations/55307265",
                                      "https://www.inaturalist.org/observations/42675127"


)] = "other"



#ADDING ELEVATION
countries <- unique(gull$countryCode)
alt<- raster::getData('worldclim', var='alt', res=10)
bio <- raster::getData('worldclim', var='bio', res=10)
gull_data_df <- list()
for(country.tag in 1: length(countries)){
  all_data <- gull%>%
    filter(countryCode == countries[country.tag])
States <- raster::getData("GADM", country = countries[country.tag], level = 2)
USborder <- rgeos::gUnaryUnion(States, id = States$ISO)
boundary_points <- fortify(USborder)
sp_df <- SpatialPointsDataFrame(cbind(as.numeric(all_data$decimalLongitude),
                                      as.numeric(all_data$decimalLatitude)),
                                all_data,proj4string = CRS("+proj=longlat +datum=WGS84"))
df_data <- raster::intersect(sp_df, States)
gull_data <- df_data@data


# sp_df <- SpatialPointsDataFrame(cbind(as.numeric(titmouse$longitude),
#                                       as.numeric(titmouse$latitude)),
#                                 titmouse,proj4string = CRS("+proj=longlat +datum=WGS84"))
altitude <- raster::extract(alt, df_data, df=T)
bio_data <- raster::extract(bio, df_data, df = T)

gull_data_df[[country.tag]] <- cbind(gull_data, altitude, bio_data)%>%
  data.frame()
}

gull_df <- do.call("rbind", gull_data_df)

save(gull_df, file = "gull_df.RData")

#######################
# Calculate distance to road
if(!require("tidyverse")) install.packages("tidyverse")
if(!require("osmdata")) install.packages("osmdata")
if(!require("sf")) install.packages("sf")
if(!require("dodgr")) remotes::install_git("https://git.sr.ht/~mpadge/dodgr")
if(!require("geosphere")) install.packages("geosphere")
if(!require("classInt")) install.packages("classInt")
if(!require("extrafont")) install.packages("extrafont")
if(!require("ggmap")) install.packages("ggmap")
if(!require("dplyr")) install.packages("dplyr")

library(tidyverse, quietly=T) # data processing
library(osmdata, quietly=T) # load osm data
library(sf, quietly=T) # use spatial vector data
library(dodgr, quietly=T) # driving distance
library(geosphere, quietly=T) # aerial distance
library(classInt, quietly=T) # legend
library(extrafont, quietly=T) # font
library(ggmap)
library(dplyr)

#load filtered data

distance_estimate <- function(i, data){

  # message("Obtaining the states from the data")
  states <- unique(data$level2Name)
  states <- gsub("\\s*\\([^\\)]+\\)","",as.character(states))
  states <- states[states != "Rhein-Neckar-Kreis"]
  # define Belgrade's bounding box based on ggmap's bounding box
  message(paste("Extracting maps for", states[i]))
  bg_map <- ggmap::get_map(getbb(states[i]),
                           maptype = "toner-lite",
                           source = "stamen",
                           color="bw",
                           force=T)

  bg_bbox <- attr(bg_map, 'bb')

  message(paste("map coordinates for", states[i]))
  bbox <- c(xmin=bg_bbox[,2],
            ymin= bg_bbox[,1],
            xmax= bg_bbox[,4],
            ymax=bg_bbox[,3])


  #get states's paved roads
  message(paste("Getting roads for", states[i]))
  bg_way <- opq(bbox = bbox, timeout = 1240, memsize = 1004857600) %>%
    add_osm_feature(
      key = 'highway') %>%
    osmdata_sf(quiet = T)

  message(paste("Extracting the lines", states[i]))

  bg_r <- bg_way$osm_lines

  message(paste("Formatting data from", states[i], " to Polygons"))
  pk <- data %>%
    filter(level2Name == states[i]) %>%
    st_as_sf(., coords = c("decimalLongitude", "decimalLatitude"),
             crs = st_crs(bg_r))

  message(paste("Calculating distance for ", states[i]))
  dists_0 <- st_distance(pk,bg_r)
  dist2road <- apply(dists_0,1,min)

  message(paste("Returning distance for ", states[i]))
  ret <- data %>%
    filter(level2Name == states[i]) %>%
    mutate(distance = dist2road)

  return(ret)

}

states_length <- length(unique(gull_df$level2Name))
gull_data_with_distance <- formatted_data <- lapply(as.list(seq(86,100)), function(x){
    ret <- distance_estimate(x, gull_df)
  })
save(gull_data_with_distance, file = "gull_data_for_nimble.RData")



