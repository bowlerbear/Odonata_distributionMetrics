library(tidyverse)
library(ggthemes)
library(sf)
library(cowplot)
library(ggpirate)

theme_set(theme_few())
options(scipen=10000)

#This scripts grabs the models and all the other stuff we need for later scripts

### mtbq info ####

load("splines/prep/mtbqsDF.RData")
names(mtbqsDF)[2] <- "MTB"
mtbsDF <- subset(mtbqsDF,!duplicated(MTB))

### model outputs  ####

modelSummaries <- readRDS("outputs/modelSummaries.rds")
PA_matrix <- readRDS("outputs/PA_matrix.rds")

### vectors ###

allspecies <- sort(unique(modelSummaries$Species))
allYears <- sort(unique(modelSummaries$Year))

### GIS data ####

utmProj <- "+proj=utm +zone=33 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
mtbs <- st_read(dsn="splines/prep/MTBQ_shapefile",layer="MTB_25832")
#Projected CRS: ETRS89 / UTM zone 32N
equalProj <- st_crs(mtbs)
area <- st_area(mtbs)
meanArea <- mean(area)
totalArea <- as.numeric(sum(area))

#get germany outline
germanOutline <- readRDS("gadm36_DEU_0_sp.rds") %>%
                    st_as_sf(germanOutline) %>%
                    st_transform(.,equalProj)

#add on lon and lat
modelSummaries$lon <- mtbsDF$lon[match(modelSummaries$MTB, mtbsDF$MTB)]
modelSummaries$lat <- mtbsDF$lat[match(modelSummaries$MTB, mtbsDF$MTB)]
modelSummaries <- subset(modelSummaries, !is.na(lon) & !is.na(lat))

### subset data #####

#subset to first and last year
modelSummaries_Limits <- subset(modelSummaries, Year %in% c(1990,2016))

### end ####