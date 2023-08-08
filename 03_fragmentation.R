#install.packages("landscapemetrics")
library(landscapemetrics)
library(raster)

### get data ####

source("01_getModels.R")
source("00_functions.R")

### annual ###########

fragAnnual <- lapply(allspecies,function(x){
  applyFragStats(x, modelSummaries_Limits, summary = "annual")}) %>%
  reduce (rbind) 
saveRDS(fragAnnual, file="outputs/clumpiAnnual.rds")

### change in clumpiness ####

fragAnnual <- lapply(selectSpecies,function(x){
  applyFragStats(x, modelSummaries_Limits, summary = "change")}) %>%
  reduce (rbind) 
saveRDS(fragAnnual, file="outputs/clumpiAnnualChange.rds")
