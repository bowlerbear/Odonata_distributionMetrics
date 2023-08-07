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

### summary statistics #######

fragAnnual <- readRDS("outputs/clumpiAnnual.rds")

#mean aggregation in each time point
fragAnnual %>%
  group_by(Year) %>%
  summarise(meanMetric = median(medianMetric),
                          lower = quantile(medianMetric, 0.025),
                            upper= quantile(medianMetric, 0.975))

#how many become more clumpi
fragAnnual %>%
  group_by(Species) %>%
  summarise(clumpiIncrease = medianMetric[Year==2016]>medianMetric[Year==1990],
            clumpiChange = abs(medianMetric[Year==2016]-medianMetric[Year==1990])) %>%
  ungroup() %>%
  summarise(mean(clumpiIncrease),
            mean(clumpiChange))
