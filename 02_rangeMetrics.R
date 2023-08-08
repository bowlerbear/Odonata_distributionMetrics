### get data ####

source("01_getModels.R")
source("00_functions.R")

### area ####

areaAnnual <- lapply(allspecies, function(x){
  applyRangeArea(x,modelSummaries_Limits, summary="annual")}) %>%
  reduce(rbind) 
saveRDS(areaAnnual,file="outputs/areaAnnual.rds")

areaChanges <- lapply(allspecies, function(x){
  applyRangeArea(x,modelSummaries_Limits)})%>%
  reduce(rbind)
saveRDS(areaChanges,file="outputs/areaChanges.rds")

#cal some summary statistics
nrow(subset(areaChanges, medianChange>0))
nrow(subset(areaChanges, medianChange<0))

summary(subset(areaChanges, medianChange>0)$medianChange)
summary(subset(areaChanges, medianChange<0)$medianChange)

### extent ####

#### concave hull ####

hullChanges <- lapply(selectSpecies, function(x){
  applyConcaveMan(x,modelSummaries_Limits)})%>%
  reduce(rbind) %>%
  mutate(species = fct_reorder(species, desc(medianChange))) 
saveRDS(hullChanges,file="outputs/concavehullChanges.rds")

#cal some statistics
nrow(hullChanges)
nrow(subset(hullChanges, medianChange>0))
nrow(subset(hullChanges, medianChange<0))

summary(subset(hullChanges, medianChange>0)$medianChange)
summary(subset(hullChanges, medianChange<0)$medianChange)

### saturation #####

saturationAnnual <- lapply(allspecies,function(x){
  applySaturation(x, modelSummaries_Limits, summary = "annual")
}) %>% reduce(rbind)
saveRDS(saturationAnnual, file="outputs/saturationAnnual.rds")

saturationChanges <- lapply(selectSpecies,function(x){
  applySaturation(x, modelSummaries_Limits, summary = "change")
}) %>% reduce(rbind)
saveRDS(saturationChanges, file="outputs/saturationChanges.rds")

### end ####

