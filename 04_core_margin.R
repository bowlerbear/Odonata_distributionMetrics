#script to pull out core and marginal areas
library(tidyverse)
library(raster)
library(ggthemes)
theme_set(theme_few())

### get data ####

source("01_getModels.R")
source("00_functions.R")

### define core or marginal ####

coreDF <- selectSpecies2 %>%
  map(.,getCoreRegions) %>%
  reduce(rbind)

saveRDS(coreDF, file="outputs/coreDF.rds")

### sumarise core ####

coreSummary <- coreDF %>%
                group_by(Species) %>%
                summarize(nuMarginal = sum(Core=="absent" | Core=="marginal"),
                        nuCore = sum(Core=="core"),
                        sum = length(Core),
                        ratio = nuMarginal/nuCore)%>%
                rename(species = "Species")

saveRDS(coreSummary, file="outputs/coreSummary.rds")

### annual occu ####

coreDF<- readRDS("outputs/coreDF.rds")

coreAnnual <- selectSpecies2 %>%
   map(~getCoreCalc(., summary="annual")) %>%
   reduce(rbind)
saveRDS(coreAnnual, file="outputs/coreAnnual.rds")

### occu change  ####

coreAnnualChanges <- selectSpecies2 %>%
                 map(~getCoreCalc(., summary="annualchange")) %>%
                 reduce(rbind)
saveRDS(coreAnnualChanges, file="outputs/coreAnnualChanges.rds")

### occu change stats ####

coreChanges <- selectSpecies2 %>%
                 map(~getCoreCalc(., summary="change")) %>%
                 reduce(rbind)
saveRDS(coreChanges, file="outputs/coreChanges.rds")

#add on size of each region
coreSummary <- coreSummary %>%
  pivot_longer(cols=nuMarginal:ratio,
               names_to="Core",
               values_to="size") %>%
  filter(Core %in% c("nuMarginal","nuCore")) %>%
  mutate(Core = ifelse(Core=="nuMarginal", "marginal", "core")) %>%
  rename(Species = species)
  
#add on to core changes
coreChanges <- inner_join(coreChanges, coreSummary)

#fix directions when they are the same
coreChanges$direction[coreChanges$obs_change==coreChanges$null_change]<-0

#summarise them
coreChanges_national <- summariseCoreChange(coreChanges, level="national")
#excess is gain_obs - gain_null is the core

coreChanges_regional <- summariseCoreChange(coreChanges, level="regional")
#chi is positive if we gain more than expected
#chi is negative is we lose more than expected

#e.g. A affinis gained less than expected in core and more than expected in marginal

saveRDS(coreChanges_national, file="outputs/coreChanges_national.rds")
saveRDS(coreChanges_regional, file="outputs/coreChanges_regional.rds")
