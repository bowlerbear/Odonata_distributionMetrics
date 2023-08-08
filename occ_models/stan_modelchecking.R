library(tidyverse)
library(ggthemes)
source("00_functions.R")

### data ############################################################

myfolder <- "splines/prep"

adultData <- readRDS(paste(myfolder,"adultData.rds",sep="/"))
mtbsDF_extended <- readRDS(paste(myfolder,"MTB_stateextendedpoints.rds",sep="/")) %>%
                      filter(type!="extension") 

names(mtbsDF_extended)[1] <- "MTB"

mtbsDF <- mtbsDF_extended %>%
              filter(!is.na(MTB)) %>%
              filter(!duplicated(MTB))

environData <- readRDS(paste(myfolder, "environData.rds", sep="/")) %>%
                  rename(MTB = Value) %>% 
                  mutate(MTB = as.numeric(MTB))

mtbsDF <- inner_join(mtbsDF, environData, by="MTB")

#### site prop ##########################################################   

modelDirectory <- "~/Dropbox/sMon/=Odonata_stan_environ_spline/42905465"

stanFiles <- list.files(modelDirectory, full.names=TRUE) %>% 
  str_subset("m_fit") %>%
  str_subset("_prop") 

propData <- lapply(stanFiles, readRDS) %>%
  bind_rows()

summary(propData$Obs-propData$Sim)

(g1 <- ggplot(propData) +
  geom_point(aes(x=Obs, y=Sim))+
  geom_abline(intercept=0, slope=1, linetype="dashed")+
  xlab("Observed prop. of sites with detected species")+
  ylab("Simulated prop. of sites with detected species"))
saveRDS(g1, file="plots/bpv_prop.rds")

propData %>%
  select(Species,Obs,Sim) %>%
  pivot_longer(Obs:Sim) %>%
ggplot()+
  geom_point(aes(x=Species, y=value, colour=name))+
  scale_color_viridis_d("Type")+
  ylab("Proportion of sites with detected species")+
  coord_flip()

#tends to slightly overpredict the proportion of occupied sites....

#### end ################################################################

