library(tidyverse)
source("01_getModels.R")
source("00_functions.R")

### initial area ####

areaAnnual <- readRDS("outputs/areaAnnual.rds") %>% filter(Year==1990) %>%
  rename(initialArea = medianArea,species = Species)

areaChanges <- readRDS("outputs/areaChanges.rds") %>% 
  inner_join(.,areaAnnual) %>%
  filter(species %in% selectSpecies)

#### relationship with AOO ####

gInitial <- ggplot(data = areaChanges,
                   aes(x = initialArea/1000000, y = medianChange)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = lowerChange,ymax = upperChange)) + 
  geom_errorbarh(aes(xmin = lowerArea/1000000, xmax = upperArea/1000000))+
  xlab("Initial area (km2)") + ylab("Change in area (LRR)")+
  geom_hline(linetype="dashed",yintercept=0)

#### absolute  change #####

areaChange_abs <- readRDS("outputs_29/areaChanges_abs.rds") %>%
  setNames(gsub("Change", "Change_abs", names(.)))

areaAnnual <- readRDS("outputs/areaAnnual.rds") %>% filter(Year==1990) %>%
  rename(initialArea = medianArea,species = Species) %>%
  inner_join(., areaChange_abs, by="species")

gInitialAbsolute <- ggplot(data = areaAnnual,
                           aes(x = initialArea/1000000, y =medianChange_abs/1000000)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = lowerChange_abs/1000000, ymax = upperChange_abs/1000000)) + 
  geom_errorbarh(aes(xmin = lowerArea/1000000, xmax = upperArea/1000000))+
  xlab("Initial area (km2)") + ylab("Change in area (km2)")+
  geom_hline(linetype="dashed",yintercept=0)

#### saving ####

plot_grid(gInitial, gInitialAbsolute, 
          nrow=2, 
          labels=c("a","b"))

ggsave("plots/Fig.S_Initial_area_effects.png",width=4.5, height=6)

#only ones with small inital area have high change
summary(lm(medianChange ~ initialArea, data=areaChanges))
#slight negative effect

### extent change ####

#### initial area ####
extentChanges <- readRDS("outputs/concavehullChanges.rds") %>%
  inner_join(.,areaAnnual)

ggplot(data = extentChanges,
       aes(x = initialArea, y = medianChange)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = lowerChange,ymax = upperChange)) + 
  geom_errorbarh(aes(xmin = lowerArea, xmax = upperArea))+
  geom_hline(linetype="dashed",yintercept=0)+
  geom_vline(linetype="dashed",xintercept=0)+
  scale_x_log10()

summary(lm(medianChange ~ initialArea, data=extentChanges))
#negative effect

#### initial extent ####
extentAnnual <- readRDS("outputs/concavehullAnnual.rds") %>%
  filter(Year==1990) %>%
  rename(initialArea = medianArea,
         species = Species)

extentChanges <- readRDS("outputs/concavehullChanges.rds") %>%
  inner_join(.,extentAnnual)

ggplot(data = extentChanges,
       aes(x = initialArea, y = medianChange)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = lowerChange,ymax = upperChange)) + 
  geom_errorbarh(aes(xmin = lowerArea, xmax = upperArea))+
  geom_hline(linetype="dashed",yintercept=0)+
  geom_vline(linetype="dashed",xintercept=0)+
  scale_x_log10()

### clumpi change ####

#### initial area ####

fragChanges <- readRDS("outputs/clumpiChanges.rds") %>%
  filter(class=="change") %>%
  inner_join(.,areaAnnual, by="species")
  
ggplot(data = fragChanges,
       aes(x = initialArea, y = medianChange)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = lowerChange,ymax = upperChange)) + 
  geom_errorbarh(aes(xmin = lowerArea, xmax = upperArea))+
  geom_hline(linetype="dashed",yintercept=0)+
  geom_vline(linetype="dashed",xintercept=0)+
  scale_x_log10()

#large change with smaller initial area

#### initial clump####

fragAnnual <- readRDS("outputs/clumpiAnnual.rds") %>%
                filter(Year==1990) %>%
                rename(species = Species)

fragChanges <- readRDS("outputs/clumpiChanges.rds") %>%
  filter(class=="change") %>%
  inner_join(.,fragAnnual, by="species")

ggplot(data = fragChanges,
       aes(x = medianMetric, y = medianChange)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = lowerChange,ymax = upperChange)) + 
  geom_errorbarh(aes(xmin = lowerMetric, xmax = upperMetric))+
  geom_hline(linetype="dashed",yintercept=0)+
  geom_vline(linetype="dashed",xintercept=0)+
  scale_x_log10()
#less of a relationship
