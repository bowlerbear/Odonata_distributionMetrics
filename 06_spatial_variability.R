source("00_functions.R")
source("01_getModels.R")

### Latitudinal extents ####

latitudinalAnnual <- lapply(selectSpecies,function(x){
  applyRangeExtent(x, modelSummaries_Limits,summary = "annual")
}) %>% reduce(rbind)
saveRDS(latitudinalAnnual, file="outputs/latitudinalAnnual.rds")

latitudinalChanges <- lapply(selectSpecies,function(x){
  applyRangeExtent(x, modelSummaries_Limits, summary = "change")
}) %>% reduce(rbind)

saveRDS(latitudinalChanges, file="outputs/latitudinalChanges.rds")

latitudinalPosAnnual <- lapply(allspecies,function(x){
  applyRangeMean(x, modelSummaries_Limits,summary = "annual")
}) %>% reduce(rbind)
saveRDS(latitudinalPosAnnual, file="outputs/latitudinalPosAnnual.rds")

latitudinalPosChanges <- lapply(allspecies,function(x){
  applyRangeMean(x, modelSummaries_Limits, summary = "change")
}) %>% reduce(rbind)

saveRDS(latitudinalPosChanges, file="outputs/latitudinalPosChanges.rds")

#### plotting #####

latitudinalPosAnnual <- readRDS("outputs/latitudinalPosAnnual.rds") %>% 
  rename(species=Species) %>% 
  filter(species %in% selectSpecies) 
latitudinalPosAnnual$Direction <- areaChanges$Direction[match(latitudinalPosAnnual$species,areaChanges$species)]

#also get extent change
latitudinalChanges <- readRDS("outputs/latitudinalChanges.rds") %>% 
  filter(species %in% selectSpecies) 

#reorganize
latitudinalChanges_N <- latitudinalChanges %>% 
  dplyr::select(species, median_Max, sd_Max, lower_Max, upper_Max) %>%
  dplyr::rename(median = median_Max, sd = sd_Max, lower = lower_Max, upper = upper_Max) %>%
  tibble::add_column(limit="Northern")

latitudinalChanges_S <- latitudinalChanges %>% 
  dplyr::select(species, median_Min, sd_Min, lower_Min, upper_Min) %>%
  dplyr::rename(median = median_Min, sd = sd_Min, lower = lower_Min, upper = upper_Min) %>%
  tibble::add_column(limit="Southern")

latitudinalChanges <- bind_rows(latitudinalChanges_N, latitudinalChanges_S)

#### relationship with AOCC ####

allChanges <- latitudinalChanges %>%
  inner_join(.,areaChanges,
             by=c("species"))

(gAOCC <- allChanges %>%
    ggplot(aes(y= median/1000, x = medianChange, colour=limit)) + 
    geom_point() +
    geom_errorbarh(aes(xmin = lowerChange, xmax = upperChange)) + 
    geom_errorbar(aes(ymin = lower/1000, ymax = upper/1000)) +
    stat_smooth(method="gam", se=FALSE)+
    geom_hline(yintercept=0,linetype="dashed")+
    geom_vline(linetype="dashed",xintercept=0)+
    scale_colour_brewer("Limit in Germany",type="qual")+
    xlab("ln Change in area") + ylab("Latitude change (km)"))

#### N and S extents ####

gBoxLat <- latitudinalChanges  %>%
  ggplot(., aes(x = limit, y = median/1000))+
  geom_pirate(aes(colour = limit),bars=FALSE)+
  xlab("Limit in Germany") + ylab("Latitude change (km)")+
  scale_colour_brewer(type="qual")+
  geom_hline(yintercept=0, linetype="dashed")

#### plots ####

plot_grid(gBoxLat,
          gAOCC, 
          nrow=2,
          labels=c("A","B"))

ggsave("plots/Fig.S_Lat_effects.png",
       width = 4.5, height = 5.5)

### Numbers of colonization and extinctions ####

 sdChanges <- lapply(allspecies,function(x){
   applyChangeSD(x, modelSummaries_Limits)
 }) %>%
   reduce(rbind)%>%
   rename(species="Species")

 summary(sdChanges)

 saveRDS(sdChanges, file="outputs/sdChanges.rds")

#### plots ####

areaChanges <- readRDS("outputs/areaChanges.rds") %>% filter(species %in% selectSpecies)
areaChanges$Direction <- ifelse(areaChanges$medianChange>0,"Winners","Losers")
sdChanges <- readRDS("outputs/sdChanges.rds") %>% filter(species %in% selectSpecies)
allChanges <- inner_join(sdChanges,areaChanges, by=c("species"), suffix=c("_sd","_area"))

#for each species get proportion 
propSum <- allChanges %>%
  dplyr::select(species, medianStable, medianCol, medianExt) %>%
  pivot_longer(starts_with("median"),
               names_to="type",
               values_to="value") %>%
  inner_join(.,areaChanges) %>%
  arrange(medianChange)
propSum$type <- factor(propSum$type, levels=c("medianCol", "medianStable", "medianExt"))

gBar <- ggplot(propSum)+
  geom_col(aes(x=fct_reorder(species, medianChange), y=value, 
               fill=type))+
  xlab("Species") + ylab("Number of grids")+
  scale_fill_brewer("Occupancy change",palette = "PiYG", direction=-1,
                    labels=c("colonisation",
                             "always occupied",
                             "extinction")) +
  xlab("Species") + ylab("Number of grids")+
  theme(axis.text.x = element_blank()) +
  geom_vline(linetype="dashed",xintercept=28)

ggsave("plots/Fig.S_ColExt.png",width=4, height=3)

### Margin AOO/EOO #####

margins <- read.csv("splines/prep/specieslist_odonata_margins.csv", sep=";") %>%
                rename(species = Species) %>%
                inner_join(., allChanges, by="species") %>%
                mutate(margin = ifelse(Range_margin=="stable","no","yes")) %>%
                mutate(Range_margin = ifelse(Range_margin=="stable","none",Range_margin))%>%
                mutate(Range_margin = factor(Range_margin, levels=c("none",
                                                                    "northern","southern")))

#### boxplots ####

(gBoxMargArea <- margins %>%
              ggplot(., aes(x = Range_margin, y = medianChange_area))+
              geom_pirate(aes(colour = Range_margin),bars=FALSE)+
              xlab("Range_margin") + ylab("In Change in area")+
              scale_colour_brewer(type="qual")+
              xlab("Range margin within Germany") + ylab("ln Change in area") +
              geom_hline(yintercept=0, linetype="dashed"))

(gBoxMargExtent <- margins %>%
              ggplot(., aes(x = Range_margin, y = medianChange_extent))+
              geom_pirate(aes(colour = Range_margin),bars=FALSE)+
              xlab("Range_margin") + ylab("ln Change in extent")+
              scale_colour_brewer(type="qual")+
              xlab("Range margin within Germany") + ylab("ln Change in extent") +
              geom_hline(yintercept=0, linetype="dashed"))

gBoxMarg <- cowplot::plot_grid(gBoxMargArea,gBoxMargExtent, ncol=2)

#### scatterplot ####

(gScatterMarg <- ggplot(data = margins,
       aes(x = medianChange_area,
           y = medianChange_extent)) + 
      scale_colour_brewer(type="qual")+
      geom_point(aes(colour=Range_margin)) + 
      geom_hline(linetype="dashed",yintercept=0)+
      geom_vline(linetype="dashed",xintercept=0)+
      geom_abline(intercept=0, slope=1, linetype="dashed")+
      xlab("ln Change in area") + ylab("ln Change in extent"))

cowplot::plot_grid(gBoxMarg,
                   gScatterMarg,
                   labels=c("A","B"),
                   ncol=1)

ggsave("plots/Fig.S_Margin_effects_aooeoo.png",width=6, height=6)

### Margin Clumpi #####

margins <- read.csv("splines/prep/specieslist_odonata_margins.csv", sep=";") %>%
  rename(species = Species) %>%
  inner_join(., fragAnnualChanges, by="species") %>%
  mutate(margin = ifelse(Range_margin=="stable","no","yes")) %>%
  mutate(Range_margin = ifelse(Range_margin=="stable","none",Range_margin))%>%
  mutate(Range_margin = factor(Range_margin, levels=c("none",
                                                      "northern","southern")))

#### boxplot ####

gFragBox <- margins %>%
  ggplot(aes(x = Direction, y = medianChange))+
  geom_jitter(aes(colour = Range_margin),width=0.2)+
  xlab("Direction of change") + ylab("Change in Aggregation")+
  geom_hline(yintercept=0, linetype="dashed")+
  scale_colour_brewer(type="qual")

#save with saturation plot below

### Margin Saturation #####

margins <- read.csv("splines/prep/specieslist_odonata_margins.csv", sep=";") %>%
  rename(species = Species) %>%
  inner_join(., saturationChanges, by="species") %>%
  mutate(margin = ifelse(Range_margin=="stable","no","yes")) %>%
  mutate(Range_margin = ifelse(Range_margin=="stable","none",Range_margin))%>%
  mutate(Range_margin = factor(Range_margin, levels=c("none",
                                                      "northern","southern")))

#### boxplot ####

gSatBox <- margins %>%
  ggplot(aes(x = Direction, y = medianChange))+
  geom_jitter(aes(colour = Range_margin),width=0.2)+
  xlab("Direction of change") + ylab("Change in Saturation")+
  geom_hline(yintercept=0, linetype="dashed")+
  scale_colour_brewer(type="qual")

#### saving #########

cowplot::plot_grid(gFragBox,
                   gSatBox,
                   labels=c("A","B"),
                   ncol=1)

ggsave("plots/Fig.S_Margin_effects_satfrag.png",width=6, height=7)

### end ##############
