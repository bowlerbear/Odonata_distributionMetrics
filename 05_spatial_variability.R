source("00_functions.R")
source("01_getModels.R")

### latitudinal extents ####

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

### Numbers of colonization and extinctions ####

 sdChanges <- lapply(allspecies,function(x){
   applyChangeSD(x, modelSummaries_Limits)
 }) %>%
   reduce(rbind)%>%
   rename(species="Species")

 summary(sdChanges)

 saveRDS(sdChanges, file="outputs/sdChanges.rds")

#### plotting ####

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

ggsave("plots/Fig.ColExt.png",width=4, height=3)

### Margin AOO/EOO #####

margins <- read.csv("splines/prep/specieslist_odonata_margins.csv", sep=";") %>%
                rename(species = Species) %>%
                inner_join(., allChanges, by="species") %>%
                mutate(margin = ifelse(Range_margin=="stable","no","yes")) %>%
                mutate(Range_margin = ifelse(Range_margin=="stable","none",Range_margin))%>%
                mutate(Range_margin = factor(Range_margin, levels=c("none",
                                                                    "northern","southern")))

table(margins$Range_margin)
#northern southern   stable 
#22        4       47
#26/73 - 36%

#### test ########

brm1 <- brm(abs(medianChange_extent) | mi(sdChange_extent) ~ Range_margin,
            data = margins, 
            save_mevars = TRUE)
summary(brm1)

brm1 <- brm(abs(medianChange_area) | mi(sdChange_area) ~ Range_margin,
            data = margins, 
            save_mevars = TRUE)
summary(brm1)

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

#### Scatterplot ####

(gScatterMarg <- ggplot(data = margins,
       aes(x = medianChange_area,
           y = medianChange_extent)) + 
      scale_colour_brewer(type="qual")+
      geom_point(aes(colour=Range_margin)) + 
      geom_hline(linetype="dashed",yintercept=0)+
      geom_vline(linetype="dashed",xintercept=0)+
      geom_abline(intercept=0, slope=1, linetype="dashed")+
      xlab("ln Change in area") + ylab("ln Change in extent"))

#### savings ###################

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

table(margins$Range_margin)
#northern southern   stable 
#22        4       47
#26/73 - 36%

#### test ########

brm1 <- brm(abs(medianChange) | mi(sdChange) ~ Range_margin,
            data = margins)
summary(brm1)

brm1 <- brm(abs(medianChange) | mi(sdChange) ~ Range_margin + Direction,
            data = margins)
summary(brm1)

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

table(margins$Range_margin)
#northern southern   stable 
#22        4       47
#26/73 - 36%

#### test ########

brm1 <- brm(abs(medianChange) | mi(sdChange) ~ Range_margin,
            data = margins)
summary(brm1)

brm1 <- brm(abs(medianChange) | mi(sdChange) ~ Range_margin + Direction,
            data = margins)
summary(brm1)

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
