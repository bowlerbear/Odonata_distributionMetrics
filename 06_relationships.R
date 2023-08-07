library(tidyverse)
library(ggthemes)
library(cowplot)
library(ggpirate)
theme_set(theme_few())
source("00_functions.R")
options(scipen=10000)

### Fig 3 ####

areaChanges <- readRDS("outputs/areaChanges.rds")%>% filter(species %in% selectSpecies)
areaChanges$Direction <- ifelse(areaChanges$medianChange>0,"Winners","Losers")
hullChanges <- readRDS("outputs/concavehullChanges.rds")%>%filter(species %in% selectSpecies)

allChanges <- inner_join(hullChanges,areaChanges,
                         by=c("species"),
                         suffix = c("_extent","_area")) 

#### nu species ###############

#get number that increase or decrease
allChanges2 <- allChanges %>%
                  mutate(extentSig = ifelse(lowerChange_extent <0 & upperChange_extent>0, "non-signif", "signif"),
                         areaSig =  ifelse(lowerChange_area <0 & upperChange_area>0, "non-signif", "signif"),
                         extentDirection = ifelse(medianChange_extent>0, "increase", "decrease"),
                         areaDirection = ifelse(medianChange_area>0, "increase", "decrease"),
                         extentChange = paste(extentSig, extentDirection),
                         areaChange = paste(areaSig, areaDirection)) 
              
allArea <- allChanges2 %>%
              select(species, areaChange) %>%
              group_by(areaChange) %>%
              count() %>% rename(Change = areaChange) %>% 
              add_column(Type="Area")

allExtent <- allChanges2 %>%
              select(species, extentChange) %>%
              group_by(extentChange) %>%
              count() %>% rename(Change = extentChange) %>% 
              add_column(Type="Extent")

gCounts <- bind_rows(allArea, allExtent) %>%
            mutate(Change = factor(Change, levels=c("signif decrease", "non-signif decrease",
                                                    "non-signif increase", "signif increase"))) %>%
            ggplot()+
            geom_col(aes(x=Type, y=n, fill=Change))+
            scale_fill_brewer("",palette = "PiYG", direction =1) +
            xlab("Metric") + ylab("Number of species")+
            theme(legend.position = "right") +
            guides(fill=guide_legend(nrow=2,byrow=TRUE))


#### species plots ####
allChanges$species[allChanges$species == "Gomphus flavipes"] <- "Stylurus flavipes"

ggplot(allChanges)+
  geom_pointrange(aes(x=fct_reorder(species, medianChange_area), y=medianChange_area,
                      ymin = lowerChange_area,
                      ymax = upperChange_area))+
  geom_hline(yintercept=0, linetype="dashed")+
  coord_flip() +
  xlab("species") + ylab("ln Change in range area") +
  theme(axis.text.y = element_text(size=6))

ggsave("plots/Fig.S_Species_level_area.png",width=5,height=6.5)

ggplot(allChanges)+
  geom_pointrange(aes(x=fct_reorder(species, medianChange_extent), y=medianChange_extent,
                      ymin = lowerChange_extent,
                      ymax = upperChange_extent))+
  geom_hline(yintercept=0, linetype="dashed")+
  coord_flip() +
  xlab("species") + ylab("ln Change in range extent") +
  theme(axis.text.y = element_text(size=6))

ggsave("plots/Fig.S_Species_level_extent.png",width=5,height=6.5)

#### boxplots ####

(gBoxAOO <- ggplot(areaChanges, aes(x = Direction, y = medianChange))+
   geom_pirate(aes(colour = Direction),bars=FALSE)+
   xlab("Direction of change") + ylab("ln Change in area")+
   geom_hline(yintercept=0, linetype="dashed")+
   scale_colour_manual(values=c("black","black")))

hullChanges$Direction <- areaChanges$Direction[match(hullChanges$species,areaChanges$species)]

(gBoxEOO <- ggplot(hullChanges, aes(x = Direction, y = medianChange))+
    geom_pirate(aes(colour = Direction), bars=FALSE)+
    xlab("Direction of change") + ylab("ln Change in exent")+
    geom_hline(yintercept=0, linetype="dashed")+
    scale_colour_manual(values=c("black","black")))

#### AOO vs EOO ####

(fig3c <- ggplot(data = allChanges,
                 aes(x = medianChange_area,y = medianChange_extent)) + 
    geom_point() + 
    geom_smooth(method="gam",se=FALSE, colour="grey") +
    geom_errorbar(aes(ymin = lowerChange_extent,ymax = upperChange_extent)) + 
    geom_errorbarh(aes(xmin = lowerChange_area, xmax = upperChange_area))+
    geom_hline(linetype="dashed",yintercept=0)+
    geom_vline(linetype="dashed",xintercept=1.59, color="red")+
    geom_vline(linetype="dashed",xintercept=0)+
    xlab("ln Change in area") + ylab("ln Change in extent"))


(fig3c.inset <- ggplot(data = allChanges,
                 aes(x = medianChange_area,y = medianChange_extent)) + 
    geom_point(size=rel(0.5)) + 
    geom_smooth(method="lm",se=FALSE, colour="grey") +
    #geom_errorbar(aes(ymin = lowerChange_extent,ymax = upperChange_extent),size=rel(0.2)) + 
    #geom_errorbarh(aes(xmin = lowerChange_area, xmax = upperChange_area),size=rel(0.2))+
    scale_colour_viridis_c()+
    xlim(-1.5,1.6)+ylim(-2,2)+
    geom_hline(linetype="dashed",yintercept=0,size=rel(0.5))+
    geom_vline(linetype="dashed",xintercept=0,size=rel(0.5))+
    xlab("") + ylab("") +
    theme(axis.title = element_text(size=rel(0.5)),
          axis.text = element_text(size=rel(0.5))))

(plot.with.inset <-
  ggdraw() +
  draw_plot(fig3c) +
  draw_plot(fig3c.inset, x = 0.62, y = 0.16, width = .35, height = .35))

#### plot ####

#plot together
top <- plot_grid(gCounts, ncol=1, labels=c("a"))
middle <- plot_grid(gBoxAOO,gBoxEOO,ncol=2, labels=c("b","c"))
bottom <- plot.with.inset 
  
plot_grid(top,middle,bottom,nrow=3, labels=c("","","d"))

ggsave("plots/Fig.3.png",width=7.5, height=10)

### Fig. 4 ####

areaChanges <- readRDS("outputs/areaChanges.rds")
areaChanges$Direction <- ifelse(areaChanges$medianChange>0,"Winners","Losers")

#fragmentation
fragAnnualChanges <- readRDS("outputs/clumpiAnnualChange.rds") %>% filter(species %in% selectSpecies2)
fragAnnualChanges$Direction <- areaChanges$Direction[match(fragAnnualChanges$species,areaChanges$species)]

#saturation
saturationChanges <- readRDS("outputs/saturationChanges_absScale.rds") %>% filter(species %in% selectSpecies2)
saturationChanges$Direction <- areaChanges$Direction[match(saturationChanges$species,areaChanges$species)]

#### counts #####

countClumpi <- fragAnnualChanges %>%
          mutate(Sig = ifelse(lowerChange <0 & upperChange>0, "non-signif", "signif"),
                Direction = ifelse(medianChange>0, "increase", "decrease"),
              Change = paste(Sig, Direction)) %>%
          dplyr::select(species, Change) %>%
          group_by(Change) %>%
          count() %>% add_column(Type="Aggregation")

countSat <- saturationChanges %>%
  mutate(Sig = ifelse(lowerChange <0 & upperChange>0, "non-signif", "signif"),
         Direction = ifelse(medianChange>0, "increase", "decrease"),
         Change = paste(Sig, Direction)) %>%
  dplyr::select(species, Change) %>%
  group_by(Change) %>%
  count() %>% add_column(Type="Saturation")

gCounts <- bind_rows(countClumpi, countSat) %>%
  mutate(Change = factor(Change, levels=c("signif decrease", "non-signif decrease",
                                          "non-signif increase", "signif increase"))) %>%
  ggplot()+
  geom_col(aes(x=Type, y=n, fill=Change))+
  scale_fill_brewer("",palette = "PiYG", direction =1) +
  xlab("Metric") + ylab("Number of species")+
  theme(legend.position = "right") +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

#### boxplots ####

gFragBox <- fragAnnualChanges %>%
  ggplot(aes(x = Direction, y = medianChange))+
  geom_pirate(aes(colour = Direction),bars=FALSE)+
  xlab("Direction of change") + ylab("Change in Aggregation")+
  geom_hline(yintercept=0, linetype="dashed")+
  scale_colour_manual(values=c("black","black"))

#saturation
gSatBox <- saturationChanges %>%
  ggplot(aes(x = Direction, y = medianChange))+
  geom_pirate(aes(colour = Direction),bars=FALSE)+
  xlab("Direction of change") + ylab("Change in Saturation")+
  geom_hline(yintercept=0, linetype="dashed")+
  scale_colour_manual(values=c("black","black"))

#### relationship with AOCC ####

# fragmention changes
allChanges <- fragAnnualChanges %>% ungroup() %>%
  inner_join(.,areaChanges,by=c("species"),suffix = c("_frag","_area"))

(gFragAOCC <- allChanges %>%
    ggplot(aes(x = medianChange_area, y = medianChange_frag)) + 
    geom_point() + 
    geom_errorbarh(aes(xmin = lowerChange_area,xmax = upperChange_area)) + 
    geom_errorbar(aes(ymin = lowerChange_frag, ymax = upperChange_frag)) +
    stat_smooth(method="lm", se=FALSE, colour="grey")+
    geom_hline(yintercept=0,linetype="dashed")+
    geom_vline(linetype="dashed",xintercept=0)+
    xlab("ln Change in range area") + ylab("Aggregation change"))

#identify species in bottom left
subset(allChanges, medianChange_frag < (-0.15) & medianChange_area<0)

# saturation changes
allChanges <- saturationChanges %>% ungroup() %>%
  inner_join(.,areaChanges,by=c("species"),suffix = c("_sat","_area"))

(gSatAOCC <- allChanges %>%
    ggplot(aes(x = medianChange_area, y = medianChange_sat)) + 
    geom_point() + 
    geom_errorbarh(aes(xmin = lowerChange_area,xmax = upperChange_area)) + 
    geom_errorbar(aes(ymin = lowerChange_sat, ymax = upperChange_sat)) +
    stat_smooth(method="lm", se=FALSE, colour="grey")+
    geom_hline(yintercept=0,linetype="dashed")+
    geom_vline(linetype="dashed",xintercept=0)+
    xlab("ln Change in range area") + ylab("Saturation change"))

#### plot ####

middle <- plot_grid(gFragBox, gSatBox, nrow=1, labels=c("b","c"))
bottom <- plot_grid(gFragAOCC, gSatAOCC, nrow=1,labels=c("d","e"))

plot_grid(gCounts,
          middle,
          bottom,
          nrow=3,
          labels=c("a","",""))

ggsave("plots/Fig.4.png",
       width = 6.5, height = 8)

### Fig 5 #####

areaChanges <- readRDS("old/outputs_first/areaChanges.rds")%>% filter(species %in% selectSpecies2)
areaChanges$Direction <- ifelse(areaChanges$medianChange>0,"increase","decrease")

coreChanges_regional <- readRDS("old/outputs_first/coreChanges_regional.rds")%>% rename(species = Species) %>% filter(species %in% selectSpecies2)
coreChanges_regional$Direction <- areaChanges$Direction[match(coreChanges_regional$species, areaChanges$species)]

coreChanges_national <- readRDS("old/outputs_first/coreChanges_national.rds")%>% rename(species = Species) %>% filter(species %in% selectSpecies2)
coreChanges_national$Direction <- areaChanges$Direction[match(coreChanges_national$species, areaChanges$species)]

#### count ####

#regional-Chi as a measure of differences from null

#increasing species
countCore <- coreChanges_regional %>%
                  filter(Core=="core") %>%
                  mutate(Sig = ifelse(lowerCoreChi <0 & upperCoreChi>0, "non-signif", "signif"),
                  Rel = ifelse(medianCoreChi>0, "greater gains in core", "greater gains in marginal"),
                  Change = paste(Sig, Rel)) %>%
                  dplyr::select(species, Sig, Rel, Change, Core) %>%
                  inner_join(., areaChanges) %>%
                  group_by(Core, Sig, Rel, Change, Direction) %>%
                  filter(Direction=="increase") %>%
                  count() 
             
(gCounts_Inc <- countCore  %>%
  ggplot()+
  geom_col(aes(x=Rel, y=n, fill=Sig))+
  scale_fill_brewer("",palette = "Greens", direction =1) +
  xlab("Deviation from null") + ylab("Number of increasing species")+
  coord_flip()+
  theme(legend.position = "top")) 

#decreasing species
countCore <- coreChanges_regional %>%
    filter(Core=="core") %>%
    mutate(Sig = ifelse(lowerCoreChi <0 & upperCoreChi>0, "non-signif", "signif"),
           Rel = ifelse(medianCoreChi<0, "greater loss in core", "greater loss in marginal"),
           Change = paste(Sig, Rel)) %>%
    dplyr::select(species, Sig, Rel, Change, Core) %>%
    inner_join(., areaChanges) %>%
    group_by(Core, Sig, Rel, Change, Direction) %>%
    filter(Direction=="decrease") %>%
    count() 
  
(gCounts_Dec <- countCore  %>%
      ggplot()+
      geom_col(aes(x=Rel, y=n, fill=Sig)) +
      scale_fill_brewer("",palette = "PiYG", direction =-1) +
      xlab("Deviation from null") + 
      ylab("Number of decreasing species") +
      coord_flip()+
      theme(legend.position = "top"))

#### relationships with AOO ####

coreAnnualChanges <- readRDS("outputs/coreAnnualChanges.rds")

coreChanges_full <- coreAnnualChanges %>%
  rename(species = Species) %>%
  inner_join(.,areaChanges,
             by=c("species"),
             suffix = c("_core","_total")) %>%
  mutate(Core = fct_relevel(Core, "core", "marginal")) %>%
  rename(Region = Core)

(gAOO <- ggplot(data = coreChanges_full,
              aes(x = medianChange_core, y = medianChange_total)) + 
    geom_point(aes(fill=Region), shape=21) + 
    geom_errorbar(aes(ymin = lowerChange_total,ymax = upperChange_total, colour=Region)) + 
    geom_errorbarh(aes(xmin = lowerChange_core, xmax = upperChange_core, colour=Region))+
    geom_hline(linetype="dashed",yintercept=0)+
    geom_vline(linetype="dashed",xintercept=0)+
    geom_smooth(method="lm", aes(colour=Region), se=FALSE)+
    #scale_fill_brewer("Region",type="div")+
    #scale_colour_brewer("Region",type="div")+
    scale_fill_manual(values=c("black","grey"))+
    scale_colour_manual(values=c("black","grey"))+
    xlab("Change in regional range area") + ylab("Change in total range area")+
    theme_few()+
    theme(legend.position = "right"))

#### plot ####

gCounts <- plot_grid(gCounts_Dec, gCounts_Inc, nrow=1)

cowplot::plot_grid(gAOO,gCounts,
                   nrow=2,
                   labels=c("a","b"))

ggsave("plots/Fig.5.png",width=8,height=7)

### Fig. Lat effects ####

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

#### initial effects ############

#initial latitude and change in area
latitudinalPosAnnualS <- latitudinalPosAnnual %>%
                            filter(Year==1990) %>%
                              inner_join(.,areaChanges,
                              by=c("species","Direction")) 

gInitial <- ggplot(data = latitudinalPosAnnualS,
                   aes(x = median_Mean/1000000, y = medianChange)) + 
  geom_point() + 
  #geom_text(aes(label=species))+
  geom_errorbar(aes(ymin = lowerChange,ymax = upperChange)) + 
  geom_errorbarh(aes(xmin = lower_Mean/1000000, xmax = upper_Mean/1000000))+
  xlab("Initial mean latitude (1000 km)") + ylab("ln Change in area")+
  geom_hline(linetype="dashed",yintercept=0)

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

### end ######