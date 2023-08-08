library(tidyverse)
library(ggthemes)
source("00_functions.R")

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

#### model directory ####################################################

modelDirectory <- "~/Dropbox/sMon/=Odonata_stan_environ_spline/42901613"

stanFiles <- list.files(modelDirectory) %>% 
                str_subset("m_fit") %>%
                str_subset("_39")

### detection prob #####

readP <- function(file){
  
  temp <- as.data.frame(readRDS(paste(modelDirectory,file,sep="/")))
  
  #get rows for the parameter of interest
  temp$Param <- row.names(temp)
  temp <- subset(temp, grepl("mean_p",temp$Param))
  
  #get species name
  temp$File <- file
  temp$Species <- sapply(temp$File, function(x){strsplit(x, "_")[[1]][5]})
  return(temp)
  
}

modelP<- lapply(stanFiles, readP) %>% bind_rows() %>%
            janitor::clean_names()

ggplot(modelP) +
      geom_pointrange(aes(x = species, y = mean,
                          ymin = x2_5_percent,
                          ymax = x75_percent)) +
  xlab("Species") + ylab("Mean detection probability")+
  theme_classic()+
  theme(axis.text.x = element_text(size=5))+
  coord_flip() 

ggsave("plots/Fig.S_Meandetection.png",width=5,height=7)

### detection covariates ####

readPcovs <- function(file){
  
  temp <- as.data.frame(readRDS(paste(modelDirectory,file,sep="/")))
  
  #get rows for the parameter of interest
  temp$Param <- row.names(temp)
  temp <- subset(temp, grepl("beta_p",temp$Param))
  
  #get species name
  temp$File <- file
  temp$Species <- sapply(temp$File, function(x){strsplit(x, "_")[[1]][5]})
  
  #just first set
  temp <- temp[1:6,]
  
  #add names
  temp$Variable <- c("intercept","short","single","year","yday","yday2")
  
  return(temp)
  
}

modelPcovs <- lapply(stanFiles, readPcovs) %>% bind_rows() %>%
  janitor::clean_names()

ggplot(modelPcovs) +
  geom_pointrange(aes(x = species, y = mean,
                      ymin = x2_5_percent,
                      ymax = x75_percent)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
  facet_wrap(~variable) +
  coord_flip()

### national predictions #######################################################

#predictions to all MTBs

#function to apply to each file
readStanModel <- function(file){
  
  temp <- as.data.frame(readRDS(paste(modelDirectory,file,sep="/")))
  
  #get rows for the parameter of interest
  temp$Param <- row.names(temp)
  temp <- subset(temp, grepl("psi",temp$Param))
  
  #get species name
  temp$File <- file
  temp$Species <- strsplit(file, "_")[[1]][5]
  
  return(temp)
}

modelSummaries <- stanFiles %>%
  map_df(readStanModel) %>%
  mutate(siteIndex = parse_number(Param)) %>%
  filter(Species %in% selectSpecies)

summary(modelSummaries$Rhat)
mean(!is.na(modelSummaries$Rhat))

#get site and year information 
siteInfo_NAs <- readRDS("splines/prep/siteInfo_environData_NAs.rds") %>% 
    dplyr::select(MTB, x_MTB, y_MTB, siteIndex, Year) 

#merge
modelSummaries <- modelSummaries %>%
                      inner_join(.,siteInfo_NAs, by=c("siteIndex")) %>%
                      arrange(Species, siteIndex)

nuMTBs <- length(unique(siteInfo_NAs$MTB))

length(unique(modelSummaries$Species))

#### time-series ##############################################################

#analyse nationwide time series for each species

myspecies <- sort(unique(modelSummaries$Species))

annualTS <- modelSummaries %>%
            dplyr::group_by(Species,Year) %>%
            dplyr::summarise(total = sum(mean)/nuMTBs)

ggplot(annualTS)+
  geom_line(aes(x=Year, y=total))+
  facet_wrap(~Species)+
  theme_bw()

#### posterior draws ##################################################################

#pull models
stanFiles <- list.files(modelDirectory) %>% 
  str_subset("draws") 

#remove species not of interest
modelSpecies <- gsub("draws_psi_spacetime_","",stanFiles)
modelSpecies <- gsub("_39.rds","",modelSpecies)
stanFiles <- stanFiles[modelSpecies %in% selectSpecies]

#complete site and year information 
siteInfo_NAs <- readRDS("splines/prep/siteInfo_environData_NAs.rds") %>%
  dplyr::select(!c(Species,SpeciesOrig)) %>%
  dplyr::filter(type!="extension")

readPAmatrix <- function(file){
  
  temp <- readRDS(paste(modelDirectory,file,sep="/")) %>%
    t() %>%
    as.data.frame() %>%
    add_column(File = file)
  
  #add row names
  temp$File <- file
  temp$siteIndex <- siteInfo_NAs$siteIndex
  temp$Year <- siteInfo_NAs$Year
  temp$State <- siteInfo_NAs$State
  temp$MTB <- siteInfo_NAs$MTB
  temp$Species <- strsplit(file, "_")[[1]][4]
  
  #subset to first and last year
  temp <- subset(temp, Year %in% c(1990,2016))
  
  return(temp)
}

PA_matrix  <- stanFiles %>% 
                map(readPAmatrix) %>%
                bind_rows() %>%
                arrange(Species, siteIndex)

#check it aligns with the model summaries
modelSummaries_Limits <- modelSummaries %>%
                          filter(Year %in% c(1990,2016)) 

all(PA_matrix$Species == modelSummaries_Limits$Species)
all(PA_matrix$siteIndex == modelSummaries_Limits$siteIndex)
all(PA_matrix$Year == modelSummaries_Limits$Year)

### saving ####################################################################

PA_matrix <- PA_matrix  %>% select(starts_with("V"))
saveRDS(PA_matrix, file="outputs/PA_matrix.rds")

saveRDS(modelSummaries, file="outputs/modelSummaries.rds")

#### spatial maps #############################################################

allspecies <- sort(unique(modelSummaries$Species))

for(s in unique(modelSummaries$Species)){
  
  myMax <- max(modelSummaries$mean[modelSummaries$Species==s])
  
for(i in 1:length(myYears)){
  
    ggplot(subset(modelSummaries, Species == s & Year == myYears[i]))+
      geom_point(aes( x=x_MTB, y = y_MTB, colour = mean), size = 3)+
      scale_color_viridis_c("Occupancy",option = "A", direction = -1, limits=c(0,myMax))+
      theme_void()
    
ggsave(file=paste0("plots/species/Map_",s,"_",myYears[i],".png"), width=5.5, height=6)    
       
  }
}

#### first and last year ######

modelSummaries_Limits <- subset(modelSummaries, Year %in% c(1990,2016))

for(s in unique(modelSummaries_Limits$Species)){
  
  myMax <- max(modelSummaries_Limits$mean[modelSummaries_Limits$Species==s])
  
    ggplot(subset(modelSummaries_Limits, Species == s))+
      geom_point(aes( x=x_MTB, y = y_MTB, colour = mean), size = 3)+
      scale_color_viridis_c("Occupancy",option = "A", 
                            direction = -1, limits=c(0, myMax)) +
        facet_wrap(~Year) + theme_void()
    
    ggsave(file=paste0("plots/species/model_firstlast_emviron/Map_",s,"_39.png"), width=7.5, height=4.5)    

}

### end ####