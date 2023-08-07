browser <- list

selectSpecies <- c("Aeshna affinis","Aeshna caerulea","Aeshna cyanea","Aeshna grandis",
                   "Aeshna isoceles","Aeshna juncea","Aeshna mixta","Aeshna subarctica",         
                   "Aeshna viridis","Anax imperator","Anax parthenope","Brachytron pratense",
                   "Calopteryx splendens","Calopteryx virgo","Ceriagrion tenellum","Chalcolestes viridis",
                   "Coenagrion armatum","Coenagrion hastulatum","Coenagrion lunulatum","Coenagrion mercuriale",
                   "Coenagrion ornatum","Coenagrion puella","Coenagrion pulchellum","Coenagrion scitulum",
                   "Cordulegaster bidentata","Cordulegaster boltonii","Cordulia aenea","Crocothemis erythraea",
                   "Enallagma cyathigerum","Epitheca bimaculata","Erythromma lindenii","Erythromma najas",
                   "Erythromma viridulum","Gomphus flavipes","Gomphus pulchellus","Gomphus vulgatissimus",
                   "Ischnura elegans","Ischnura pumilio","Lestes barbarus","Lestes dryas","Lestes sponsa", 
                   "Lestes virens","Leucorrhinia albifrons","Leucorrhinia caudalis","Leucorrhinia dubia",
                   "Leucorrhinia pectoralis","Leucorrhinia rubicunda","Libellula depressa","Libellula fulva",
                   "Libellula quadrimaculata","Nehalennia speciosa","Onychogomphus forcipatus","Ophiogomphus cecilia",
                   "Orthetrum albistylum","Orthetrum brunneum","Orthetrum cancellatum","Orthetrum coerulescens",
                   "Platycnemis pennipes","Pyrrhosoma nymphula","Somatochlora alpestris","Somatochlora arctica",
                   "Somatochlora flavomaculata","Somatochlora metallica","Sympecma fusca","Sympecma paedisca",
                   "Sympetrum danae","Sympetrum depressiusculum","Sympetrum fonscolombii","Sympetrum meridionale",
                   "Sympetrum pedemontanum","Sympetrum sanguineum","Sympetrum striolatum","Sympetrum vulgatum")


selectSpecies2 <- c("Aeshna caerulea","Aeshna cyanea","Aeshna grandis",
                   "Aeshna isoceles","Aeshna juncea","Aeshna mixta","Aeshna subarctica",         
                   "Aeshna viridis","Anax imperator","Anax parthenope","Brachytron pratense",
                   "Calopteryx splendens","Calopteryx virgo","Ceriagrion tenellum","Chalcolestes viridis",
                   "Coenagrion armatum","Coenagrion hastulatum","Coenagrion lunulatum","Coenagrion mercuriale",
                   "Coenagrion ornatum","Coenagrion puella","Coenagrion pulchellum",
                   "Cordulegaster bidentata","Cordulegaster boltonii","Cordulia aenea","Crocothemis erythraea",
                   "Enallagma cyathigerum","Epitheca bimaculata","Erythromma lindenii","Erythromma najas",
                   "Erythromma viridulum","Gomphus flavipes","Gomphus pulchellus","Gomphus vulgatissimus",
                   "Ischnura elegans","Ischnura pumilio","Lestes barbarus","Lestes dryas","Lestes sponsa", 
                   "Lestes virens","Leucorrhinia albifrons","Leucorrhinia caudalis","Leucorrhinia dubia",
                   "Leucorrhinia pectoralis","Leucorrhinia rubicunda","Libellula depressa","Libellula fulva",
                   "Libellula quadrimaculata","Nehalennia speciosa","Onychogomphus forcipatus","Ophiogomphus cecilia",
                   "Orthetrum albistylum","Orthetrum brunneum","Orthetrum cancellatum","Orthetrum coerulescens",
                   "Platycnemis pennipes","Pyrrhosoma nymphula","Somatochlora alpestris","Somatochlora arctica",
                   "Somatochlora flavomaculata","Somatochlora metallica","Sympecma fusca","Sympecma paedisca",
                   "Sympetrum danae","Sympetrum depressiusculum","Sympetrum fonscolombii","Sympetrum meridionale",
                   "Sympetrum pedemontanum","Sympetrum sanguineum","Sympetrum striolatum","Sympetrum vulgatum")


lapply_with_error <- function(X,FUN,...){    
  lapply(X, function(x, ...) tryCatch(FUN(x, ...),
                                      error=function(e) NULL))
}

### occu prop #####

getOP <- function(species, modelSummaries_Limits){
  
  speciesData <- subset(modelSummaries_Limits,Species==species)
  
  speciesData %>%
    dplyr::group_by(Species,Year) %>%
    dplyr::summarise(occProp = mean(PA)) %>%
    dplyr::ungroup() %>%
    tidyr::complete(Year = allYears) %>%
    dplyr::mutate(occProp = ifelse(is.na(occProp),0,occProp)) %>% 
    dplyr::mutate(Species = species)
  
}


applyOP<- function(species, modelSummaries_Limits, summary = "change"){
  
  #apply to each realization
  temp <- lapply(1:ncol(PA_matrix),function(i){
    modelSummaries_Limits$PA <- PA_matrix[,i] 
    out <- getOP(species, modelSummaries_Limits)
    out$simNu <- i
    return(out)
  })
  temp <- do.call(rbind,temp)
  
  
  #summarise change
  if(summary == "change"){
    temp %>%
      select(Year,Species,occProp,simNu) %>%
      tidyr::pivot_wider(everything(),names_from = Year,values_from=occProp) %>%
      janitor::clean_names() %>%  
      tidyr::complete() %>%
      dplyr::mutate(x1990 = ifelse(is.na(x1990),1/2978,x1990)) %>%
      dplyr::mutate(x2016 = ifelse(is.na(x2016),1/2978,x2016)) %>%
      dplyr::mutate(x1990 = ifelse(x1990==0,1/2978,x1990)) %>%
      dplyr::mutate(x2016 = ifelse(x2016==0,1/2978,x2016)) %>%
      dplyr::mutate(change = boot::logit(x2016) - boot::logit(x1990)) %>%
      dplyr::group_by(species) %>%
      dplyr::summarise(medianChange = median(change), 
                       sdChange = sd(change),
                       lowerChange = quantile(change, 0.025),
                       upperChange = quantile(change, 0.975))
    
  }else if (summary == "annual") {
    
    #summarise annual
    temp %>%
      dplyr::group_by(Species,Year) %>%
      dplyr::summarise(medianArea = median(occProp), 
                       sdArea = sd(occProp),
                       lowerArea = quantile(occProp, 0.025),
                       upperArea = quantile(occProp, 0.975))
    
  }
  
}


### range area (AOO) ####

getRangeArea <- function(species, modelSummaries_Limits){
  
  speciesData <- subset(modelSummaries_Limits,Species==species)
  
  dat <- subset(speciesData, PA == 1)[,c("Species","PA","x_MTB","y_MTB","Year","MTB")]
  
  dat %>%
    dplyr::group_by(Species,Year) %>%
    dplyr::summarise(nuGrids = length(unique(MTB))) %>%
    dplyr::ungroup() %>%
    tidyr::complete(Year = allYears) %>%
    dplyr::mutate(Species = species) %>%
    dplyr::mutate(nuGrids = ifelse(is.na(nuGrids),0,nuGrids)) %>%
    dplyr::mutate(sumGrids = nuGrids * as.numeric(meanArea))
  
}


applyRangeArea <- function(species, modelSummaries_Limits, summary = "change"){
  
  #apply to each realization
  temp <- lapply(1:ncol(PA_matrix),function(i){
    modelSummaries_Limits$PA <- PA_matrix[,i] 
    out <- getRangeArea(species, modelSummaries_Limits)
    out$simNu <- i
    return(out)
  })
  temp <- do.call(rbind,temp)
  
  
  #summarise change
  if(summary == "change"){
  temp %>%
    select(Year,Species,sumGrids,simNu) %>%
    tidyr::pivot_wider(everything(),names_from = Year,values_from=sumGrids) %>%
    janitor::clean_names() %>%  
    tidyr::complete() %>%
    dplyr::mutate(x1990 = ifelse(is.na(x1990),as.numeric(meanArea),x1990)) %>%
    dplyr::mutate(x2016 = ifelse(is.na(x2016),as.numeric(meanArea),x2016)) %>%
    dplyr::mutate(x1990 = ifelse(x1990==0,as.numeric(meanArea),x1990)) %>%
    dplyr::mutate(x2016 = ifelse(x2016==0,as.numeric(meanArea),x2016)) %>%
    dplyr::mutate(change = log(x2016/x1990)) %>%
    dplyr::group_by(species) %>%
    dplyr::summarise(medianChange = median(change), 
                     sdChange = sd(change),
                     lowerChange = quantile(change, 0.025),
                     upperChange = quantile(change, 0.975))
    
  }else if(summary == "abs_change"){
    temp %>%
      select(Year,Species,sumGrids,simNu) %>%
      tidyr::pivot_wider(everything(),names_from = Year,values_from=sumGrids) %>%
      janitor::clean_names() %>%  
      tidyr::complete() %>%
      dplyr::mutate(x1990 = ifelse(is.na(x1990),as.numeric(meanArea),x1990)) %>%
      dplyr::mutate(x2016 = ifelse(is.na(x2016),as.numeric(meanArea),x2016)) %>%
      dplyr::mutate(x1990 = ifelse(x1990==0,as.numeric(meanArea),x1990)) %>%
      dplyr::mutate(x2016 = ifelse(x2016==0,as.numeric(meanArea),x2016)) %>%
      dplyr::mutate(change = (x2016 - x1990)) %>%
      dplyr::group_by(species) %>%
      dplyr::summarise(medianChange = median(change), 
                       sdChange = sd(change),
                       lowerChange = quantile(change, 0.025),
                       upperChange = quantile(change, 0.975))
    
  }else if (summary == "annual") {
    
  #summarise annual
   temp %>%
     dplyr::group_by(Species,Year) %>%
     dplyr::summarise(medianArea = median(sumGrids), 
                      sdArea = sd(sumGrids),
                      lowerArea = quantile(sumGrids, 0.025),
                      upperArea = quantile(sumGrids, 0.975))
    
  }
  
}

### latitudinal extents ####


getRangeMean <- function(species, modelSummaries_Limits){
  
  speciesData <- subset(modelSummaries_Limits,Species==species)
  
  dat <- subset(speciesData, PA == 1)[,c("Species","PA","x_MTB","y_MTB","Year","MTB")]
  #dat <- subset(speciesData, mean>0.1)[,c("x_MTB","y_MTB","Year","MTB")]
  
  dat %>%
    dplyr::group_by(Species,Year) %>%
    dplyr::summarise(mean_Y = mean(y_MTB)) %>%
    dplyr::ungroup()
  
}

getRangeExtents <- function(species, modelSummaries_Limits){
  
  speciesData <- subset(modelSummaries_Limits,Species==species)
  
  dat <- subset(speciesData, PA == 1)[,c("Species","PA","x_MTB","y_MTB","Year","MTB")]
  #dat <- subset(speciesData, mean>0.1)[,c("x_MTB","y_MTB","Year","MTB")]
  
  dat %>%
    dplyr::group_by(Species,Year) %>%
    dplyr::summarise(max_Y = max(y_MTB), 
                     min_Y = min(y_MTB)) %>%
    dplyr::ungroup()
  
}

applyRangeExtent <- function(species, modelSummaries_Limits, summary="change"){
  
  #apply to each realization
  temp <- lapply(1:ncol(PA_matrix),function(i){
    modelSummaries_Limits$PA <- PA_matrix[,i] 
    out <- getRangeExtents(species, modelSummaries_Limits)
    out$sim <- i
    return(out)
  })
  temp <- do.call(rbind,temp)
  
  if((sum(temp$Year==1990)) == (sum(temp$Year==2016))){
    
  #summarise change
  if(summary =="change"){
  temp %>%
    tidyr::pivot_wider(everything(),names_from = Year,values_from=c(max_Y,min_Y)) %>%
    janitor::clean_names() %>%  
    tidyr::complete() %>%
    dplyr::filter(complete.cases(.)) %>%
    dplyr::mutate(change_maxY = (max_y_2016 - max_y_1990), 
                  change_minY = (min_y_2016 - min_y_1990))  %>%
    dplyr::group_by(species) %>%
    dplyr::summarise(median_Max = median(change_maxY), 
                     sd_Max = sd(change_maxY),
                     lower_Max = quantile(change_maxY, 0.025),
                     upper_Max = quantile(change_maxY, 0.975),
                     median_Min = median(change_minY), 
                     sd_Min = sd(change_minY),
                     lower_Min = quantile(change_minY, 0.025),
                     upper_Min = quantile(change_minY, 0.975))
    
  }else if(summary =="annual"){
    
  #summarise annual
   temp %>%
     dplyr::group_by(Species,Year) %>%
     dplyr::summarise(median_Max = median(max_Y), 
                      lower_Max = quantile(max_Y, 0.025),
                      upper_Max = quantile(max_Y, 0.975),
                      median_Min = median(min_Y), 
                      lower_Min = quantile(min_Y, 0.025),
                      upper_Min = quantile(min_Y, 0.975))
  }
  }
  
}

applyRangeMean <- function(species, modelSummaries_Limits, summary="change"){
  
  #apply to each realization
  temp <- lapply(1:ncol(PA_matrix),function(i){
    modelSummaries_Limits$PA <- PA_matrix[,i] 
    out <- getRangeMean(species, modelSummaries_Limits)
    out$sim <- i
    return(out)
  })
  temp <- do.call(rbind,temp)
  
  if((sum(temp$Year==1990)>50) &  (sum(temp$Year==2016)>50)){
    
    #summarise change
    if(summary =="change"){
      temp %>%
        tidyr::pivot_wider(everything(),names_from = Year,values_from=c(mean_Y)) %>%
        janitor::clean_names() %>%  
        tidyr::complete() %>%
        dplyr::filter(complete.cases(.)) %>%
        dplyr::mutate(change_Y = (x2016 - x1990))  %>%
        dplyr::group_by(species) %>%
        dplyr::summarise(median_Mean = median(change_Y), 
                         sd_Mean = sd(change_Y),
                         lower_Mean = quantile(change_Y, 0.025),
                         upper_Mean = quantile(change_Y, 0.975))
      
    }else if(summary =="annual"){
      
      #summarise annual
      temp %>%
        dplyr::group_by(Species,Year) %>%
        dplyr::summarise(median_Mean = median(mean_Y), 
                         sd_Mean = sd(mean_Y),
                         lower_Mean = quantile(mean_Y, 0.025),
                         upper_Mean = quantile(mean_Y, 0.975))
    }
  }
}


### extent (EOCC) ####

getConvexHull <- function(species, modelSummaries_Limits){
  
  speciesData <- subset(modelSummaries_Limits,Species==species)
  
  #put hull around data for each year
  #dat <- subset(speciesData, mean>0.1)[,c("x_MTB","y_MTB","Year","MTB")]
  dat <- subset(speciesData, PA == 1)[,c("x_MTB","y_MTB","Year","MTB")]
  
  allYears <- sort(unique(speciesData$Year))
  
  rangeHull <- sapply(allYears,function(x){
    
    ydat <- dat[dat$Year==x,c("x_MTB","y_MTB")]
    
    if(nrow(ydat)>0 & length(unique(ydat$x_MTB))>2){
      ch <- chull(ydat)
      coords <- ydat[c(ch, ch[1]), ] 
      plot(ydat, pch=19, main = species)
      lines(coords, col="red")
      sp_poly <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(coords)), 
                                                       ID=1)),
                                     proj4string=sp::CRS("+proj=utm +zone=32 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
      
      #sp_poly_cut <- raster::intersect(sp_poly,germanOutline)#make sense or not??
      rangeSize <- raster::area(sp_poly)
      
      return(rangeSize)
    }else {
      return(0)
    }})
  
  data.frame(Species = species, Year = allYears, rangeHull = rangeHull) 
}


getAlphaHull <- function(species, modelSummaries_Limits){
  
  speciesData <- subset(modelSummaries_Limits,Species==species)
  
  #put hull around data for each year
  #dat <- subset(speciesData, mean>0.1)[,c("x_MTB","y_MTB","Year","MTB")]
  dat <- subset(speciesData, PA == 1)[,c("x_MTB","y_MTB","Year","MTB")]
  
  allYears <- sort(unique(speciesData$Year))
  
  rangeAlpha <- sapply(allYears,function(x){
    
    ydat <- dat[dat$Year==x,c("x_MTB","y_MTB")]
    
    if(nrow(ydat)>0 & length(unique(ydat$x_MTB))>2){
      
      dist <- 200000
      ah <- alphahull::ahull(ydat, alpha = dist)
      #plot(ah, main = species)
      alphahull::areaahull(ah)
      
    }else {
      return(0)
    }})
  
  data.frame(Species = species, Year = allYears, rangeAlpha = rangeAlpha) 
  
}


getMCP <- function(species, modelSummaries_Limits){
  
  require(adehabitatHR)
  require(sp)
  
  speciesData <- subset(modelSummaries_Limits,Species==species)
  
  #put hull around data for each year
  #dat <- subset(speciesData, mean>0.1)[,c("x_MTB","y_MTB","Year","MTB")]
  dat <- subset(speciesData, PA == 1)[,c("x_MTB","y_MTB","Year","MTB")]
  
  allYears <- sort(unique(speciesData$Year))
  
  rangeHull <- sapply(allYears,function(x){
    
    ydat <- dat[dat$Year==x,c("x_MTB","y_MTB")]
    
    if(nrow(ydat)>0 & length(unique(ydat$x_MTB))>4){
      
      mycoords <- ydat[,c("x_MTB","y_MTB")]
      coordinates(mycoords) <- c("x_MTB","y_MTB")
      proj4string(mycoords) <- CRS("+proj=utm +zone=32 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
      ch <- mcp(mycoords, percent=99, unout="m2")
      rangeSize <- ch@data$area
      return(rangeSize)
    }else {
      return(0)
    }})
  
  data.frame(Species = species, Year = allYears, rangeHull = rangeHull) 
}

getConcaveMan <- function(species, modelSummaries_Limits){
  
  speciesData <- subset(modelSummaries_Limits,Species==species)
  
  #put hull around data for each year
  #dat <- subset(speciesData, mean>0.1)[,c("x_MTB","y_MTB","Year","MTB")]
  dat <- subset(speciesData, PA == 1)[,c("x_MTB","y_MTB","Year","MTB")]
  
  allYears <- sort(unique(speciesData$Year))
  
  rangeMan <- sapply(allYears,function(x){
    
    ydat <- dat[dat$Year==x,c("x_MTB","y_MTB")]
    
    if(nrow(ydat)>0 & length(unique(ydat$x_MTB))>2){
      
      ydat <- sf::st_as_sf(ydat, coords =c("x_MTB", "y_MTB"),crs = 25832)
      sp_poly <- concaveman::concaveman(ydat)
      #germanOutline_sf <- st_transform(germanOutline, st_crs(ydat))
      #sp_poly_cut <- st_intersection(sp_poly, germanOutline_sf)
      #plot(sp_poly_cut)
      return(as.numeric(st_area(sp_poly)))#m2 units
    }else {
      return(0)
    }})
  
  data.frame(Species = species, Year = allYears, rangeMan = rangeMan)
  
}

applyConcaveMan <- function(species, modelSummaries_Limits, summary = "change"){
  
  #apply to each realization
  temp <- lapply(1:ncol(PA_matrix),function(i){
    
    modelSummaries_Limits$PA <- PA_matrix[,i] 
    out <- getConcaveMan(species, modelSummaries_Limits)
    out$sim <- i
    return(out)
  })
  temp <- do.call(rbind,temp)
  
  #summarise change
  if(summary == "change"){
  temp %>%
    tidyr::pivot_wider(everything(),names_from = Year, values_from = rangeMan) %>%
    janitor::clean_names() %>%  
    dplyr::mutate(change = log((x2016+as.numeric(meanArea))/(x1990+as.numeric(meanArea))))  %>%
    dplyr::group_by(species) %>%
    dplyr::summarise(medianChange = median(change), 
                     sdChange = sd(change), 
                     lowerChange = quantile(change, 0.025),
                     upperChange = quantile(change, 0.975))
    
  }else if(summary == "annual"){
    
    #summarise annual
     temp %>%
       dplyr::group_by(Species,Year) %>%
       dplyr::summarise(medianArea = median(rangeMan), 
                        sdArea = sd(rangeMan), 
                        lowerArea = quantile(rangeMan, 0.025),
                        upperArea = quantile(rangeMan, 0.975))
    
  }
  
}


applyAlphaHull <- function(species, modelSummaries_Limits, summary = "change"){
  
  #apply to each realization
  temp <- lapply(1:ncol(PA_matrix),function(i){
    
    modelSummaries_Limits$PA <- PA_matrix[,i] 
    out <- getAlphaHull(species, modelSummaries_Limits)
    #print(i)
    out$sim <- i
    return(out)
  })
  temp <- do.call(rbind,temp)
  
  #summarise change
  if(summary == "change"){
    temp %>%
      tidyr::pivot_wider(everything(),names_from = Year, values_from = rangeAlpha) %>%
      janitor::clean_names() %>%  
      dplyr::mutate(change = log((x2016+10)/(x1990+10)))  %>%
      dplyr::group_by(species) %>%
      dplyr::summarise(medianChange = median(change,na.rm=T), 
                       sdChange = sd(change,na.rm=T), 
                       lowerChange = quantile(change, 0.025,na.rm=T),
                       upperChange = quantile(change, 0.975,na.rm=T))
    
  }else if(summary == "annual"){
    
    #summarise annual
    temp %>%
      dplyr::group_by(Species,Year) %>%
      dplyr::summarise(medianArea = median(rangeAlpha,na.rm=T),
                       sdArea = sd(rangeAlpha,na.rm=T),
                       lowerArea = quantile(rangeAlpha, 0.025,na.rm=T),
                       upperArea = quantile(rangeAlpha, 0.975,na.rm=T))
    
  }
  
}


applyMCPHull <- function(species, modelSummaries_Limits, summary = "change"){
  
  #apply to each realization
  temp <- lapply(1:ncol(PA_matrix),function(i){
    
    modelSummaries_Limits$PA <- PA_matrix[,i] 
    out <- getMCP(species, modelSummaries_Limits)
    out$sim <- i
    return(out)
  })
  temp <- do.call(rbind,temp)
  
  #summarise change
  if(summary == "change"){
    temp %>%
      tidyr::pivot_wider(everything(),names_from = Year, values_from = rangeHull) %>%
      janitor::clean_names() %>%  
      dplyr::mutate(change = log10((x2016+100000)/(x1990+100000)))  %>%
      dplyr::group_by(species) %>%
      dplyr::summarise(medianChange = median(change),
                       sdChange = sd(change),
                       lowerChange = quantile(change, 0.025),
                       upperChange = quantile(change, 0.975))
    
  }else if(summary == "annual"){
    
    #summarise annual
    temp %>%
      dplyr::group_by(Species,Year) %>%
      dplyr::summarise(medianArea = median(rangeHull),
                       sdArea = sd(rangeHull),
                       lowerArea = quantile(rangeHull, 0.025),
                       upperArea = quantile(rangeHull, 0.975))
    
  }
  
}


compareHulls <- function(species, modelSummaries_Limits, myYear=2016){
  
  #pick species
  speciesData <- subset(modelSummaries_Limits,Species==species)
  
  #generate binary data
  speciesData$PA <- sapply(speciesData$mean, function(x) rbinom(1,1,x))
  
  #pick year
  speciesData <- subset(speciesData, Year == myYear)
  
  #get presences
  dat <- subset(speciesData, PA == 1)[,c("x_MTB","y_MTB","Year","MTB")]
  
    if(nrow(dat)>0 & length(unique(dat$x_MTB))>2){
      
      require(sp)
      require(sf)
          
      #convex hull
      ydat <- dat[,c("x_MTB","y_MTB")]
      ch <- chull(ydat)
      coords <- ydat[c(ch, ch[1]), ] 
      par(mfrow=c(2,2))
      plot(ydat, pch=19, main = "convex hull")
      lines(coords, col="red")
      
      #alpha hull
      ydat <- dat[,c("x_MTB","y_MTB")]
      dist <- 200000
      ah <- alphahull::ahull(ydat, alpha = dist)
      plot(ah, main = "alpha hull")
      #such that an edge of a disk of radius 1/a can be drawn 
      #between any two edge members of a set 
      #of points and still contain all the points.
      
      #minimum convex hull
      ydat <- dat[,c("x_MTB","y_MTB")]
      mycoords <- ydat[,c("x_MTB","y_MTB")]
      coordinates(mycoords) <- c("x_MTB","y_MTB")
      proj4string(mycoords) <- CRS("+proj=utm +zone=32 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
      ch <- adehabitatHR::mcp(mycoords, percent=99, unout="m2")
      plot(ydat, pch=19, main = "mch")
      plot(ch,add=T)
      
      #concaveman
      ydat <- as.matrix(dat[,c("x_MTB","y_MTB")])
      sp_poly <- concaveman::concaveman(ydat)
      plot(ydat, pch=19, main = "concave man")
      lines(sp_poly)
    
    } else{
      print("insufficient data")
    }
  }



### spatial var ####

getChangeSD <- function(species, modelSummaries_Limits){
  
  speciesData <- subset(modelSummaries_Limits,Species==species)
  e <- 1
  
  speciesData %>%
    dplyr::group_by(Species) %>%
    dplyr::summarise(nuStable = sum(PA[Year==1990]==PA[Year==2016]),
                     nu1990 = sum(PA[Year==1990]),
                     nu2016 = sum(PA[Year==2016]),
                     nuStable0 = sum(PA[Year==1990]==0 & PA[Year==2016]==0),
                     nuStable1 = sum(PA[Year==1990]==1 & PA[Year==2016]==1),
                     nuIncrease = sum(PA[Year==1990]==0 & PA[Year==2016]==1),
                     nuDecrease = sum(PA[Year==1990]==1 & PA[Year==2016]==0),
                     totalGridChange = (nuIncrease + nuDecrease),
                     total = length(PA[Year==1990]),
                     totalChange = sum(PA[Year==2016])-sum(PA[Year==1990])) %>%
    #Site turnover is defined as the probability that a 
    #randomly chosen occupied site is newly occupied
    dplyr::mutate(propCol = (nuIncrease)/(total - nu1990 + 1),
                  propExt = (nuDecrease)/(nu1990+1),
                  Turnover = (nuIncrease)/(nu2016+1))#proportion of occupied patchs that are new
  
}


applyChangeSD <-function(species, modelSummaries_Limits, summary = "change"){
  
  #apply to each realization
  temp <- lapply(1:ncol(PA_matrix),function(i){
    
    modelSummaries_Limits$PA <- PA_matrix[,i] 
    out <- getChangeSD(species, modelSummaries_Limits)
    out$simNu <- i
    return(out)
  }) %>%
    reduce(rbind) %>%
    dplyr::group_by(Species) %>%
    dplyr::summarise(medianTurnover = quantile(Turnover, 0.5),
                     sdTurnover = sd(Turnover),
                     lowerTurnover = quantile(Turnover, 0.025),
                     upperTurnover = quantile(Turnover, 0.975),
                     medianExtp = quantile(propExt, 0.5),
                     sdExtp = sd(propExt),
                     lowerExtp = quantile(propExt, 0.025),
                     upperExpt = quantile(propExt, 0.975),
                     medianColp = quantile(propCol, 0.5),
                     sdColp = sd(propCol),
                     lowerColp = quantile(propCol, 0.025),
                     upperColp = quantile(propCol, 0.975),
                     medianStable = median(nuStable1),
                     medianCol = quantile(nuIncrease, 0.5),
                     sdCol = sd(nuIncrease),
                     lowerCol = quantile(nuIncrease, 0.025),
                     upperCol = quantile(nuIncrease ,0.975),
                     medianExt = quantile(nuDecrease, 0.5),
                     sdExt = sd(nuDecrease),
                     lowerExt = quantile(nuDecrease, 0.025),
                     upperExt = quantile(nuDecrease ,0.975))

  
}


### clumpiness ####

#calculate the local occupancy change and get the clumpiness of the change

getFragChange <- function(species, modelSummaries_Limits){
  
  #data for 1990
  speciesSummary <- modelSummaries_Limits %>%
    filter(Species==species) %>%
    filter(Year == 1990)
  
  speciesRaster <- SpatialPixelsDataFrame(points = as.matrix(speciesSummary[,c("lon","lat")]),
                                          data = data.frame(PA=speciesSummary$PA),
                                          tolerance = 0.6,
                                          proj4string = crs("+proj=longlat +datum=WGS84"))
  
  #make into a raster
  r1 <- raster(speciesRaster)
  
  #data for 2016
  speciesSummary <- modelSummaries_Limits %>%
    filter(Species==species) %>%
    filter(Year == 2016)
  
  speciesRaster <- SpatialPixelsDataFrame(points = as.matrix(speciesSummary[,c("lon","lat")]),
                                          data = data.frame(PA=speciesSummary$PA),
                                          tolerance = 0.6,
                                          proj4string = crs("+proj=longlat +datum=WGS84"))
  
  #make into a raster
  r2 <- raster(speciesRaster)
  
  
  #make stack  
  speciesStack <- stack(list(r1,r2))
  #speciesStack <- projectRaster(speciesStack, crs=utmProj, method="ngb")
  speciesChange <- calc(speciesStack, function(x) x[2]-x[1])
  
  #calc change metrics
  changes <- calculate_lsm(speciesChange, 
                what = c("lsm_c_pland",
                         "lsm_c_clumpy"),
                full_name = TRUE) %>%
    mutate(class = ifelse(class==0, "no change",
                          ifelse(class==1, "increase", "decrease"))) %>%
    add_column(Species = species)
  
  #also absolute change
  speciesAbsChange <- calc(speciesChange, function(x) ifelse(x==0,0,1))
  
  calculate_lsm(speciesAbsChange, 
                           what = c("lsm_c_pland",
                                    "lsm_c_clumpy"),
                           full_name = TRUE) %>%
    filter(class == 1) %>%
    mutate(class = "change") %>%
    add_column(Species = species) %>%
    bind_rows(.,changes)
  
}


applyFragChange <- function(species, modelSummaries_Limits){
  
  #apply to each realization
  temp <- lapply(1:ncol(PA_matrix),function(i){
    
    modelSummaries_Limits$PA <- PA_matrix[,i] 
    out <- getFragChange(species, modelSummaries_Limits)
    out$simNu <- i
    return(out)
  })
  temp <- do.call(rbind,temp)
  
  #summarise across sims
  temp %>%
      dplyr::filter(metric == "clumpy" & !is.na(value)) %>%
      dplyr::select(Species,class,value,simNu) %>%
      dplyr::group_by(Species,class) %>%
      dplyr::summarise(medianChange = median(value), 
                       sdChange = sd(value), 
                       lowerChange = quantile(value, 0.025),
                       upperChange = quantile(value, 0.975))
    
}


#get the clumpiness of the annual distributions, and then the change in the clumpiness
getFragStats <- function(species, modelSummaries_Limits){
  
  #data for 1990
  speciesSummary <- modelSummaries_Limits %>%
    filter(Species==species) %>%
    filter(Year == 1990)
  
  speciesRaster <- SpatialPixelsDataFrame(points = as.matrix(speciesSummary[,c("lon","lat")]),
                                          data = data.frame(PA=speciesSummary$PA),
                                          tolerance = 0.6,
                                          proj4string = crs("+proj=longlat +datum=WGS84"))
  
  #if all values are zero, insert one occupied site
  if(all(speciesRaster$PA==0)){
    ind <- sample(1:length(speciesRaster$PA),1)
    speciesRaster$PA[ind] <- 1
  }
  
  #make into a raster
  r1 <- raster(speciesRaster)
  
  out1 <- calculate_lsm(r1, 
                what = c("lsm_c_clumpy"),
                full_name = TRUE) %>%
    filter(class == 1) %>%
    filter(name == "clumpiness index") %>%
    add_column(Year = 1990) %>%
    add_column(Species = species)
  
  #data for 2016
  speciesSummary <- modelSummaries_Limits %>%
    filter(Species==species) %>%
    filter(Year == 2016)
  
  speciesRaster <- SpatialPixelsDataFrame(points = as.matrix(speciesSummary[,c("lon","lat")]),
                                          data = data.frame(PA=speciesSummary$PA),
                                          tolerance = 0.6,
                                          proj4string = crs("+proj=longlat +datum=WGS84"))
  
  if(all(speciesRaster$PA==0)){
    ind <- sample(1:length(speciesRaster$PA),1)
    speciesRaster$PA[ind] <- 1
  }
  
  #make into a raster
  r2 <- raster(speciesRaster)
  
  #make stack  
  #speciesStack <- stack(list(r1,r2))
  #speciesStack <- projectRaster(speciesStack, crs=utmProj, method="ngb")
  
  #calc metrics
  out2 <- calculate_lsm(r2, 
                        what = c("lsm_c_clumpy"),
                        full_name = TRUE) %>%
    filter(class == 1) %>%
    filter(name == "clumpiness index") %>%
    add_column(Year = 2016) %>%
    add_column(Species = species)
  
  rbind(out1,out2)
  
}

applyFragStats <- function(species, modelSummaries_Limits, summary = "change"){
  
  #apply to each realization
  temp <- lapply(1:ncol(PA_matrix),function(i){
    modelSummaries_Limits$PA <- PA_matrix[,i] 
    out <- getFragStats(species, modelSummaries_Limits)
    out$simNu <- i
    return(out)
  })
  temp <- do.call(rbind,temp)
  
  #summarise change
  if(summary == "change"){
    temp %>%
      dplyr::filter(class == 1 & metric == "clumpy") %>%
      dplyr::select(Species,Year,value,simNu) %>%
      tidyr::pivot_wider(everything(),names_from = Year, values_from = value) %>%
      janitor::clean_names() %>%  
      mutate(x1990 = ifelse(is.infinite(x1990), 0,x1990)) %>% #insert zero when missing
      mutate(x1990 = ifelse(is.na(x1990), 0,x1990)) %>% #insert zero when missing
      mutate(x1990 = ifelse(x1990<0, 0,x1990)) %>% # insert zero when missing
      mutate(x2016 = ifelse(is.infinite(x2016), 0,x2016)) %>% #insert zero when missing
      mutate(x2016 = ifelse(is.na(x2016), 0,x2016)) %>% #insert zero when missing
      mutate(x2016 = ifelse(x2016<0, 0,x2016)) %>% # insert zero when missing
      dplyr::mutate(change = x2016 - x1990)  %>%
      dplyr::group_by(species) %>%
      dplyr::summarise(medianChange = median(change),
                       sdChange = sd(change),
                       lowerChange = quantile(change, 0.025),
                       upperChange = quantile(change, 0.975))
    
  }else if(summary == "annual"){
    
    #summarise annual
    temp %>%
      dplyr::group_by(Species,Year,class, metric,name,type,function_name) %>%
      mutate(value = ifelse(is.na(value), 0,value)) %>% #insert zero when missing
      mutate(value = ifelse(value<0, 0,value)) %>% # insert zero when missing
      dplyr::summarise(medianMetric = median(value), 
                       sdMetric = sd(value), 
                       lowerMetric = quantile(value, 0.025),
                       upperMetric = quantile(value, 0.975))
    
  }
  
}


### core regions ####

getCoreRegions <- function(myspecies){
  
  speciesSummary <- modelSummaries %>%
    filter(Species == myspecies) %>%
    filter(Year %in% 1990:1995) %>% #2002:2007
    dplyr::group_by(Species, MTB, lon, lat) %>%
    dplyr::summarise(meanpsi = mean(mean)) %>%
    dplyr::mutate(normalpsi = (meanpsi-min(.$meanpsi))/(max(.$meanpsi)-min(.$meanpsi)))
  
  speciesPixels <- SpatialPixelsDataFrame(points = as.matrix(speciesSummary[,c("lon","lat")]),
                                          data = data.frame(PA=speciesSummary$normalpsi),
                                          tolerance = 0.6,
                                          proj4string = crs("+proj=longlat +datum=WGS84"))
  
  #make into a raster
  speciesRaster <- raster(speciesPixels)
  
  #for each raster cell, get mean class of surroundings
  ## in a 3-by-3 window surrounding each focal cell 
  rmean <- focal(speciesRaster, 
                 w = matrix(1, ncol=3, nrow=3), 
                 fun=mean, na.rm=TRUE)
  #plot(rmean)
  
  #define core as all those surrounded by 50% presences
  speciesRasterCore <- rmean
  speciesRasterCore[speciesRasterCore >= 0.5] <- 1 #indicator for core sites
  speciesRasterCore[speciesRasterCore < 0.5] <- 0.5 #indicator for marginal sites
  speciesRasterCore[is.na(speciesRaster)] <- NA
  #speciesRasterCore[speciesRaster < 0.1] <- 0 #absent sites
  plot(speciesRasterCore)
  
  
  # #plot example
  # library(tmap)
  #  t1 <- tm_shape(speciesRaster)+
  #    tm_raster(style="cont",legend.show=FALSE)
  #  t2 <- tm_shape(speciesRasterCore)+
  #    tm_raster(legend.show=FALSE)
  #  tmap_arrange(t1,t2,ncol=2)
  
  
  #convert back into a data frame
  coreDF <- as.data.frame(speciesRasterCore, xy=TRUE)
  coreDF$Obs <- as.data.frame(speciesRaster)[,1] #original observation
  coreDF$cellNu <- cellFromXY(speciesRaster, coreDF[,c("x","y")])
  coreDF <- subset(coreDF, !is.na(layer)) #remove sites beyond german border
  
  #clean marginal information
  coreDF$Core <- ifelse(coreDF$layer==1,"core",
                        ifelse(coreDF$layer==0.5, "marginal", "absent"))
  
  #add on other data
  coreDF$Species <- myspecies
  
  
  #map to mtbqs
  mtbPixels <- SpatialPixelsDataFrame(points = as.matrix(speciesSummary[,c("lon","lat")]),
                                      data = data.frame(MTB = as.numeric(speciesSummary$MTB)),
                                      tolerance = 0.6,
                                      proj4string = crs("+proj=longlat +datum=WGS84"))
  
  mtbRaster <- raster(mtbPixels)
  mtbDF <- as.data.frame(mtbRaster, xy=TRUE)
  mtbDF$cellNu <- cellFromXY(mtbRaster, mtbDF[,c("x","y")])
  coreDF <- left_join(coreDF, mtbDF[,3:4], by="cellNu")
  
  return(coreDF)
  
}

# getCoreCalc <- function(myspecies, summary="annual"){
# 
#   coreDF_species <- subset(coreDF, Species== myspecies)
#   modelSummaries_Limits$Core <- coreDF_species$Core[match(modelSummaries_Limits$MTB,
#                                                           coreDF_species$MTB)]
#   modelSummaries_Limits$Core[modelSummaries_Limits$Core=="absent"] <- "marginal"
# 
#   #for each realization do the following
#   temp <- lapply(1:dim(PA_matrix)[2], function(i){
#     modelSummaries_Limits %>%
#       add_column(PA = PA_matrix[,i]) %>%
#       filter(!is.na(Core)) %>%
#       filter(Species==myspecies) %>%
#       group_by(Core, Year) %>%
#       summarize(occ = sum(PA), total= length(PA)) %>%
#       mutate(prop = occ/total) %>%
#       add_column(simNu = i)
#   }) %>% bind_rows() %>% ungroup()
# 
# 
#   if(summary=="annual"){
# 
#     temp %>%
#       dplyr::group_by(Core,Year) %>%
#       dplyr::summarize(medianOcc = quantile(occ,0.5),
#                        lowerOcc = quantile(occ,0.025),
#                        upperOcc = quantile(occ,0.975),
#                        medianUnocc = quantile(total-occ,0.5),
#                        lowerUnocc = quantile(total-occ,0.025),
#                        upperUnocc = quantile(total-occ,0.975),
#                        medianTotal = quantile(total,0.5),
#                        lowerTotal = quantile(total,0.025),
#                        upperTotal = quantile(total,0.975),
#                        medianProp = quantile(prop,0.5),
#                        sdProp = sd(prop),
#                        lowerProp = quantile(prop,0.025),
#                        upperProp = quantile(prop,0.975)) %>%
#       add_column(Species = myspecies)
# 
#   }else if(summary=="change"){
# 
#     temp %>%
#       dplyr::select(Core, Year, prop, simNu) %>%
#       dplyr::group_by(Core,Year, simNu) %>%
#       tidyr::pivot_wider(names_from="Year", values_from="prop") %>%
#       janitor::clean_names(case = "title") %>%
#       dplyr::mutate(change = boot::logit(X2016+0.01) - boot::logit(X1990+0.01)) %>%
#       dplyr::mutate(change = ifelse(is.infinite(change),0, change)) %>%
#       dplyr::mutate(change = ifelse(is.na(change),0, change)) %>%
#       dplyr::group_by(Core) %>%
#       dplyr::summarize(medianChange = quantile(change,0.5),
#                        sdChange = sd(change),
#                        lowerChange = quantile(change,0.025),
#                        upperChange = quantile(change,0.975)) %>%
#       add_column(Species = myspecies)
# 
#   }
# 
# }

getCoreCalc <- function(myspecies, summary="annual"){
  
  coreDF_species <- subset(coreDF, Species== myspecies)
  
  modelSummaries_Limits$Core <- coreDF_species$Core[match(modelSummaries_Limits$MTB, 
                                                          coreDF_species$MTB)]
  modelSummaries_Limits$Core[modelSummaries_Limits$Core=="absent"] <- "marginal"
  
  #for each realization do the following
  temp <- lapply(1:dim(PA_matrix)[2], function(i){
    modelSummaries_Limits %>%
      add_column(PA = PA_matrix[,i]) %>%
      filter(!is.na(Core)) %>%
      filter(Species==myspecies) %>%
      group_by(Core, Year) %>%
      summarize(occ = sum(PA), total= length(PA)) %>%
      mutate(prop = occ/total) %>%
      add_column(simNu = i)
  }) %>% bind_rows() %>% ungroup() 
  
  
  if(summary=="annual"){
    
    temp %>%
      dplyr::group_by(Core,Year) %>%
      dplyr::summarize(medianTotal = median(total),
                      medianCoreOcc = quantile(occ,0.5),
                       sdCoreOcc = sd(occ),
                       lowerCoreOcc = quantile(occ,0.025),
                       upperCoreOcc = quantile(occ,0.975),
                       medianCoreProp = quantile(prop,0.5),
                       sdCoreProp = sd(prop),
                       lowerCoreProp = quantile(prop,0.025),
                       upperCoreProp = quantile(prop,0.975)) %>%
      add_column(Species = myspecies)
    
    
    } else if(summary=="annualchange"){
      
      temp %>%
            dplyr::select(Core, Year, prop, simNu) %>%
            dplyr::group_by(Core,Year, simNu) %>%
            tidyr::pivot_wider(names_from="Year", values_from="prop") %>%
            janitor::clean_names(case = "title") %>%
            dplyr::mutate(change = boot::logit(X2016+0.01)-boot::logit(X1990+0.01)) %>%
            dplyr::mutate(change = ifelse(is.infinite(change),0, change)) %>%
            dplyr::mutate(change = ifelse(is.na(change),0, change)) %>%
            dplyr::group_by(Core) %>%
            dplyr::summarize(medianChange = quantile(change,0.5),
                             sdChange = sd(change),
                             lowerChange = quantile(change,0.025),
                             upperChange = quantile(change,0.975)) %>%
            add_column(Species = myspecies)
      
  }else if(summary=="change"){
    
    #get total change
    changeTemp <- temp %>% 
      dplyr::select(Core, Year, occ, total, simNu) %>%
      dplyr::group_by(Core, Year, simNu) %>%
      tidyr::pivot_wider(names_from="Year", values_from=c('occ','total')) %>%
      janitor::clean_names() %>%
      rename(Core = core) %>%
      dplyr::group_by(sim_nu) %>%
      dplyr::mutate(totalChange = sum(occ_2016) - sum(occ_1990)) %>%
      ungroup()
      
    e <- 0.0001
    
      #take section with increase
      changeIncrease <- changeTemp %>% 
                          filter(totalChange>0) %>%
                          dplyr::group_by(sim_nu, Core) %>%
        dplyr::mutate(free_2016 = total_2016 - occ_2016,
                      free_1990 = total_1990 - occ_1990) %>%
        dplyr::group_by(sim_nu) %>%
        dplyr::mutate(free_total = sum(free_1990)) %>%
         ungroup() %>%
        dplyr::mutate(freeRatio = free_1990/free_total,
                                 gain_obs = occ_2016 - occ_1990,
                                 gain_null = totalChange * freeRatio,
                                 chi = (gain_obs - gain_null)^2/(gain_null + e),
                                 direction = ifelse(gain_obs > gain_null, 1, -1),
                                 chi = chi * direction,
                                total_excess = gain_obs - gain_null) %>% # abs diff
                                dplyr::select(sim_nu, Core, 
                                              gain_obs, gain_null,
                                              chi, direction,
                                              total_excess, totalChange)  %>%
                                rename(obs_change = gain_obs, null_change = gain_null)
    
      
      #or with decreases
      changeDecrease <- changeTemp %>% 
                          filter(totalChange<0) %>%
        dplyr::group_by(sim_nu, Core) %>%
        dplyr::mutate(loss_2016 = occ_2016 - occ_1990) %>%
        dplyr::group_by(sim_nu) %>%
        dplyr::mutate(occ_total = sum(occ_1990)) %>%
                          ungroup() %>%
        dplyr::mutate(occRatio = occ_1990/occ_total,
                                  loss_obs = occ_2016 - occ_1990,
                                  loss_null = totalChange * occRatio,
                                  chi = (loss_obs - loss_null)^2/(abs(loss_null) + e),
                                  direction = ifelse(loss_obs < loss_null, -1, 1),
                                  chi = chi * direction,
                                total_excess = loss_obs-loss_null) %>% #absolute change
                        dplyr::select(sim_nu, Core, loss_obs, loss_null,
                                        chi, direction, total_excess, totalChange) %>%
                              rename(obs_change = loss_obs, null_change = loss_null)
      
      #combine both
      bind_rows(changeIncrease, changeDecrease) %>%
        add_column(Species=myspecies)
      
      
  }
  
}


summariseCoreChange <- function(allChanges, level="national"){

if(level=="national"){ #proportion difference

  allChanges %>%
    filter(Core=="core") %>%
    group_by(Species) %>%
    dplyr::summarize(totalChange = mean(totalChange),
                     medianExcess = quantile(total_excess,0.5),
                     sdExcess = sd(total_excess),
                     lowerExcess = quantile(total_excess,0.025),
                     upperExcess = quantile(total_excess,0.975)) %>%
    ungroup()
                     #medianTotChi = quantile(tot_chi,0.5),
                     #sdTotChi = sd(tot_chi),
                     #lowerTotChi = quantile(tot_chi,0.025),
                     #upperTotChi = quantile(tot_chi,0.975),
                     #medianOdds = quantile(abs(total_odds),0.5),
                     #sdOdds = sd(abs(total_odds)),
                     #lowerOdds = quantile(abs(total_odds),0.025),
                     #upperOdds = quantile(abs(total_odds),0.975)) %>%

}else if (level=="regional"){ #change in core and marginal

  allChanges %>%
    group_by(Core, Species) %>%
    dplyr::summarize(totalChange = mean(totalChange),
                     obs = mean(obs_change),
                     null = mean(null_change),
                     medianCoreChi = quantile(chi,0.5),
                     sdCoreChi = sd(chi),
                     lowerCoreChi = quantile(chi,0.025),
                     upperCoreChi = quantile(chi,0.975),
                     medianDirection = quantile(direction,0.5),
                     sdDirection = sd(direction),
                     lowerDirection = quantile(direction,0.025),
                     upperDirection = quantile(direction,0.975)) %>%
    ungroup()
}
  
}



getColExt <- function(myspecies){
  
  #add on core region info
  coreDF_species <- subset(coreDF, Species== myspecies)
  modelSummaries_Limits$Core <- coreDF_species$Core[match(modelSummaries_Limits$MTB, 
                                                          coreDF_species$MTB)]
  modelSummaries_Limits$Core[modelSummaries_Limits$Core=="absent"] <- "marginal"
  
  #colonizations
  colTemp <- lapply(1:dim(PA_matrix)[2], function(i){
    
    emptySites <- modelSummaries_Limits %>%
      add_column(PA = PA_matrix[,i]) %>%
      filter(!is.na(Core)) %>%
      filter(Species==myspecies) %>%
      filter(Year==1990 & PA==0) %>%
      pull(MTB)
    
    modelSummaries_Limits %>%
      add_column(PA = PA_matrix[,i]) %>%
      filter(!is.na(Core)) %>%
      filter(Year==2016) %>%
      filter(Species==myspecies) %>%
      filter(MTB %in% emptySites) %>%
      group_by(Core) %>%
      summarize(occ = sum(PA), total= length(PA)) %>%
      mutate(prop = occ/total) %>%
      add_column(simNu = i)
  }) %>% 
    bind_rows() %>% 
    ungroup() %>%
    add_column(Type = "Colonization")
  
  
  #extinctions
  extTemp <- lapply(1:dim(PA_matrix)[2], function(i){
    
    occuSites <- modelSummaries_Limits %>%
      add_column(PA = PA_matrix[,i]) %>%
      filter(!is.na(Core)) %>%
      filter(Species==myspecies) %>%
      filter(Year==1990 & PA==1) %>%
      pull(MTB)
    
    modelSummaries_Limits %>%
      add_column(PA = PA_matrix[,i]) %>%
      filter(!is.na(Core)) %>%
      filter(Year==2016) %>%
      filter(Species==myspecies) %>%
      filter(MTB %in% occuSites) %>%
      group_by(Core) %>%
      summarize(occ = sum(PA), total= length(PA)) %>%
      mutate(prop = 1-occ/total) %>%
      add_column(simNu = i)
  }) %>% 
    bind_rows() %>% 
    ungroup() %>%
    add_column(Type = "Extinction")
  
  bind_rows(colTemp,extTemp) %>%
    group_by(Core,Type) %>%
    summarise(medianProp = median(prop),
              sdProp = sd(prop),
              lowerProp = quantile(prop,0.025),
              upperProp = quantile(prop, 0.975)) %>%
    group_by(Type) %>%
    mutate(ratio = boot::logit(medianProp[Core=="core"])-boot::logit(medianProp[Core=="marginal"])) %>%
    ungroup() %>%
    add_column(species = myspecies)
  
}

### saturation ####

applySaturation <- function(species, modelSummaries_Limits, summary = "change"){
  
  #apply to each realization
  temp <- lapply(1:ncol(PA_matrix),function(i){
    
    modelSummaries_Limits$PA <- PA_matrix[,i] 
    
    #get hull
    out <- getConcaveMan(species, modelSummaries_Limits)
    
    #add on number of occupied grids
    out2 <- getRangeArea(species, modelSummaries_Limits) %>%
              filter(Year %in% c(1990,2016))
    out$sumGrids <- out2$sumGrids
    
    #deal with zeros - assume one grid occupied
    out$rangeMan[is.na(out$rangeMan)] <- as.numeric(meanArea)
    out$sumGrids[is.na(out$sumGrids)] <- as.numeric(meanArea)
    out$rangeMan[out$rangeMan==0] <- as.numeric(meanArea)
    out$sumGrids[out$sumGrids==0] <- as.numeric(meanArea)
    
    out$saturation <- out$sumGrids/out$rangeMan
    out$saturation[is.infinite(out$saturation)] <- 0
    out$sim <- i
    return(out[,c("Species","Year","saturation","sim")])
  })
  temp <- do.call(rbind,temp)
  
  #summarise change
  if(summary == "change"){
    temp %>%
      tidyr::pivot_wider(everything(),names_from = Year, values_from = saturation) %>%
      janitor::clean_names() %>%  
      dplyr::mutate(change = (x2016 - x1990))  %>%
      #dplyr::mutate(change = (boot::logit(x2016) - boot::logit(x1990-0.01)))  %>%
      dplyr::group_by(species) %>%
      dplyr::summarise(medianChange = median(change,na.rm=T), 
                       sdChange = sd(change,na.rm=T), 
                       lowerChange = quantile(change, 0.025,na.rm=T),
                       upperChange = quantile(change, 0.975,na.rm=T))
    
  }else if(summary == "annual"){
    
    #summarise annual
    temp %>%
      dplyr::group_by(Species,Year) %>%
      dplyr::summarise(medianArea = median(saturation),
                       sdArea = sd(saturation),
                       lowerArea = quantile(saturation, 0.025),
                       upperArea = quantile(saturation, 0.975))
    
  }
  
}



### ecoregion analysis ####

applyEcoregion <- function(species, modelSummaries_Limits, summary = "change"){
  
  #apply to each realization
  temp <- lapply(1:ncol(PA_matrix),function(i){
    
    modelSummaries_Limits$PA <- PA_matrix[,i]
    
    modelSummaries_Limits %>%
      filter(Species==species) %>%
      dplyr::group_by(Naturraum,Year) %>%
      dplyr::summarise(nu=sum(PA),total=length(PA)) %>%
      add_column(Species=species, simNu=i)
    
  })
  
  temp <- do.call(rbind,temp)
  
  
  #summarise
  if(summary=="change"){
    temp %>%
      dplyr::group_by(Naturraum,Species,simNu) %>%
      dplyr::summarize(change = (nu[Year==2016]-nu[Year==1990]+1)/(total[Year==1990])) %>%
      dplyr::group_by(Species,Naturraum) %>%
      dplyr::summarise(medianChange = quantile(change,0.5),
                       sdChange = sd(change),
                       lowerChange = quantile(change,0.025),
                       upperChange = quantile(change,0.975))
    
    
  }else if(summary=="annual"){
    temp %>%
      dplyr::group_by(Species,Naturraum,simNu) %>%
      dplyr::summarize(prop1990 = nu[Year==1990]/total[Year==1990],
                       prop2016 = nu[Year==2016]/total[Year==2016]) %>%
      dplyr::group_by(Species,Naturraum) %>%
      dplyr::summarise(median1990 = quantile(prop1990,0.5),
                       lower1990 = quantile(prop1990,0.025),
                       upper1990 = quantile(prop1990,0.975),
                       median2016 = quantile(prop2016,0.5),
                       lower2016 = quantile(prop2016,0.025),
                       upper2016 = quantile(prop2016,0.975))
    
  }
  
}

