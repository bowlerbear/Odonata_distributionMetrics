library(tidyverse)
library(rstan)
library(brms)
library(mgcv)

#HPC
myfolder <- "/data/idiv_ess/Odonata"

#Local
#myfolder <- "splines/prep"

# load in the full dataset ---------------------------------

adultData <- readRDS(paste(myfolder,"adultData.rds",sep="/"))

mtbsDF_extended <- readRDS(paste(myfolder,"MTB_stateextendedpoints.rds",sep="/"))
mtbsDF_extended <- subset(mtbsDF_extended, type!="extension")
names(mtbsDF_extended)[1] <- "MTB"
mtbsDF <- subset(mtbsDF_extended, !is.na(MTB))
mtbsDF <- subset(mtbsDF, !duplicated(MTB))

#and environ data
environData <- readRDS(paste(myfolder, "environData.rds", sep="/")) %>%
                rename(MTB = Value) %>% mutate(MTB = as.numeric(MTB))
#mtbsDF$MTB[!mtbsDF$MTB %in% environData$Value]#all present

mtbsDF <- inner_join(mtbsDF, environData, by="MTB")

# subset years/months --------------------------------------------

adultData <- subset(adultData, Year>=1985  & Year<2018)
adultData <- subset(adultData, Month %in% 4:10)

#for testing, take one week
#adultData <- subset(adultData, week == 39)

### choose task settings ####

task.id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

#models
models <- 39

#species
allSpecies <- sort(unique(adultData$Species))
allSpecies <- allSpecies[allSpecies!="Onychogomphus uncatus"]
allSpecies <- allSpecies[allSpecies!="Coenagrion hylas"]

#all combinations
taskDF <- expand.grid(allSpecies,models)
nrow(taskDF)#79

#focal for this task
myspecies <- as.character(taskDF$Var1[task.id])
#myspecies <- "Sympetrum danae" #test case
mymodel <- taskDF$Var2[task.id]


# define a visit -------------------------------------------

adultData$visit <- paste(adultData$MTB, adultData$Date, sep="_")

### subset MTBs #####

adultData$MTB <- sapply(as.character(adultData$MTB_Q),function(x){
  len <- nchar(x)
  substr(x,1,(len-1))})

#check all are official
unique(adultData$MTB[!adultData$MTB %in% mtbsDF$MTB])
adultData <- subset(adultData, MTB %in% mtbsDF$MTB)

#add on state info
adultData$State <- mtbsDF$State[match(adultData$MTB, mtbsDF$MTB)]

# organise species and visits ------------------------------------

#get occurence matrix  - detection and non-detection
occMatrix <- reshape2::acast(adultData,visit~Species,value.var="Anzahl_min",fun=function(x)length(x[x!=0]))
occMatrix[occMatrix>0] <- 1

#get list length
visit_df <- adultData %>% 
    dplyr::group_by(visit, Date, MTB) %>%
    dplyr::summarise(nuSpecies=length(unique(Species)),
                    nuRecords=length(Species))


#add on some indices
visit_df$Date <- as.Date(visit_df$Date)
visit_df$Year <- lubridate::year(visit_df$Date)
visit_df$yday <- lubridate::yday(visit_df$Date)
visit_df$yearIndex <- as.numeric(factor(visit_df$Year))

#get other effort variables
visit_df$singleList <- ifelse(visit_df$nuSpecies==1,1,0)
visit_df$shortList <- ifelse(visit_df$nuSpecies%in%2:3,1,0)
visit_df$longList <- ifelse(visit_df$nuSpecies>3,1,0)
visit_df$LL <- ifelse(visit_df$nuSpecies==1,"single",
                      ifelse(visit_df$nuSpecies%in%2:3,"short","long"))
visit_df$NR <- mtbsDF$State[match(visit_df$MTB, mtbsDF$MTB)]

#order data
visit_df <- arrange(visit_df, visit)

# get data for a species --------------------------------------

#check all aligns
all(row.names(occMatrix)==visit_df$visit)
visit_df$Species <- occMatrix[,myspecies]
table(visit_df$Species)

# arrange data -------------------------------------------------

visit_df$siteIndex <- as.numeric(as.factor(paste0(visit_df$MTB,visit_df$Year)))
visit_df <- arrange(visit_df, siteIndex, visit)

# mtb summary --------------------------------------------------

site_df <- visit_df %>%
  group_by(siteIndex, MTB, Year) %>%
  summarise(Species = max(Species, na.rm=T))

table(site_df$Species)

#few checks
nrow(site_df)==length(unique(visit_df$siteIndex))#should be TRUE

# expand to fill mtbq grid across germany ------------------------

full_df <- mtbsDF
nuSites <- nrow(full_df)
full_df <- full_df %>% slice(rep(1:n(), each = length(unique(visit_df$Year))))
full_df$Year <- rep(sort(unique(visit_df$Year)), nuSites)
full_df$yearIndex <- rep(1:length(unique(visit_df$Year)), nuSites)

#sort species data
full_df$Species <- full_df$SpeciesOrig  <- site_df$Species[match(interaction(full_df$MTB,full_df$Year),
                                                                           interaction(site_df$MTB,site_df$Year))]
#add in temporary zeros
full_df$Species[is.na(full_df$Species)] <- 0
table(full_df$Species)

#keep original coordinates
full_df$x_MTB <- full_df$x
full_df$y_MTB <- full_df$y

#arrange file
full_df$siteIndex <- as.numeric(as.factor(paste0(full_df$MTB,full_df$Year)))
full_df <- arrange(full_df, siteIndex)
#saveRDS(siteInfo_NAs, file="splines/prep/siteInfo_environData_NAs.rds")

# fit model to get spline set up ---------------------------------

model_data_complete <- model_data <- make_standata(bf(Species ~ t2(x, y, yearIndex, 
                                                             d = c(2,1), 
                                                             bs=c("ds","cr"),
                                                             k = c(23,9))),
                                             data = full_df,
                                             family = bernoulli())

names(model_data_complete) <- sapply(names(model_data_complete), 
                                 function(x) paste0("complete_",x))

# subset code to sites with data ------------------------------------------

#define sites with data
presentData <- !is.na(full_df$SpeciesOrig)
table(presentData)

#subset remaining data
model_data$Y <- model_data_complete$complete_Y[presentData] #response
model_data$X <- model_data_complete$complete_X[presentData,] #linear predictor - only an intercept atm
model_data$Xs <- model_data_complete$complete_Xs[presentData,] #spline stuff
model_data$Zs_1_1 <- model_data_complete$complete_Zs_1_1[presentData,]
model_data$Zs_1_2 <- model_data_complete$complete_Zs_1_2[presentData,]
model_data$Zs_1_3 <- model_data_complete$complete_Zs_1_3[presentData,]

#covariate data
covMatrix_complete <- full_df[,c("meanTemp", "totalPrecip")] %>%
                      mutate(across(everything(), scale)) %>%
                      add_column(Intercept = 1) %>%
                      as.matrix()

covMatrix <- covMatrix_complete[presentData,]

# Check alignment ------------------------------------------

obs_df <- full_df[presentData,]
df_visited <- subset(visit_df, !duplicated(siteIndex))
all(obs_df$MTB==df_visited$MTB)
all(obs_df$Year==df_visited$Year)

# Occupancy states ------------------------------------------

# define a design matrix for site-level occupancy
n_site <- length(unique(visit_df$siteIndex))

#lets fit intercept only model (beyond the spline) 
#X_psi <- matrix(c(rep(1, n_site)))

#add on climate and land cover covariates
X_psi_complete <- covMatrix_complete
X_psi <- covMatrix

#total number of covariates
m_psi <- ncol(X_psi)    # m_psi is the number of columns in the site level design matrix

# Survey data --------------------------

# determine number of surveys per site
n_survey <- as.numeric(tapply(visit_df$visit,visit_df$siteIndex,length))
total_surveys <- sum(n_survey)
total_surveys

# define a survey-level design matrix for detection probabilities
#intercept only
#X_p <- matrix(c(rep(1, total_surveys)))
#m_p <- ncol(X_p)

#or covariate model
#effort terms and detection trend term
X_p_new <- data.frame(model.matrix(~LL, visit_df))
#phenology/year terms
#X_p_new$yearIndex <- as.numeric(visit_df$yearIndex)
X_p_new$yday <- as.numeric(scale(visit_df$yday))
X_p_new$yday2 <- as.numeric(scale(visit_df$yday^2))

#combine all
X_p <- as.matrix(cbind(X_p_new)) 
m_p <- ncol(X_p)#12

# indices ---------------------------------------------------------

# get start and end indices to extract slices of y for each site
start_idx <- rep(0, n_site)
end_idx <- rep(0, n_site)
for (i in 1:n_site) {
  if (n_survey[i] > 0) {
    site_indices <- which(visit_df$siteIndex == i)
    start_idx[i] <- site_indices[1]
    end_idx[i] <- site_indices[n_survey[i]]
  }
}

# create vector of whether any positive observations were seen at each site
any_seen <- rep(0, n_site)
total_seen <- rep(0, n_site)
for (i in 1:n_site) {
  if (n_survey[i] > 0) {
    any_seen[i] <- max(visit_df$Species[start_idx[i]:end_idx[i]])
    total_seen[i] <- sum(visit_df$Species[start_idx[i]:end_idx[i]])
  }
}

n_year <- length(unique(visit_df$yearIndex))
total_seen_year <- rep(0, n_year)
for (i in 1:n_year) {
    total_seen_year[i] <- sum(visit_df$Species[visit_df$yearIndex==i])
}

summary(any_seen)
summary(n_survey)
summary(start_idx)
summary(end_idx)

# as above but per mtb ----------------------------------------------------

visit_df$mtbIndex <- as.numeric(as.factor(visit_df$MTB))
n_mtbsurvey <- as.numeric(tapply(visit_df$visit,visit_df$mtbIndex,length))
summary(n_mtbsurvey)
total_mtbsurveys <- sum(n_mtbsurvey)

n_mtb = length(unique(visit_df$mtbIndex))
start_midx <- rep(0, n_mtb)
end_midx <- rep(0, n_mtb)

for (i in 1:n_mtb) {
  if (n_mtbsurvey[i] > 0) {
    site_indices <- which(visit_df$mtbIndex == i)
    start_midx[i] <- site_indices[1]
    end_midx[i] <- site_indices[n_mtbsurvey[i]]
  }
}

# create vector of whether any positive observations were seen at each site
any_seen_mtb <- rep(0, n_mtb)
total_seen_mtb <- rep(0, n_mtb)
for (i in 1:n_mtb) {
  if (n_mtbsurvey[i] > 0) {
    any_seen_mtb[i] <- max(visit_df$Species[start_midx[i]:end_midx[i]])
    total_seen_mtb[i] <- sum(visit_df$Species[start_midx[i]:end_midx[i]])
  }
}

summary(total_seen_mtb)

# Bundle data for Stan ----------------------------------------------------

stan_d <- list(n_site = n_site, 
               m_psi = m_psi, #number of parameters
               X_psi = X_psi, #design matrix for occupancy 
               X_psi_complete = X_psi_complete, #design matrix for occupancy 
               total_surveys = nrow(visit_df), 
               m_p = m_p, #number of parameters
               X_p = X_p, #design matrix for detection
               site = visit_df$siteIndex, #site indices
               y = as.numeric(visit_df$Species), #observation data
               start_idx = start_idx, #start indices
               end_idx = end_idx, # end indices
               start_midx = start_midx, #start indices
               end_midx = end_midx, # end indices
               any_seen = any_seen, #anything seen at a site
               total_seen = total_seen,
               total_seen_mtb = total_seen_mtb,
               detYear = visit_df$yearIndex,
               detState = as.numeric(as.factor(visit_df$NR)), 
               n_survey = n_survey,
               n_year = n_year,
               n_mtb = n_mtb,
               n_state = length(unique(visit_df$NR)),
               year = visit_df$yearIndex) 

# add on spline elements to the data list
stan_d <- c(stan_d, model_data)

#add on spline elements for prediction
stan_d <- c(stan_d, model_data_complete)

# Fit model ---------------------------------------------------------------

modelfile<- "complete_space_time_v20e.stan"

#select model
m_init <- stan_model(paste(myfolder,modelfile,sep="/"))

#get cores
# try to get SLURM_CPUS_PER_TASK from submit script, otherwise fall back to 1
cpus_per_task = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
rstan_options(auto_write = FALSE)
options(mc.cores = cpus_per_task)

m_fit <- sampling(m_init, 
                  chains = 4,
                  #warmup = 1000,
                  iter = 1500,
                  thin = 2,
                  data = stan_d, 
                  pars = c('mean_p','beta_p','psi_comp_rep'))

saveRDS(as.data.frame(summary(m_fit)$summary),
        file=paste0("m_fit_summary_spacetime_",myspecies,"_",mymodel,".rds"))

# Draws ---------------------------------------------------

mydraws <- extract(m_fit, permuted = TRUE, inc_warmup = FALSE)$psi_comp_rep
ndraws <- dim(mydraws)[1]
saveRDS(mydraws[(ndraws-500):ndraws,],
        file=paste0("draws_psi_spacetime_",myspecies,"_",mymodel,".rds"))
