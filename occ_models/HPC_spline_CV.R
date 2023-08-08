#code to work out number of knots for each species
library(tidyverse)
library(pROC)
library(PresenceAbsence)
library(dismo)
library(ROCR)
library(rstan)
library(brms)
library(mgcv)


#HPC
myfolder <- "/data/idiv_ess/Odonata"

#Local
#myfolder <- "splines/prep"

# load in the full dataset ------------------------------------

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

# subset years/months ----------------------------------------------

adultData <- subset(adultData, Year>=1985  & Year<2018)
adultData <- subset(adultData, Month %in% 4:10)

#for testing, take one week
#adultData <- subset(adultData, week == 39)

# choose task settings ---------------------------------------------

task.id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

#folds 1:5
folds <- 1:5

#models
models <- 39

#species
allSpecies <- sort(unique(adultData$Species))
allSpecies <- allSpecies[allSpecies!="Onychogomphus uncatus"]
allSpecies <- allSpecies[allSpecies!="Coenagrion hylas"]

#all combinations
taskDF <- expand.grid(allSpecies, models, folds)
nrow(taskDF)#395

#focal for this task
myspecies <- as.character(taskDF$Var1[task.id])
#myspecies <- "Sympetrum danae" #test case
mymodel <- taskDF$Var2[task.id]

### define a visit ####

adultData$visit <- paste(adultData$MTB, adultData$Date, sep="_")

# subset MTBs --------------------------------------------------

adultData$MTB <- sapply(as.character(adultData$MTB_Q),function(x){
  len <- nchar(x)
  substr(x,1,(len-1))})

#check all are official
unique(adultData$MTB[!adultData$MTB %in% mtbsDF$MTB])
adultData <- subset(adultData, MTB %in% mtbsDF$MTB)

#add on state info
adultData$State <- mtbsDF$State[match(adultData$MTB, mtbsDF$MTB)]

# organise species and visits -----------------------------------

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
visit_df$stateIndex <- as.numeric(factor(visit_df$NR))

#order data
visit_df <- arrange(visit_df, visit)

# get data for a species --------------------------------------

#check all aligns
all(row.names(occMatrix)==visit_df$visit)
visit_df$Species <- occMatrix[,myspecies]
table(visit_df$Species)

# arrange data ------------------------------------------

visit_df$siteIndex <- as.numeric(as.factor(paste0(visit_df$MTB,visit_df$Year)))
visit_df <- arrange(visit_df, siteIndex, visit)

# environ data -------------------------------------------

full_df <- visit_df %>%
              mutate(MTB = as.numeric(MTB)) %>%
              inner_join(., mtbsDF, by="MTB")

# site data ---------------------------------------------

site_df <- full_df %>%
              group_by(siteIndex) %>%
              mutate(Species = max(Species)) %>%
              ungroup() %>%
              filter(!duplicated(siteIndex))
              
# fit model ---------------------------------------------

model_data_complete <- model_data <- make_standata(bf(Species ~ t2(x, y, yearIndex, 
                                                                     d = c(2,1), 
                                                                     bs=c("ds","cr"),
                                                                     k = c(23,9))),
                                                     data = site_df,
                                                     family = bernoulli())


names(model_data_complete) <- sapply(names(model_data_complete), 
                                 function(x) paste0("complete_",x))


# get folds ---------------------------------

fold_df <- readRDS(paste(myfolder, "species_MTB_folds.rds", sep="/")) %>%
              filter(Species == myspecies)

#pick fold for this task
myfold <- as.character(taskDF$Var3[task.id])

trainData <- !site_df$MTB %in% fold_df$MTB[fold_df$folds==myfold]
testData <- site_df$MTB %in% fold_df$MTB[fold_df$folds==myfold]

trainvisitData <- !full_df$MTB %in% fold_df$MTB[fold_df$folds==myfold]
testvisitData <- full_df$MTB %in% fold_df$MTB[fold_df$folds==myfold]

# subset code to sites with training data -------------------------

#training dataset
model_data$y <- model_data$Y[trainData] #response
model_data$X <- model_data$X[trainData,] #linear predictor - only an intercept atm
model_data$Xs <- model_data$Xs[trainData,] #spline stuff
model_data$Zs_1_1 <- model_data$Zs_1_1[trainData,]
model_data$Zs_1_2 <- model_data$Zs_1_2[trainData,]
model_data$Zs_1_3 <- model_data$Zs_1_3[trainData,]

#testing dataset
model_data_complete$complete_y <- model_data_complete$complete_Y[testData] #response
model_data_complete$complete_X<- model_data_complete$complete_X[testData,] #linear predictor - only an intercept atm
model_data_complete$complete_Xs <- model_data_complete$complete_Xs[testData,] #spline stuff
model_data_complete$complete_Zs_1_1 <- model_data_complete$complete_Zs_1_1[testData,]
model_data_complete$complete_Zs_1_2 <- model_data_complete$complete_Zs_1_2[testData,]
model_data_complete$complete_Zs_1_3 <- model_data_complete$complete_Zs_1_3[testData,]

# Occupancy states ------------------------------------------

# define a design matrix for site-level occupancy
train_df <- site_df[trainData,]
train_df$siteIndex <- as.numeric(as.factor(train_df$siteIndex))
n_site <- nrow(train_df)

test_df <- site_df[testData,]
test_df$siteIndex <- as.numeric(as.factor(test_df$siteIndex))
n_site_complete <- nrow(test_df)

#lets fit intercept only model (beyond the spline) 
#X_psi <- matrix(c(rep(1, n_site)))

#add on climate and land cover covariates
covMatrix_complete <- site_df[,c("meanTemp", "totalPrecip")] %>%
                          mutate(across(everything(), scale)) %>%
                          add_column(Intercept = 1) %>%
                          as.matrix()

X_psi_complete <- covMatrix_complete[testData,]
X_psi <- covMatrix_complete[trainData,]

#total number of covariates
m_psi <- ncol(X_psi)    # m_psi is the number of columns in the site level design matrix

# Detection covariates -----------------

trainvisit_df <- full_df[trainvisitData,]
testvisit_df <- full_df[testvisitData,]
trainvisit_df$siteIndex <- as.numeric(as.factor(trainvisit_df$siteIndex))
testvisit_df$siteIndex <- as.numeric(as.factor(testvisit_df$siteIndex))

# define a survey-level design matrix for detection probabilities
#intercept only
#X_p <- matrix(c(rep(1, total_surveys)))
#m_p <- ncol(X_p)

#or covariate model
#effort terms and detection trend term
X_p_new <- data.frame(model.matrix(~LL, full_df))

#phenology/year terms
#X_p_new$yearIndex <- as.numeric(full_df$yearIndex)
X_p_new$yday <- as.numeric(scale(full_df$yday))
X_p_new$yday2 <- as.numeric(scale(full_df$yday^2))

#subset to training and test datasets
X_p_complete <- X_p_new[testvisitData,]
X_p <- X_p_new[trainvisitData,]

#combine all
m_p <- ncol(X_p)

# Surveys  --------------------------

# determine number of surveys per site
n_survey <- as.numeric(tapply(trainvisit_df$visit,trainvisit_df$siteIndex,length))
total_surveys <- sum(n_survey)
total_surveys

total_surveys_complete <- length(unique(testvisit_df$visit))

# indices ---------------------------------------------------------

# get start and end indices to extract slices of y for each site
start_idx <- rep(0, n_site)
end_idx <- rep(0, n_site)
for (i in 1:n_site) {
  if (n_survey[i] > 0) {
    site_indices <- which(trainvisit_df$siteIndex == i)
    start_idx[i] <- site_indices[1]
    end_idx[i] <- site_indices[n_survey[i]]
  }
}

# create vector of whether any positive observations were seen at each site
any_seen <- rep(0, n_site)
total_seen <- rep(0, n_site)
for (i in 1:n_site) {
  if (n_survey[i] > 0) {
    any_seen[i] <- max(trainvisit_df$Species[start_idx[i]:end_idx[i]])
    total_seen[i] <- sum(trainvisit_df$Species[start_idx[i]:end_idx[i]])
  }
}

summary(any_seen)
summary(n_survey)
summary(start_idx)
summary(end_idx)

# Bundle data for Stan ----------------------------------------------------

stan_d <- list(n_site = n_site, 
               n_site_complete = n_site_complete, 
               n_year = length(unique(visit_df$yearIndex)),
               n_state = length(unique(visit_df$stateIndex)),
               m_psi = m_psi, #number of parameters
               X_psi = X_psi, #design matrix for occupancy 
               X_psi_complete = X_psi_complete, #design matrix for occupancy 
               total_surveys = nrow(trainvisit_df), 
               total_surveys_complete = nrow(testvisit_df), 
               y = trainvisit_df$Species,
               m_p = m_p, #number of parameters
               X_p = X_p, #design matrix for detection
               X_p_complete = X_p_complete, #design matrix for detection
               site = trainvisit_df$siteIndex, #site indices
               site_complete = testvisit_df$siteIndex, #site indices
               detYear = trainvisit_df$yearIndex,
               detYear_complete = testvisit_df$yearIndex,
               detState = trainvisit_df$stateIndex,
               detState_complete = testvisit_df$stateIndex,
               start_idx = start_idx, #start indices
               end_idx = end_idx, # end indices
               any_seen = any_seen, #anything seen at a site
               total_seen = total_seen,
               n_survey = n_survey) 

# add on spline elements to the data list
stan_d <- c(stan_d, model_data)

#add on spline elements for prediction
stan_d <- c(stan_d, model_data_complete)

# Fit model ---------------------------------------------------------------

modelfile<- "complete_space_time_v20e_CV.stan"

#select model
m_init <- stan_model(paste(myfolder,modelfile,sep="/"))

#get cores
# try to get SLURM_CPUS_PER_TASK from submit script, otherwise fall back to 1
cpus_per_task = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
rstan_options(auto_write = FALSE)
options(mc.cores = cpus_per_task)

m_fit <- sampling(m_init, 
                  chains = 4,
                  #iter = 50,
                  iter = 1500,
                  thin = 2,
                  data = stan_d, 
                  pars = c('y_pred'))

# quick calc -------------------------------------------------

#add predictions to the data frame of observations
temp <- summary(m_fit)$summary
temp <- subset(temp, grepl("y_pred", row.names(temp)))
testvisit_df$Preds <- temp[,1]

TSS = function(thresh.dat.opt){
  PresenceAbsence::sensitivity(cmx, st.dev=F) + 
  PresenceAbsence::specificity(cmx, st.dev=F) - 1
}

evalSDM <- function(observation, predictions, 
                    thresh.method='MaxSens+Spec', 
                    req.sens=0.85, req.spec = 0.85, FPC=1, FNC=1){
  
  thresh.dat <- data.frame(ID=seq_len(length(observation)), 
                           obs = observation,
                           pred = predictions)
  
  thresh <- PresenceAbsence::optimal.thresholds(DATA = thresh.dat, 
                                                req.sens = req.sens, 
                                                req.spec = req.spec, 
                                                FPC=FPC, 
                                                FNC=FNC)
  
  thresh.dat.opt <- PresenceAbsence::cmx(DATA= thresh.dat, 
                                  threshold=thresh[thresh$Method==thresh.method,2])
  
  data.frame(AUC = PresenceAbsence::auc(thresh.dat, st.dev=F),
             TSS = TSS(cmx.opt))
}

out <- evalSDM(testvisit_df$Species, testvisit_df$Preds, req.sens=0.85, req.spec = 0.85)

saveRDS(out, file=paste0("CVmean_",myspecies,"_",mymodel,"_",myfold,".rds"))
