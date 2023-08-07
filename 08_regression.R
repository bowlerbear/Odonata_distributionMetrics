#for each section, run the start scripts of 07

library(brms)

### aoo vs eoo change ####

#### statistics ######

nrow(subset(allChanges, lowerChange_area>0))
nrow(subset(allChanges, upperChange_area<0))
nrow(subset(allChanges, lowerChange_extent>0))
nrow(subset(allChanges, upperChange_extent<0))

subset(allChanges, lowerChange_area>0 & lowerChange_extent>0 & 
         medianChange_area<1.6 & medianChange_area>0)

#### regression ####

#differences in magnitude of winners and losers in extent
brm1 <- brm(abs(medianChange_extent) | mi(sdChange_extent) ~ Direction,
            data = allChanges, save_mevars = TRUE)
summary(brm1)

#differences in magnitude of winners and losers in area
brm1 <- brm(abs(medianChange_area) | mi(sdChange_area) ~ Direction,
            data = allChanges, save_mevars = TRUE)
summary(brm1)

#interaction
brm1 <- brm(medianChange_extent | mi(sdChange_extent) ~ Direction * me(medianChange_area, sdChange_area),
            data = allChanges, save_mevars = TRUE)
summary(brm1)

#### breakpoint #####

bform3 <- bf(
  medianChange_extent| mi(sdChange_extent) ~ Intercept + slope1 *  medianChange_area * step(change - medianChange_area) +  # Section 1
    (slope1 * change + slope2 * (medianChange_area - change)) * step(medianChange_area - change),  # Section 2
  Intercept + slope1 + slope2 ~ 1, # Fixed intercept and slopes
  change ~ 1,
  nl = TRUE
)

# Priors
bprior3 <- prior(normal(0, 5), nlpar = "Intercept") +
  prior(normal(0, 2), nlpar = "slope1") +
  prior(normal(0, 2), nlpar = "slope2") +
  prior(uniform(0, 10), nlpar = "change")  # Within observed range

# Initial values
inits3 = list(list(
  slope1 = 0.14,
  slope2 = 1.5,
  Intercept = 0.04
))

# Fit it!
fit3 <- brm(bform3, 
            data = allChanges, 
            prior = bprior3, 
            chains = 1, 
            init = inits3)

summary(fit3)#1.59


subset(allChanges, medianChange_area>1.59)
allChanges <- allChanges %>% arrange(desc(medianChange_area))

#test relationship within this area
brm1 <- brm(medianChange_extent | mi(sdChange_extent) ~ me(medianChange_area, sdChange_area),
            data = subset(allChanges,medianChange_area<1.59), 
            save_mevars = TRUE)
summary(brm1)

brm1 <- brm(medianChange_extent | mi(sdChange_extent) ~ me(medianChange_area, sdChange_area),
            data = subset(allChanges,medianChange_area>=1.59), 
            save_mevars = TRUE)
summary(brm1)

#get predicted extent based on first model
allChanges$predExtent <- predict(brm1,newdata = allChanges)[,1]
allChanges$predExtentupper <- predict(brm1,newdata = allChanges)[,4]
qplot(predExtent, medianChange_extent, data = allChanges) + 
  geom_abline(intecept=0, slope=1)

subset(allChanges, medianChange_area>1.59)
subset(allChanges, lowerChange_extent > predExtentupper)

### clumpines/saturation ####

#### statistics ####

#fragmentaion
nrow(subset(fragAnnualChanges, lowerChange>0))
nrow(subset(fragAnnualChanges, upperChange<0))

nrow(subset(fragAnnualChanges, Direction=="Losers"))
nrow(subset(fragAnnualChanges, upperChange<0 & Direction=="Losers"))
nrow(subset(fragAnnualChanges, Direction=="Winners"))
nrow(subset(fragAnnualChanges, upperChange<0 & Direction=="Winners"))

#saturation
nrow(subset(saturationChanges, Direction=="Winners"))
nrow(subset(saturationChanges, lowerChange>0 & Direction=="Winners"))

nrow(subset(saturationChanges, Direction=="Losers"))
nrow(subset(saturationChanges, upperChange<0 & Direction=="Losers"))

#### regression ####

#fragchanges
brm1 <- brm(medianChange | mi(sdChange) ~ 1,
            data = subset(fragAnnualChanges, Direction=="Winners"), save_mevars = TRUE)
summary(brm1)

brm1 <- brm(medianChange | mi(sdChange) ~ 1,
            data = subset(fragAnnualChanges, Direction=="Losers"), save_mevars = TRUE)
summary(brm1)

brm1 <- brm(medianChange | mi(sdChange) ~ Direction,
            data = fragAnnualChanges, save_mevars = TRUE)
summary(brm1)

#saturation
brm1 <- brm(medianChange | mi(sdChange) ~ Direction,
            data = saturationChanges, save_mevars = TRUE)

brm1 <- brm(medianChange | mi(sdChange) ~ 1,
            data = subset(saturationChanges, Direction=="Winners"), save_mevars = TRUE)
summary(brm1)

brm1 <- brm(medianChange | mi(sdChange) ~ 1,
            data = subset(saturationChanges, Direction=="Losers" & !is.na(sdChange)), save_mevars = TRUE)
summary(brm1)

subset(saturationChanges, Direction=="Winners" & medianChange <0)
subset(saturationChanges, Direction=="Winners" & upperChange >0) %>%
  arrange(desc(medianChange))

#relationship with area

# fragmention

brm1 <- brm(medianChange_frag | mi(sdChange_frag) ~ me(medianChange_area, sdChange_area),
            data = allChanges, save_mevars = TRUE)
summary(brm1)

brm1 <- brm(medianChange_frag | mi(sdChange_frag) ~ Direction_area * me(medianChange_area, sdChange_area),
            data = allChanges, save_mevars = TRUE)
summary(brm1)#ns

#saturation

brm1 <- brm(medianChange_sat | mi(sdChange_sat) ~ me(medianChange_area, sdChange_area),
            data = allChanges, save_mevars = TRUE)
summary(brm1)

brm1 <- brm(medianChange_sat | mi(sdChange_sat) ~ Direction_area * me(medianChange_area, sdChange_area),
            data = allChanges, save_mevars = TRUE)
summary(brm1)

### core-margin #####

coreAnnualChanges <- readRDS("outputs/coreAnnualChanges.rds")

#### statistics ####

coreChanges_wide <- coreAnnualChanges %>%
                      pivot_wider(names_from="Core",
                                  values_from=c("medianChange",
                                                "sdChange",
                                                "lowerChange",
                                                "upperChange"))

cor(coreChanges_wide$medianChange_core, coreChanges_wide$medianChange_marginal)

nrow(coreChanges_wide)
nrow(subset(coreChanges_wide, lowerChange_core >0 & lowerChange_marginal>0))
nrow(subset(coreChanges_wide, upperChange_core <0 & upperChange_marginal<0))
nrow(subset(coreChanges_wide, upperChange_core <0 & lowerChange_marginal>0))
nrow(subset(coreChanges_wide, lowerChange_core >0 & upperChange_marginal<0))

### latitude #####

head(latitudinalChanges)

brm1 <- brm(median | mi(sd) ~ 1, data = subsey(latitudinalChanges, limit=="Northern"), save_mevars = TRUE)
summary(brm1)

brm1 <- brm(median | mi(sd) ~ 1, data = subsey(latitudinalChanges, limit=="Southern"), save_mevars = TRUE)
summary(brm1)


head(allChanges)

brm1 <- brm(median | mi(sd) ~ me(medianChange, sdChange), data = subsey(allChanges, limit=="Northern"), save_mevars = TRUE)
summary(brm1)

brm1 <- brm(median | mi(sd) ~ me(medianChange, sdChange), data = subsey(allChanges, limit=="Southern"), save_mevars = TRUE)
summary(brm1)
