##########################################################################################
### (2) Model fitting and selection R code for:

### "Timing and intensity of weather events shape nestling development strategies 
###  in three alpine breeding songbirds"

### DR de Zwaan, A Drake, JL Greenwood, and K Martin (2020)
### Frontiers in Ecology and Evolution	

### doi: 10.3389/fevo.2020.570034

### Run with R version 3.6.3
########################################################################################## 

### Set working directory

setwd("C:/PATH")


### Load required packages

library(plyr)
library(dplyr)
library(mgcViz)
library(mgcv)
library(lme4)
library(lmerTest)
library(olsrr)
library(piecewiseSEM)
library(ggplot2)
library(tidymv)
library(doBy)
library(broom)


### Read in HOLA, DEJU, and SAVS data produced in 01_de Zwaan et al_2020_Sliding window R code

HOLA_data <- read.csv("HOLA_weather_summarized_final.csv")
DEJU_data <- read.csv("DEJU_weather_summarized_final.csv")
SAVS_data <- read.csv("SAVS_weather_summarized_final.csv")



###########################################################################################
### (1) Horned lark (HOLA)
###########################################################################################


######################
### Prepare data
######################

### Scale weather variables

head(HOLA_data)

HOLA_sc <- scale(HOLA_data[,c(14:32)], center=TRUE, scale=TRUE)

### Recombine with non-standardized variables

HOLA_sc <- as.data.frame(HOLA_sc)

HOLA_sc$DFE <- HOLA_data$DFE
HOLA_sc$total_storms_w <- HOLA_data$total_storms_w
HOLA_sc$total_storms_t <- HOLA_data$total_storms_t
HOLA_sc$year <- HOLA_data$year
HOLA_sc$nest_ID <- HOLA_data$nest_id
HOLA_sc$treatment <- HOLA_data$HOLA_treatment
HOLA_sc$nestling_ID <- HOLA_data$nestling_id
HOLA_sc$measurement_age <- HOLA_data$measurement_age
HOLA_sc$brood_size <- HOLA_data$brood_size
HOLA_sc$wing <- HOLA_data$wing
HOLA_sc$tarsus <- HOLA_data$tarsus
HOLA_sc$mass <- HOLA_data$mass


### Convert nest id into factor

HOLA_sc$nest_ID <- factor(HOLA_sc$nest_ID)


### Convert storm events to a factor

HOLA_sc$storms_w_f <- factor(HOLA_sc$total_storms_w)

HOLA_sc$storms_t_f <- factor(HOLA_sc$total_storms_t)



### Test correlation between avg temp and sub 10 hrs

cor.test(HOLA_sc$avg_temp_m, HOLA_sc$sub_10_m) # -0.93

### Note: correlation so high that we will continue with only average temperature, not sub 10 hrs.



##############################
### (A) HOLA: wing model
##############################

### Note:
# See model selection methods for details. Terms penalized out of the model are removed
# by commenting out. 
# The number behind certain variables indicates the model it was included in.

# If all terms are considered linear, then a linear model is fit below.

##############
### HOLA wing GAM
##############

HOLA_wing <- gam(wing ~ measurement_age + brood_size + treatment + DFE 
                        + s(avg_temp_w, k=4)    
                        #+ s(min_temp_w, k=4) ### model 2
                        #+ s(max_temp_w, k=4) ### model 3
                        + s(sub_5_short_w, k=4)
                        #+ s(sub_5_long_w, k=4)
                        #+ s(total_early_precip_w, k=4)
                        #+ s(total_late_precip_w, k=4)
                        + storms_w_f
                        + s(nest_ID, bs = "re"),
                           data= HOLA_sc, select=TRUE, method = 'REML')

summary(HOLA_wing) ## 63.1% deviance explained

AIC(HOLA_wing, HOLA_wing2, HOLA_wing3)

### Note: all temp variables fit equally well in GAM format

### Note: all terms are linear so fit a linear model.


### Linear wing model

HOLA_wing_lm <- lmer(wing ~ measurement_age + brood_size + treatment + DFE 
                     + storms_w_f 
                     + avg_temp_w 
                     #+ min_temp_w   ###### model 2
                     #+ max_temp_w   ###### model 3
                     + sub_5_short_w 
                     + (1|nest_ID), data=HOLA_sc)

AIC(HOLA_wing_lm, HOLA_wing_lm2, HOLA_wing_lm3) # Mean model is the best.

# Results
summary(HOLA_wing_lm)
confint(HOLA_wing_lm)
rsquared(HOLA_wing_lm) # 56.4% R2



##########
### Linear model check: wing
##########

### For linear mixed effect model:

## Standardized residuals versus fitted values
plot(HOLA_wing_lm, resid(., scaled=TRUE) ~ fitted(.), abline = 0)

### Observed versus fitted values
plot(HOLA_wing_lm, wing ~ fitted(.), abline = c(0,1))

### Test for colinearity among weather variables

HOLA_wing_linear <- lm(wing ~ DFE + brood_size + measurement_age 
                       + avg_temp_w + sub_5_short_w + storms_w_f, data=HOLA_sc)

### Colinearity

ols_coll_diag(HOLA_wing_linear)

# Low collinearity - all weather variables less than 3.



##############################
### (B) HOLA: tarsus model
##############################

################
### GAM
################

HOLA_tarsus <- gam(tarsus ~ measurement_age + brood_size + treatment + DFE + storms_t_f 
                 #+ s(avg_temp_t, k=4)  
                 #+ s(min_temp_t, k=4) #### model 2
                 #+ s(max_temp_t, k=4) #### model 3
                 + s(sub_5_short_t, k=4)
                 #+ s(sub_5_long_t, k=4)
                 + s(nest_ID, bs = "re"),
                 data= HOLA_sc, select=TRUE, method = 'REML')

summary(HOLA_tarsus) ## 64.3% deviance explained

AIC(HOLA_tarsus, HOLA_tarsus2, HOLA_tarsus3) # No difference between temperature GAM models


### Note: 
# No non-linearities so try linear model.


### Linear model

HOLA_tarsus_lm <- lmer(tarsus ~ measurement_age + brood_size + treatment + DFE
                   + sub_5_short_t 
                   + storms_t_f
                   + (1|nest_ID), data= HOLA_sc)

summary(HOLA_tarsus_lm)
confint(HOLA_tarsus_lm)
rsquared(HOLA_tarsus_lm) # R2 57.6


##########
### Linear model check: tarsus
##########


### Standardized residuals versus fitted values
plot(HOLA_tarsus_lm, resid(., scaled=TRUE) ~ fitted(.), abline = 0)

### Observed versus fitted values
plot(HOLA_tarsus_lm, tarsus ~ fitted(.), abline = c(0,1))

### Test for colinearity among linear weather variables

HOLA_tarsus_linear <- lm(tarsus ~ DFE + brood_size + measurement_age
                         + storms_t_f + sub_5_short_t, data=HOLA_sc)

### Colinearity

ols_coll_diag(HOLA_tarsus_linear)

### No collinearity - all weather variables less than 3.



##############################
### (C) HOLA: mass model
##############################

##################
### GAM
##################

HOLA_mass <- gam(mass ~ measurement_age + brood_size + treatment + DFE
                  #+ s(avg_temp_m, k=4)
                  + s(min_temp_m, k=4)
                  + s(sub_5_m, k=4)
                  #+ s(total_precip_m, k=4) ##### Model 2
                  + s(nest_ID, bs = "re"),
                  data= HOLA_sc, select=TRUE, method = 'REML')

summary(HOLA_mass) ## 63.9% deviance explained

AIC(HOLA_mass, HOLA_mass2) 

### Note: No difference if remove precipitation, so remove precipitation


### Note: no non-linearities so fit linear model


### Linear model

HOLA_mass_lm <- lmer(mass ~ measurement_age + brood_size + treatment + DFE
                        + sub_5_m 
                        + min_temp_m 
                        + (1|nest_ID), data= HOLA_sc)

summary(HOLA_mass_lm)
confint(HOLA_mass_lm)
rsquared(HOLA_mass_lm) # R2 54.4



################
### Linear model check: Mass
################


## Standardized residuals versus fitted values
plot(HOLA_mass_lm, resid(., scaled=TRUE) ~ fitted(.), abline = 0)

### Observed versus fitted values
plot(HOLA_mass_lm, mass ~ fitted(.), abline = c(0,1))

### Test for colinearity in linear terms

HOLA_mass_linear <- lm(mass ~ DFE + brood_size + measurement_age 
                       + sub_5_m + avg_temp_w, data=HOLA_sc)

### Colinearity
ols_coll_diag(HOLA_mass_linear)

# No collinearity - all weather variables less than 3.



######################################################################################################
#### (2) Dark eyed Junco ####################################################################
######################################################################################################


###########
### Prepare data
###########

### Scale weather variables

head(DEJU_data)

DEJU_sc <- scale(DEJU_data[,c(13:24)], center=TRUE, scale=TRUE)

### Recombine with other variables
DEJU_sc <- as.data.frame(DEJU_sc)

DEJU_sc$DFE <- DEJU_data$DFE

DEJU_sc$total_storms_early_m <- DEJU_data$total_storms_early_m
DEJU_sc$total_storms_middle_m <- DEJU_data$total_storms_middle_m
DEJU_sc$total_storms_late_m <- DEJU_data$total_storms_late_m


DEJU_sc$nest_ID <- DEJU_data$nest_id
DEJU_sc$measurement_age <- DEJU_data$measurement_age
DEJU_sc$brood_size <- DEJU_data$brood_size
DEJU_sc$wing <- DEJU_data$wing
DEJU_sc$tarsus <- DEJU_data$tarsus
DEJU_sc$mass <- DEJU_data$mass


### Convert nest id into a factor

DEJU_sc$nest_ID <- factor(DEJU_sc$nest_ID)


### Convert storm events to a factor as an alternative

DEJU_sc$storms_early_m_f <- factor(DEJU_sc$total_storms_early_m)

DEJU_sc$storms_middle_m_f <- factor(DEJU_sc$total_storms_middle_m)

DEJU_sc$storms_late_m_f <- factor(DEJU_sc$total_storms_late_m)



##############################
### (A) DEJU wing model
##############################

########
### GAM
########

DEJU_wing <- gam(wing ~ measurement_age + brood_size + DFE 
                 #+ s(avg_temp_w, k=4) 
                 + s(max_temp_w, k=4) 
                 #+ s(sub_5_w, k=4)
                  + s(nest_ID, bs = "re"),
                 data= DEJU_sc, select=TRUE, method = 'REML')

summary(DEJU_wing) ## 68.3% deviance explained


### Linear model

DEJU_wing_lm <- lmer(wing ~ measurement_age + brood_size + DFE 
                     + max_temp_w
                     + (1|nest_ID), data=DEJU_sc)

summary(DEJU_wing_lm)
confint(DEJU_wing_lm)
rsquared(DEJU_wing_lm) # 61.8% R2


############################
### Linear model check: DEJU wing
############################

## Standardized residuals versus fitted values
plot(DEJU_wing_lm, resid(., scaled=TRUE) ~ fitted(.), abline = 0)

### Observed versus fitted values
plot(DEJU_wing_lm, wing ~ fitted(.), abline = c(0,1))


### Test for colinearity among linear terms

DEJU_wing_linear <- lm(wing ~ measurement_age + DFE + brood_size
                       + avg_temp_w, data=DEJU_sc)

### Colinearity
ols_coll_diag(DEJU_wing_linear)

# No collinearity - weather variable less than 2.





##############################
### (B) DEJU tarsus model
##############################

#########
### GAM
#########


DEJU_tarsus <- gam(tarsus ~ measurement_age + brood_size + DFE 
                 #+ s(avg_temp_early_t, k=4) 
                 #+ s(avg_temp_late_t, k=4)
                 #+ s(min_temp_t, k=4)
                 + s(max_temp_t, k=4)
                 #+ s(sub_5_t, k=4)
                 + s(nest_ID, bs = "re"),
                 data= DEJU_sc, select=TRUE, method = 'REML')

summary(DEJU_tarsus) ## 69.1% deviance explained



### Linear model

DEJU_tarsus_lm <- lmer(tarsus ~ measurement_age + brood_size + DFE 
                       + max_temp_t
                       + (1|nest_ID), data=DEJU_sc)

summary(DEJU_tarsus_lm)
confint(DEJU_tarsus_lm)
rsquared(DEJU_tarsus_lm) # 65.0% R2


#############################
### Linear model check: DEJU tarsus
#############################

## Standardized residuals versus fitted values
plot(DEJU_tarsus_lm, resid(., scaled=TRUE) ~ fitted(.), abline = 0)

### Observed versus fitted values
plot(DEJU_tarsus_lm, tarsus ~ fitted(.), abline = c(0,1))


### Test for colinearity among linear terms

DEJU_tarsus_linear <- lm(tarsus ~ measurement_age + DFE + brood_size
                       + avg_temp_late_t, data=DEJU_sc)

### Colinearity
ols_coll_diag(DEJU_tarsus_linear)

# No collinearity - weather variable less than 2.



##############################
### (C) DEJU mass model
##############################

#########
### GAM
#########

DEJU_mass <- gam(mass ~ measurement_age + brood_size + DFE 
                   + s(max_temp_m, k=3)
                   #+ s(total_storms_early_m, k=3)
                   #+ s(total_storms_middle_m, k=3)
                   #+ s(total_storms_late_m, k=3)
                   + s(nest_ID, bs = "re"),
                   data= DEJU_sc, select=TRUE, method = 'REML')

summary(DEJU_mass) ## 67.2% deviance explained


### Note: All terms are linear so fit linear model.


### Linear model

DEJU_mass_lm <- lmer(mass ~ measurement_age + brood_size + DFE 
                       + max_temp_m
                      #+ storms_late_m_f
                       + (1|nest_ID), data=DEJU_sc)

AIC(DEJU_mass_lm, DEJU_mass_lm2) 

433.7800 - 431.9765 # 1.8 AIC

summary(DEJU_mass_lm)
summary(DEJU_mass_lm2)

confint(DEJU_mass_lm)
rsquared(DEJU_mass_lm) # 61.9% R2
### Note: Model including storms similar to top model, but not significant.


#################
### Linear model check: DEJU mass
#################


## Standardized residuals versus fitted values
plot(DEJU_mass_lm, resid(., scaled=TRUE) ~ fitted(.), abline = 0)

### Observed versus fitted values
plot(DEJU_mass_lm, mass ~ fitted(.), abline = c(0,1))


### Test for colinearity among linear terms

DEJU_mass_linear <- lm(mass ~ measurement_age + DFE + brood_size
                         + storms_late_m_f, data=DEJU_sc)

### Colinearity
ols_coll_diag(DEJU_mass_linear)

# No collinearity - weather variable less than 1.



###########################################################################################
### (3) Savannah sparrow
###########################################################################################


###########
### Prepare data
###########


### Scale weather variables
head(SAVS_data)

SAVS_sc <- scale(SAVS_data[,c(11:13,21:27)], center=TRUE, scale=TRUE)


### Recombine with other variables.
SAVS_sc <- as.data.frame(SAVS_sc)

SAVS_sc$DFE <- SAVS_data$DFE

SAVS_sc$nest_ID <- SAVS_data$nest_id
SAVS_sc$measurement_age <- SAVS_data$measurement_age
SAVS_sc$brood_size <- SAVS_data$brood_size
SAVS_sc$tarsus <- SAVS_data$tarsus
SAVS_sc$mass <- SAVS_data$mass

### Convert nest id into a factor

SAVS_sc$nest_ID <- factor(SAVS_sc$nest_ID)



##############################
### (A) SAVS tarsus model
##############################

#############
### GAM
#############

SAVS_tarsus <- gam(tarsus ~ measurement_age + brood_size + DFE 
                 + s(avg_temp_t, k=4)
                 #+ s(min_temp_t, k=4)   
                 + s(nest_ID, bs = "re"),
                 data = SAVS_sc, select=TRUE, method = 'REML')

summary(SAVS_tarsus) ## 61.8% deviance explained


### Note: All terms are linear so fit linear model.


### Linear model

SAVS_tarsus_lm <- lmer(tarsus ~ measurement_age + brood_size + DFE 
                     + avg_temp_t
                     + (1|nest_ID), data=SAVS_sc)

summary(SAVS_tarsus_lm)
confint(SAVS_tarsus_lm)
rsquared(SAVS_tarsus_lm) # 55.8% R2



####################
### Linear model check: SAVS tarsus
####################

## Standardized residuals versus fitted values
plot(SAVS_tarsus_lm, resid(., scaled=TRUE) ~ fitted(.), abline = 0)

### Observed versus fitted values
plot(SAVS_tarsus_lm, tarsus ~ fitted(.), abline = c(0,1))

### Test for colinearity among linear terms

SAVS_tarsus_linear <- lm(tarsus ~ measurement_age + DFE + brood_size
                       + avg_temp_t, data=SAVS_sc)

### Colinearity
ols_coll_diag(SAVS_tarsus_linear)

# No collinearity - weather variable less than 2.



##############################
### (B) SAVS mass model
##############################

### Subset to just nestlings with a mass measurement

SAVS_mass_sc <- subset(SAVS_sc, mass >0)

#############
### GAM
#############

SAVS_mass <- gam(mass ~ measurement_age + brood_size + DFE 
                   #+ s(avg_temp_early_m, k=4)
                   #+ s(avg_temp_late_m, k=4)
                   #+ s(min_temp_m, k=4)
                   #+ s(max_temp_early_m, k=4)
                   #+ s(max_temp_late_m, k=4)
                   #+ s(sub_5_m, k=4)
                   + s(total_precip_m, k=4)
                   + s(nest_ID, bs = "re"),
                   data = SAVS_mass_sc, select=TRUE, method = 'REML')

summary(SAVS_mass) ## 56.4% deviance explained

### Note: If precipitation is not included, then avg temp is the strongest factor.


### Linear model

SAVS_mass_lm2 <- lmer(mass ~ measurement_age + brood_size + DFE 
                        + avg_temp_early_m  #### Model 2
                        #+ total_precip_m
                        + (1|nest_ID), data=SAVS_mass_sc)

AIC(SAVS_mass_lm, SAVS_mass_lm2) ### Precipitation model is best.

summary(SAVS_mass_lm)
confint(SAVS_mass_lm)
rsquared(SAVS_mass_lm) # 50.8% R2

summary(SAVS_mass_lm2)

### Test correlation between precip and temperature.

cor.test(SAVS_mass_sc$total_precip_m, SAVS_mass_sc$avg_temp_early_m) # 0.72

### Note: best model includes precipitation, not temperature.

### Note: Both are highly correlated (0.72), so although precipitation is selected, 
# it seems like there is probably a contribution from both and it is not
# possible to parse here.


###################
### Linear model check: SAVS mass
###################

## Standardized residuals versus fitted values
plot(SAVS_mass_lm, resid(., scaled=TRUE) ~ fitted(.), abline = 0)

### Observed versus fitted values
plot(SAVS_mass_lm, mass ~ fitted(.), abline = c(0,1))

### Test for colinearity among linear terms

SAVS_mass_linear <- lm(mass ~ measurement_age + DFE + brood_size
                         + total_precip_m 
                       , data=SAVS_mass_sc)

### Colinearity
ols_coll_diag(SAVS_mass_linear)

# No collinearity - weather variable less than 1.



############################################################################################
###################### End of model fitting and selection code #############################
############################################################################################