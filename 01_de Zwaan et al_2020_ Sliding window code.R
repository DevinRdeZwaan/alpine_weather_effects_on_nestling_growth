##########################################################################################
### (1) Sliding window analysis R code for:

### "Timing and intensity of weather events shape nestling development strategies 
###  in three alpine breeding songbirds"

###  DR de Zwaan, A Drake, JL Greenwood, and K Martin (2020)
###  Frontiers in Ecology and Evolution	

### doi: 10.3389/fevo.2020.570034

### Run with R version 3.6.3
##########################################################################################  


### Set working directory

setwd("C:/PATH")

### Load required packages

library(lubridate)
library(plyr)
library(dplyr)
library(lme4)
library(climwin)


### Load nestling data

nestlings <- read.csv('alpine_songbird_nestling_data.csv')


### Split data into HOLA, DEJU, and SAVS datasets

HOLA_data <- subset(nestlings, species == "HOLA")
DEJU_data <- subset(nestlings, species == "DEJU")
SAVS_data <- subset(nestlings, species == "SAVS")


### Load weather data

weather_data <- read.csv("raw_weather_data.csv")

### Split weather data by species

HOLA_weather <- subset(weather_data, species == "HOLA")
DEJU_weather <- subset(weather_data, species == "DEJU")
SAVS_weather <- subset(weather_data, species == "SAVS")


###########################################################################################
### (1) Horned Lark
###########################################################################################


#### Convert measurement date into proper date format

HOLA_data$measured <- strptime(paste(HOLA_data$year, HOLA_data$measurement_date), format="%Y %j") 
HOLA_data$measured_date_formatted <- format(HOLA_data$measured, format= "%d/%m/%Y")


### Select only nestlings with size measurements

HOLA_data_use <- subset(HOLA_data, nestling_size_data == 1)


### Convert year and nest ID to factor

HOLA_data_use$year_f <- factor(HOLA_data_use$year)
HOLA_data_use$nest_id <- factor(HOLA_data_use$nest_id)


### HOLA sample size
length(table(HOLA_data_use$nest_id)) # Number of nests = 110
length(HOLA_data_use$species) # 361 nestlings


##################
### Incorporate weather data 
##################

### Convert weather date to correct format

HOLA_weather$date.new <- strptime(paste(HOLA_weather$year, HOLA_weather$j_date), format="%Y %j") 
HOLA_weather$date_use <- format(HOLA_weather$date.new, format= "%d/%m/%Y")

head(HOLA_weather)

### Remove data.new to avoid future conflicts

HOLA_weather <- HOLA_weather[,-11]



############################################################################################
############################# Horned Lark sliding window ###################################
############################################################################################

### Background:
### Sliding window set to maximum of 30 days prior to measurement date (7-days post-hatch).
### For an average nest (4 eggs, 1 egg laid a day, incubation on penultimate egg, and
### 12 days incubation), this means earliest time is 9 days prior to the first egg laid.
### While it is possible weather events prior to arrival at the breeding site may influence
### maternal condition and offspring development, addressing this is not possible with the 
### available data. The point of this analysis is to test the influence of local weather 
### events on the breeding grounds that influence breeding decisions and offspring 
### development. Additionally, egg formation (when nutrients and hormones are deposited
### in the egg) occurs about 4-5 days prior to the lay date, so we believe 9 days should
### capture this period adequately. 


### Standardized steps for identifying important time periods:
# For each weather variable and trait, the top period >= 3 days is selected by AIC.
# If several top periods differ by less than 2 AIC, then the period with the strongest Beta 
# coefficient is selected.
# If there are radically different time periods within the top models, each is selected
# and competed against each other in Step 2: Model fitting and selection.

### See linked paper for more details:
### doi: 10.3389/fevo.2020.570034


###############
### (A) Wing length
###############


### Test if year should be included to control for among year differences

wo_year_model <- lmer(wing ~ brood_size + measurement_age + DFE + HOLA_treatment + (1|nest_id), data=HOLA_data_use)
year_model <- lmer(wing ~ brood_size + measurement_age + DFE + HOLA_treatment + year_f + (1|nest_id), data=HOLA_data_use)

anova(wo_year_model, year_model, test="Chisq")

### Decision: No evidence that year improves the model. Do not include.


### Sliding windows

### Temperature variables

HOLA_wing_sw_temp <- slidingwin(baseline = lmer(wing ~ brood_size + measurement_age + DFE + HOLA_treatment + (1|nest_id), data=HOLA_data_use),  ### this is the base model to which the program adds the climate variables
                         xvar = list(TempMean=HOLA_weather$temp_avg_day,
                                     SubZero = HOLA_weather$temp_sub_zero_hrs,
                                     SubFive = HOLA_weather$temp_sub_5_hrs, 
                                     SubTen  = HOLA_weather$temp_sub_10_hrs),
                         type = "relative", 
                         range= c(30,0), 
                         stat = c("mean", "min", "max", "var"),
                         cmissing = "method1",
                         func = "lin",
                         cinterval = "day",
                         cdate = HOLA_weather$date_use, bdate = HOLA_data_use$measured_date_formatted) 

### Temperature results
### Note: assessed mean and variance for each tempreature variable, as well as the min and max
### for the average daily temperature within the top time periods.

HOLA_wing_sw_temp$combos

head(HOLA_wing_sw_temp[[1]]$Dataset, 20) ## Mean avg temp 23-14 (-)  (-1.03)
head(HOLA_wing_sw_temp[[2]]$Dataset, 20) ## Mean sub-zero hours (none better than the null)
head(HOLA_wing_sw_temp[[3]]$Dataset, 20) ## Mean sub 5 hours 14-1 and 7-1 (-) (-0.71)
head(HOLA_wing_sw_temp[[4]]$Dataset, 20) ## Mean sub 10 hours (23-14) (+)

head(HOLA_wing_sw_temp[[5]]$Dataset, 20) ## Min avg temp 27-16 (-)
head(HOLA_wing_sw_temp[[9]]$Dataset, 20) ## Max avg temp 6-2 (+)

head(HOLA_wing_sw_temp[[13]]$Dataset, 20) ## Variance daily temp (none better than null)
head(HOLA_wing_sw_temp[[14]]$Dataset, 20) ## Variance Sub-zero hours (none better than the null)
head(HOLA_wing_sw_temp[[15]]$Dataset, 20) ## Variance Sub 5 hours (none better than null)
head(HOLA_wing_sw_temp[[16]]$Dataset, 20) ## Variance Sub 10 hours (none better than null)



### Precipitation variables

HOLA_wing_sw_precip <- slidingwin(baseline = lmer(wing ~ brood_size + measurement_age + DFE + HOLA_treatment + (1|nest_id), data=HOLA_data_use),  ### this is the base model to which the program adds the climate variables
                                xvar = list(PrecipDays = HOLA_weather$precip_days,
                                            Storms = HOLA_weather$storms),
                                type = "relative", 
                                range= c(30,0), 
                                stat = c("mean"),
                                cmissing = "method1",
                                func = "lin",
                                cinterval = "day",
                                cdate = HOLA_weather$date_use, bdate = HOLA_data_use$measured_date_formatted) 

### Precipitation results
### Note: Assessed mean for precipitation days and storm events within each time
### period. Average and sum identify the same periods, but average has more top periods >3 days 
### so the results were easier to evaluate. 
### Sum and average should theoretically be highlighting the same stressor for a given time period.

HOLA_wing_sw_precip$combos

head(HOLA_wing_sw_precip[[1]]$Dataset, 20) ## Mean precip days (28-11) (+) 8-2 (-) (16.74)
head(HOLA_wing_sw_precip[[2]]$Dataset, 20) ## Mean storms (12-2) (-) (-45.09)



###########
### Randomizations
###########

### Mean of all temperature variables
HOLA_wing_rand_temp_mean <- randwin(repeats = 100, 
                xvar = list(TempMean=HOLA_weather$temp_avg_day, 
                            #SubZero = HOLA_weather$temp_sub_zero_hrs, # commented out because no different from null
                            SubFive = HOLA_weather$temp_sub_5_hrs, 
                            SubTen  = HOLA_weather$temp_sub_10_hrs),
                cdate = HOLA_weather$date_use, 
                bdate = HOLA_data_use$measured_date_formatted, 
                baseline = lmer(wing ~ brood_size + measurement_age + DFE + HOLA_treatment + (1|nest_id),
                                data=HOLA_data_use), 
                range = c(30, 0), 
                type = "relative", stat = c("mean"), 
                func = "lin", cinterval = "day")

### Min and max for average temp
HOLA_wing_rand_temp_min_max <- randwin(repeats = 100, 
                                    xvar = list(TempMean=HOLA_weather$temp_avg_day),
                                    cdate = HOLA_weather$date_use, 
                                    bdate = HOLA_data_use$measured_date_formatted, 
                                    baseline = lmer(wing ~ brood_size + measurement_age + DFE + HOLA_treatment + (1|nest_id),
                                                    data=HOLA_data_use), 
                                    range = c(30, 0), 
                                    type = "relative", stat = c("min", "max"), 
                                    func = "lin", cinterval = "day")

### Precip variables
HOLA_wing_rand_precip <- randwin(repeats = 100, 
                               xvar = list(PrecipDays = HOLA_weather$precip_days,
                                           Storms = HOLA_weather$storms),
                               cdate = HOLA_weather$date_use, 
                               bdate = HOLA_data_use$measured_date_formatted, 
                               baseline = lmer(wing ~ brood_size + measurement_age + DFE + HOLA_treatment + (1|nest_id),
                                               data=HOLA_data_use), 
                               range = c(30, 0), 
                               type = "relative", stat = c("mean"), 
                               func = "lin", cinterval = "day")


#########
### Compare randomizations to original sliding window
#########

# Mean avg temp
pvalue(dataset = HOLA_wing_sw_temp[[1]]$Dataset, datasetrand = HOLA_wing_rand_temp_mean[[1]], 
       metric="AIC", sample.size=361) # 0.01

# Mean sub 5 hrs
pvalue(dataset = HOLA_wing_sw_temp[[3]]$Dataset, datasetrand = HOLA_wing_rand_temp_mean[[2]], 
       metric="AIC", sample.size=361) # <0.001

# Mean sub 10 hrs
pvalue(dataset = HOLA_wing_sw_temp[[4]]$Dataset, datasetrand = HOLA_wing_rand_temp_mean[[3]], 
       metric="AIC", sample.size=361) # 0.02


# Min avg temp
pvalue(dataset = HOLA_wing_sw_temp[[5]]$Dataset, datasetrand = HOLA_wing_rand_temp_min_max[[1]], 
       metric="AIC", sample.size=361) # 0.03

# Max avg temp
pvalue(dataset = HOLA_wing_sw_temp[[9]]$Dataset, datasetrand = HOLA_wing_rand_temp_min_max[[2]], 
       metric="AIC", sample.size=361) # 0.02


# Mean precip days
pvalue(dataset = HOLA_wing_sw_precip[[1]]$Dataset, datasetrand = HOLA_wing_rand_precip[[1]], 
       metric="AIC", sample.size=361)  # 0.01

# Mean storm events
pvalue(dataset = HOLA_wing_sw_precip[[2]]$Dataset, datasetrand = HOLA_wing_rand_precip[[2]], 
       metric="AIC", sample.size=361) # 0.07


#################
### Conclusion:
### Retain mean avg temp, mean sub 5 hrs, mean sub 10 hrs, min avg temp, max avg temp,
### mean precip days, and mean storms.



###############################
### (B) Tarsus length 
###############################

### Test if year should be included to control for among year differences

tarsus_wo_year <- lmer(tarsus ~ brood_size + measurement_age + DFE + HOLA_treatment + (1|nest_id), data=HOLA_data_use)
tarsus_year <- lmer(tarsus ~ brood_size + measurement_age + DFE + HOLA_treatment + year_f + (1|nest_id), data=HOLA_data_use)

anova(tarsus_wo_year, tarsus_year, test="Chisq")

### Decision: Do not include.


### Sliding windows

# Temperature variables
HOLA_tarsus_sw_temp <- slidingwin(baseline = lmer(tarsus ~ brood_size + measurement_age + DFE + HOLA_treatment + (1|nest_id), data=HOLA_data_use),  ### this is the base model to which the program adds the climate variables
                       xvar = list(TempMean=HOLA_weather$temp_avg_day,
                                   SubZero = HOLA_weather$temp_sub_zero_hrs,
                                   SubFive = HOLA_weather$temp_sub_5_hrs, 
                                   SubTen  = HOLA_weather$temp_sub_10_hrs), 
                       type = "relative", 
                       range= c(30,0), 
                       stat = c("mean","min","max","var"), 
                       cmissing = "method1",
                       func = "lin",
                       cinterval = "day",
                       cdate = HOLA_weather$date_use, bdate = HOLA_data_use$measured_date_formatted) 

# Temperature results
HOLA_tarsus_sw_temp$combos

head(HOLA_tarsus_sw_temp[[1]]$Dataset, 20) ## Mean avg temp 8-0 (+)   (0.21)
head(HOLA_tarsus_sw_temp[[2]]$Dataset, 20) ## Mean sub zero hrs - top model is not different from null
head(HOLA_tarsus_sw_temp[[2]]$Dataset, 20) ## Mean sub 5 hrs 14-0 (-) 7-1 (-)
head(HOLA_tarsus_sw_temp[[4]]$Dataset, 20) ## mean sub 10 hrs - top model not sig different from null

head(HOLA_tarsus_sw_temp[[5]]$Dataset, 20) ## Min avg temp 8-2 (+)
head(HOLA_tarsus_sw_temp[[9]]$Dataset, 20) ## Max avg temp 5-2 (+)

head(HOLA_tarsus_sw_temp[[13]]$Dataset, 20) ## Variance avg temp (None better than null)
head(HOLA_tarsus_sw_temp[[14]]$Dataset, 20) ## Variance hrs < 0 (None better than null)
head(HOLA_tarsus_sw_temp[[15]]$Dataset, 20) ## Variance hrs < 5 (None better than null)
head(HOLA_tarsus_sw_temp[[16]]$Dataset, 20) ## Variance hrs < 10 (None better than null)


# Precipitation variables
HOLA_tarsus_sw_precip <- slidingwin(baseline = lmer(tarsus ~ brood_size + measurement_age + DFE + HOLA_treatment + (1|nest_id), data=HOLA_data_use),  ### this is the base model to which the program adds the climate variables
                                  xvar = list(PrecipDays = HOLA_weather$precip_days,
                                              Storms = HOLA_weather$storms), 
                                  type = "relative", 
                                  range= c(30,0), 
                                  stat = c("mean"), 
                                  cmissing = "method1",
                                  func = "lin",
                                  cinterval = "day",
                                  cdate = HOLA_weather$date_use, bdate = HOLA_data_use$measured_date_formatted) 

# Precipitation results
HOLA_tarsus_sw_precip$combos

head(HOLA_tarsus_sw_precip[[1]]$Dataset, 20) ## Mean precip days 9-0 (-)
head(HOLA_tarsus_sw_precip[[2]]$Dataset, 20) ## Mean storms 20-2 (-)



########
### Randomization HOLA tarsus
########

# Mean of all temperature variables
HOLA_tarsus_rand_temp_mean <- randwin(repeats = 100, 
                xvar = list(TempMean=HOLA_weather$temp_avg_day,
                            #SubZero = HOLA_weather$temp_sub_zero_hrs,
                            SubFive = HOLA_weather$temp_sub_5_hrs), 
                            #SubTen  = HOLA_weather$temp_sub_10_hrs), # Commented out because did not differ from null
                cdate = HOLA_weather$date_use, 
                bdate = HOLA_data_use$measured_date_formatted, 
                baseline = lmer(tarsus ~ brood_size + measurement_age + DFE + HOLA_treatment + (1|nest_id),
                                data=HOLA_data_use), 
                range = c(30, 0), 
                type = "relative", stat = c("mean"), 
                func = "lin", cinterval = "day")


# Min and max for avg temperature
HOLA_tarsus_rand_temp_min_max <- randwin(repeats = 100, 
                                        xvar = list(TempMean=HOLA_weather$temp_avg_day),
                                        cdate = HOLA_weather$date_use, 
                                        bdate = HOLA_data_use$measured_date_formatted, 
                                        baseline = lmer(tarsus ~ brood_size + measurement_age + DFE + HOLA_treatment + (1|nest_id),
                                                 data=HOLA_data_use), 
                                        range = c(30, 0), 
                                        type = "relative", stat = c("min", "max"), 
                                        func = "lin", cinterval = "day")


# Precipitation variables
HOLA_tarsus_rand_temp_precip <- randwin(repeats = 100, 
                                        xvar = list(PrecipDays = HOLA_weather$precip_days,
                                                     Storms = HOLA_weather$storms),
                                        cdate = HOLA_weather$date_use, 
                                        bdate = HOLA_data_use$measured_date_formatted, 
                                        baseline = lmer(tarsus ~ brood_size + measurement_age + DFE + HOLA_treatment + (1|nest_id),
                                                         data=HOLA_data_use), 
                                        range = c(30, 0), 
                                        type = "relative", stat = c("mean"), 
                                        func = "lin", cinterval = "day")


#########
### Compare randomizations to original sliding window
#########

# Mean avg temp
pvalue(dataset = HOLA_tarsus_sw_temp[[1]]$Dataset, datasetrand = HOLA_tarsus_rand_temp_mean[[1]], 
       metric="AIC", sample.size=361) # < 0.001 

# Mean sub 5 hrs
pvalue(dataset = HOLA_tarsus_sw_temp[[3]]$Dataset, datasetrand = HOLA_tarsus_rand_temp_mean[[2]], 
       metric="AIC", sample.size=361) # < 0.001

# Min avg temp
pvalue(dataset = HOLA_tarsus_sw_temp[[5]]$Dataset, datasetrand = HOLA_tarsus_rand_temp_min_max[[1]], 
       metric="AIC", sample.size=361) # < 0.001 

# Max avg temp
pvalue(dataset = HOLA_tarsus_sw_temp[[9]]$Dataset, datasetrand = HOLA_tarsus_rand_temp_min_max[[2]], 
       metric="AIC", sample.size=361) # < 0.001


# Mean precip days
pvalue(dataset = HOLA_tarsus_sw_precip[[1]]$Dataset, datasetrand = HOLA_tarsus_rand_temp_precip[[1]], 
       metric="AIC", sample.size=361) # 0.02

# Mean storm events
pvalue(dataset = HOLA_tarsus_sw_precip[[2]]$Dataset, datasetrand = HOLA_tarsus_rand_temp_precip[[2]], 
       metric="AIC", sample.size=361) # < 0.001



#########
### Conclusion: 
### Retain mean avg temp, mean sub 5 hrs, mean min temp, mean max temp, mean precip days, 
### and mean storms.



###############################
### (C) HOLA mass
###############################

### Test if year should be included to control for among year differences

mass_wo_year <- lmer(mass ~ brood_size + measurement_age + DFE + HOLA_treatment + (1|nest_id), data=HOLA_data_use)
mass_year <- lmer(mass ~ brood_size + measurement_age + DFE + HOLA_treatment + year_f + (1|nest_id), data=HOLA_data_use)

anova(mass_wo_year, mass_year, test="Chisq")

### Decision: No evidence that year improves the model. Do not include.

##########
### Sliding windows
##########

### Temperature variables
HOLA_mass_sw_temp <- slidingwin(baseline = lmer(mass ~ brood_size + measurement_age + DFE + HOLA_treatment + (1|nest_id), data=HOLA_data_use),  ### this is the base model to which the program adds the climate variables
                          xvar = list(TempMean=HOLA_weather$temp_avg_day,
                                      SubZero = HOLA_weather$temp_sub_zero_hrs,
                                      SubFive = HOLA_weather$temp_sub_5_hrs, 
                                      SubTen  = HOLA_weather$temp_sub_10_hrs),
                          type = "relative", 
                          range= c(30,0), 
                          stat = c("mean","min","max","var"),
                          cmissing = "method1",
                          func = "lin",
                          cinterval = "day",
                          cdate = HOLA_weather$date_use, bdate = HOLA_data_use$measured_date_formatted) 

# Temperature results
HOLA_mass_sw_temp$combos

head(HOLA_mass_sw_temp[[1]]$Dataset, 20) ## Mean avg temp 23-14 (-)
head(HOLA_mass_sw_temp[[2]]$Dataset, 20) ## Mean sub zero hrs - not sig different from null
head(HOLA_mass_sw_temp[[2]]$Dataset, 20) ## Mean sub 5 hrs 7-0 (-)
head(HOLA_mass_sw_temp[[4]]$Dataset, 20) ## Mean sub 10 hrs 21-15 (+)

head(HOLA_mass_sw_temp[[5]]$Dataset, 20) ## Min avg temp 27-16 (-)
head(HOLA_mass_sw_temp[[9]]$Dataset, 20) ## Max avg temp 18-15 (-)

head(HOLA_mass_sw_temp[[13]]$Dataset, 20) ## Variance avg temp (None better than null)
head(HOLA_mass_sw_temp[[14]]$Dataset, 20) ## Variance Hrs < 0 19-14 (-)
head(HOLA_mass_sw_temp[[15]]$Dataset, 20) ## Variance Hrs < 5 (None better than null)
head(HOLA_mass_sw_temp[[16]]$Dataset, 20) ## Variance Hrs < 10 (None better than null)


### Precipitation variables
HOLA_mass_sw_precip <- slidingwin(baseline = lmer(mass ~ brood_size + measurement_age + DFE + HOLA_treatment + (1|nest_id), data=HOLA_data_use),  ### this is the base model to which the program adds the climate variables
                                xvar = list(PrecipDays = HOLA_weather$precip_days,
                                            Storms = HOLA_weather$storms),
                                type = "relative", 
                                range= c(30,0), 
                                stat = c("mean"),
                                cmissing = "method1",
                                func = "lin",
                                cinterval = "day",
                                cdate = HOLA_weather$date_use, bdate = HOLA_data_use$measured_date_formatted) 

# Precipitation results
HOLA_mass_sw_precip$combos

head(HOLA_mass_sw_precip[[1]]$Dataset, 20) ## Mean precip days 28-10 (+)
head(HOLA_mass_sw_precip[[2]]$Dataset, 20) ## Mean storms 20-2 (-)



##############
### Randomization HOLA mass
##############

# Mean for all temperature variables
HOLA_mass_rand_temp <- randwin(repeats = 100, 
                xvar = list(TempMean=HOLA_weather$temp_avg_day,
                            SubFive = HOLA_weather$temp_sub_5_hrs, 
                            SubTen  = HOLA_weather$temp_sub_10_hrs),
                cdate = HOLA_weather$date_use, 
                bdate = HOLA_data_use$measured_date_formatted, 
                baseline = lmer(mass ~ brood_size + measurement_age + DFE + HOLA_treatment + (1|nest_id),
                                data=HOLA_data_use), 
                range = c(30, 0), 
                type = "relative", stat = c("mean"), 
                func = "lin", cinterval = "day")

# Variance in sub zero hrs
HOLA_mass_rand_temp_var <- randwin(repeats = 100, 
                                   xvar = list( SubZero = HOLA_weather$temp_sub_zero_hrs),
                                   cdate = HOLA_weather$date_use, 
                                   bdate = HOLA_data_use$measured_date_formatted, 
                                   baseline = lmer(mass ~ brood_size + measurement_age + DFE + HOLA_treatment + (1|nest_id),
                                                   data=HOLA_data_use), 
                                   range = c(30, 0), 
                                   type = "relative", stat = c("var"), 
                                   func = "lin", cinterval = "day")

# Min and max for avg temp
HOLA_mass_rand_temp_min_max <- randwin(repeats = 100, 
                               xvar = list(TempMean=HOLA_weather$temp_avg_day),
                               cdate = HOLA_weather$date_use, 
                               bdate = HOLA_data_use$measured_date_formatted, 
                               baseline = lmer(mass ~ brood_size + measurement_age + DFE + HOLA_treatment + (1|nest_id),
                                               data=HOLA_data_use), 
                               range = c(30, 0), 
                               type = "relative", stat = c("min", "max"), 
                               func = "lin", cinterval = "day")

# Precipitation variables
HOLA_mass_rand_precip <- randwin(repeats = 100, 
                                       xvar = list(PrecipDays = HOLA_weather$precip_days,
                                                   Storms = HOLA_weather$storms),
                                       cdate = HOLA_weather$date_use, 
                                       bdate = HOLA_data_use$measured_date_formatted, 
                                       baseline = lmer(mass ~ brood_size + measurement_age + DFE + HOLA_treatment + (1|nest_id),
                                                       data=HOLA_data_use), 
                                       range = c(30, 0), 
                                       type = "relative", stat = c("mean"), 
                                       func = "lin", cinterval = "day")


#########
### Compare randomizations to original sliding window
#########

# Mean avg temp
pvalue(dataset = HOLA_mass_sw[[1]]$Dataset, datasetrand = HOLA_mass_rand[[1]], 
       metric="AIC", sample.size=361) # 0.04

# Mean sub 5 hrs
pvalue(dataset = HOLA_mass_sw[[3]]$Dataset, datasetrand = HOLA_mass_rand[[2]], 
       metric="AIC", sample.size=361) # <0.001

# Mean sub 10 hrs
pvalue(dataset = HOLA_mass_sw[[4]]$Dataset, datasetrand = HOLA_mass_rand[[3]], 
       metric="AIC", sample.size=361)  # 0.01

# Variance sub 0 hrs
pvalue(dataset = HOLA_mass_sw_temp[[14]]$Dataset, datasetrand = HOLA_mass_rand_temp_var[[1]], 
       metric="AIC", sample.size=361) # 0.06

# Min avg temp
pvalue(dataset = HOLA_mass_sw_avg[[5]]$Dataset, datasetrand = HOLA_mass_rand_avg[[1]], 
       metric="AIC", sample.size=361) # 0.04

# Max avg temp
pvalue(dataset = HOLA_mass_sw_avg[[9]]$Dataset, datasetrand = HOLA_mass_rand_avg[[2]], 
       metric="AIC", sample.size=361) # 0.12

# Mean precip days
pvalue(dataset = HOLA_mass_sw[[1]]$Dataset, datasetrand = HOLA_mass_rand[[1]], 
       metric="AIC", sample.size=361) # 0.07

# Mean storm events
pvalue(dataset = HOLA_mass_sw[[2]]$Dataset, datasetrand = HOLA_mass_rand[[2]], 
       metric="AIC", sample.size=361) # 0.37

##########
### Conclusion: 
### Retain mean avg temp, men sub 5 hrs, mean sub 10 hrs, var sub 0 hrs,
### min avg temp, and mean precip


###########################################################################################
### HOLA: Extract weather variables from time periods identified by sliding window ########
###########################################################################################


###############
### Avg temp
###############

### HOLA wing: Best window for avg temp (23-14)

HOLA_data_use$avg_temp_start_w <- HOLA_data_use$measurement_date-23
HOLA_data_use$avg_temp_end_w <- HOLA_data_use$measurement_date-14

### HOLA wing: Best window for min temp (27-16)

HOLA_data_use$min_temp_start_w <- HOLA_data_use$measurement_date-27
HOLA_data_use$min_temp_end_w <- HOLA_data_use$measurement_date-16

### HOLA wing: Best window for max temp (18-15)

HOLA_data_use$max_temp_start_w <- HOLA_data_use$measurement_date-18
HOLA_data_use$max_temp_end_w <- HOLA_data_use$measurement_date-15


### HOLA tarsus: Best window for avg temp (8-0)

HOLA_data_use$avg_temp_start_t <- HOLA_data_use$measurement_date-8
HOLA_data_use$avg_temp_end_t <- HOLA_data_use$measurement_date-0

### HOLA tarsus: Best window for min temp (8-1)

HOLA_data_use$min_temp_start_t <- HOLA_data_use$measurement_date-8
HOLA_data_use$min_temp_end_t <- HOLA_data_use$measurement_date-1

### HOLA tarsus: Best window for max temp (5-2)

HOLA_data_use$max_temp_start_t <- HOLA_data_use$measurement_date-5
HOLA_data_use$max_temp_end_t <- HOLA_data_use$measurement_date-2


### HOLA mass: Best window for avg temp (23-14)

HOLA_data_use$avg_temp_start_m <- HOLA_data_use$measurement_date-23
HOLA_data_use$avg_temp_end_m <- HOLA_data_use$measurement_date-14

### HOLA mass: Best window for min temp (27-16)

HOLA_data_use$min_temp_start_m <- HOLA_data_use$measurement_date-27
HOLA_data_use$min_temp_end_m <- HOLA_data_use$measurement_date-16


###############
### Sub 0 hrs
###############

### HOLA mass: Best window for variance sub 0 hrs (19-14)

HOLA_data_use$var_sub0_start_m <- HOLA_data_use$measurement_date-19
HOLA_data_use$var_sub0_end_m <- HOLA_data_use$measurement_date-14


###############
### Sub 5 hrs
###############

### HOLA wing: Best window for < 5 hrs (7-1 & 14-1; short & long)

HOLA_data_use$sub_5_short_start_w <- HOLA_data_use$measurement_date-7
HOLA_data_use$sub_5_short_end_w <- HOLA_data_use$measurement_date-1

HOLA_data_use$sub_5_long_start_w <- HOLA_data_use$measurement_date-14
HOLA_data_use$sub_5_long_end_w <- HOLA_data_use$measurement_date-1


### HOLA tarsus: Best window for < 5 hrs (14-0 & 7-1; long & short)

HOLA_data_use$sub_5_long_start_t <- HOLA_data_use$measurement_date-14
HOLA_data_use$sub_5_long_end_t <- HOLA_data_use$measurement_date-0

HOLA_data_use$sub_5_short_start_t <- HOLA_data_use$measurement_date-7
HOLA_data_use$sub_5_short_end_t <- HOLA_data_use$measurement_date-1



### HOLA mass: Best window for < 5 hrs (7-0)

HOLA_data_use$sub_5_start_m <- HOLA_data_use$measurement_date-7
HOLA_data_use$sub_5_end_m <- HOLA_data_use$measurement_date-0


###############
### Sub 10 hrs
###############

### HOLA wing: Best window for sub 10 hrs (23-14)

HOLA_data_use$sub_10_start_w <- HOLA_data_use$measurement_date-23
HOLA_data_use$sub_10_end_w <- HOLA_data_use$measurement_date-14


### HOLA mass: Best window for sub 10 hrs (21-15)

HOLA_data_use$sub_10_start_m <- HOLA_data_use$measurement_date-21
HOLA_data_use$sub_10_end_m <- HOLA_data_use$measurement_date-15



###############
### Precip days
###############


### HOLA wing: Best window for precip days (28-11 & 8-2; early, late)

HOLA_data_use$precip_early_start_w <- HOLA_data_use$measurement_date-28
HOLA_data_use$precip_early_end_w <- HOLA_data_use$measurement_date-11

HOLA_data_use$precip_late_start_w <- HOLA_data_use$measurement_date-8
HOLA_data_use$precip_late_end_w <- HOLA_data_use$measurement_date-2


### HOLA mass: Best window for precip days (28-10)

HOLA_data_use$precip_start_m <- HOLA_data_use$measurement_date-28
HOLA_data_use$precip_end_m <- HOLA_data_use$measurement_date-10



###############
### Storms
###############

### HOLA wing: Best window for storms (12-2)

HOLA_data_use$storms_start_w <- HOLA_data_use$measurement_date-12
HOLA_data_use$storms_end_w <- HOLA_data_use$measurement_date-2

### HOLA tarsus: Best window for storms (20-2)

HOLA_data_use$storms_start_t <- HOLA_data_use$measurement_date-20
HOLA_data_use$storms_end_t <- HOLA_data_use$measurement_date-2


###############
### Select entire 30 day period for all weather variables for summary table
###############

HOLA_data_use$overall_start <- HOLA_data_use$measurement_date-30
HOLA_data_use$overall_end <- HOLA_data_use$measurement_date-0



###########################################
### Convert julian dates to proper posix format
###########################################

### Convert weather dates to POSIXct or else the following code will not work 
### (i.e., it won't work with same date format as sliding window)

HOLA_weather$date.int <- strptime(paste(HOLA_weather$year, HOLA_weather$j_date), format="%Y %j") 
HOLA_weather$date_new <- as.POSIXct(HOLA_weather$date.int)

head(HOLA_weather)

### Remove date.int or else loop below will not work

HOLA_weather <- HOLA_weather[,-12]


###############
### Avg temp
###############

### HOLA wing length
HOLA_data_use$avg_temp_start_w_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$avg_temp_start_w), format="%Y %j") 
HOLA_data_use$avg_temp_start_w_date <- as.POSIXct(HOLA_data_use$avg_temp_start_w_int)

HOLA_data_use$avg_temp_end_w_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$avg_temp_end_w), format="%Y %j") 
HOLA_data_use$avg_temp_end_w_date <- as.POSIXct(HOLA_data_use$avg_temp_end_w_int)

### HOLA tarsus length
HOLA_data_use$avg_temp_start_t_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$avg_temp_start_t), format="%Y %j") 
HOLA_data_use$avg_temp_start_t_date <- as.POSIXct(HOLA_data_use$avg_temp_start_t_int)

HOLA_data_use$avg_temp_end_t_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$avg_temp_end_t), format="%Y %j") 
HOLA_data_use$avg_temp_end_t_date <- as.POSIXct(HOLA_data_use$avg_temp_end_t_int)

### HOLA mass
HOLA_data_use$avg_temp_start_m_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$avg_temp_start_m), format="%Y %j") 
HOLA_data_use$avg_temp_start_m_date <- as.POSIXct(HOLA_data_use$avg_temp_start_m_int)

HOLA_data_use$avg_temp_end_m_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$avg_temp_end_m), format="%Y %j") 
HOLA_data_use$avg_temp_end_m_date <- as.POSIXct(HOLA_data_use$avg_temp_end_m_int)


###############
### Min temp
###############

### HOLA wing length
HOLA_data_use$min_temp_start_w_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$min_temp_start_w), format="%Y %j") 
HOLA_data_use$min_temp_start_w_date <- as.POSIXct(HOLA_data_use$min_temp_start_w_int)

HOLA_data_use$min_temp_end_w_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$min_temp_end_w), format="%Y %j") 
HOLA_data_use$min_temp_end_w_date <- as.POSIXct(HOLA_data_use$min_temp_end_w_int)

### HOLA tarsus length
HOLA_data_use$min_temp_start_t_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$min_temp_start_t), format="%Y %j") 
HOLA_data_use$min_temp_start_t_date <- as.POSIXct(HOLA_data_use$min_temp_start_t_int)

HOLA_data_use$min_temp_end_t_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$min_temp_end_t), format="%Y %j") 
HOLA_data_use$min_temp_end_t_date <- as.POSIXct(HOLA_data_use$min_temp_end_t_int)

### HOLA mass
HOLA_data_use$min_temp_start_m_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$min_temp_start_m), format="%Y %j") 
HOLA_data_use$min_temp_start_m_date <- as.POSIXct(HOLA_data_use$min_temp_start_m_int)

HOLA_data_use$min_temp_end_m_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$min_temp_end_m), format="%Y %j") 
HOLA_data_use$min_temp_end_m_date <- as.POSIXct(HOLA_data_use$min_temp_end_m_int)



###############
### Max temp
###############

### HOLA wing length
HOLA_data_use$max_temp_start_w_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$max_temp_start_w), format="%Y %j") 
HOLA_data_use$max_temp_start_w_date <- as.POSIXct(HOLA_data_use$max_temp_start_w_int)

HOLA_data_use$max_temp_end_w_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$max_temp_end_w), format="%Y %j") 
HOLA_data_use$max_temp_end_w_date <- as.POSIXct(HOLA_data_use$max_temp_end_w_int)

### HOLA tarsus length
HOLA_data_use$max_temp_start_t_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$max_temp_start_t), format="%Y %j") 
HOLA_data_use$max_temp_start_t_date <- as.POSIXct(HOLA_data_use$max_temp_start_t_int)

HOLA_data_use$max_temp_end_t_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$max_temp_end_t), format="%Y %j") 
HOLA_data_use$max_temp_end_t_date <- as.POSIXct(HOLA_data_use$max_temp_end_t_int)



###############
### Sub 0 hrs
###############

### HOLA mass
HOLA_data_use$var_sub0_start_m_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$var_sub0_start_m), format="%Y %j") 
HOLA_data_use$var_sub0_start_m_date <- as.POSIXct(HOLA_data_use$var_sub0_start_m_int)

HOLA_data_use$var_sub0_end_m_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$var_sub0_end_m), format="%Y %j") 
HOLA_data_use$var_sub0_end_m_date <- as.POSIXct(HOLA_data_use$var_sub0_end_m_int)



###############
### Sub 5 hrs
###############


### HOLA wing length
HOLA_data_use$sub_5_short_start_w_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$sub_5_short_start_w), format="%Y %j") 
HOLA_data_use$sub_5_short_start_w_date <- as.POSIXct(HOLA_data_use$sub_5_short_start_w_int)

HOLA_data_use$sub_5_short_end_w_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$sub_5_short_end_w), format="%Y %j") 
HOLA_data_use$sub_5_short_end_w_date <- as.POSIXct(HOLA_data_use$sub_5_short_end_w_int)

HOLA_data_use$sub_5_long_start_w_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$sub_5_long_start_w), format="%Y %j") 
HOLA_data_use$sub_5_long_start_w_date <- as.POSIXct(HOLA_data_use$sub_5_long_start_w_int)

HOLA_data_use$sub_5_long_end_w_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$sub_5_long_end_w), format="%Y %j") 
HOLA_data_use$sub_5_long_end_w_date <- as.POSIXct(HOLA_data_use$sub_5_long_end_w_int)


### HOLA tarsus
HOLA_data_use$sub_5_short_start_t_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$sub_5_short_start_t), format="%Y %j") 
HOLA_data_use$sub_5_short_start_t_date <- as.POSIXct(HOLA_data_use$sub_5_short_start_t_int)

HOLA_data_use$sub_5_short_end_t_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$sub_5_short_end_t), format="%Y %j") 
HOLA_data_use$sub_5_short_end_t_date <- as.POSIXct(HOLA_data_use$sub_5_short_end_t_int)

HOLA_data_use$sub_5_long_start_t_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$sub_5_long_start_t), format="%Y %j") 
HOLA_data_use$sub_5_long_start_t_date <- as.POSIXct(HOLA_data_use$sub_5_long_start_t_int)

HOLA_data_use$sub_5_long_end_t_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$sub_5_long_end_t), format="%Y %j") 
HOLA_data_use$sub_5_long_end_t_date <- as.POSIXct(HOLA_data_use$sub_5_long_end_t_int)

### HOLA mass
HOLA_data_use$sub_5_start_m_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$sub_5_start_m), format="%Y %j") 
HOLA_data_use$sub_5_start_m_date <- as.POSIXct(HOLA_data_use$sub_5_start_m_int)

HOLA_data_use$sub_5_end_m_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$sub_5_end_m), format="%Y %j") 
HOLA_data_use$sub_5_end_m_date <- as.POSIXct(HOLA_data_use$sub_5_end_m_int)



###############
### Sub 10 hrs
###############


### Wing length
HOLA_data_use$sub_10_start_w_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$sub_10_start_w), format="%Y %j") 
HOLA_data_use$sub_10_start_w_date <- as.POSIXct(HOLA_data_use$sub_10_start_w_int)

HOLA_data_use$sub_10_end_w_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$sub_10_end_w), format="%Y %j") 
HOLA_data_use$sub_10_end_w_date <- as.POSIXct(HOLA_data_use$sub_10_end_w_int)


### Mass
HOLA_data_use$sub_10_start_m_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$sub_10_start_m), format="%Y %j") 
HOLA_data_use$sub_10_start_m_date <- as.POSIXct(HOLA_data_use$sub_10_start_m_int)

HOLA_data_use$sub_10_end_m_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$sub_10_end_m), format="%Y %j") 
HOLA_data_use$sub_10_end_m_date <- as.POSIXct(HOLA_data_use$sub_10_end_m_int)


###############
### Precip
###############


### Wing length
HOLA_data_use$precip_early_start_w_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$precip_early_start_w), format="%Y %j") 
HOLA_data_use$precip_early_start_w_date <- as.POSIXct(HOLA_data_use$precip_early_start_w_int)

HOLA_data_use$precip_early_end_w_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$precip_early_end_w), format="%Y %j") 
HOLA_data_use$precip_early_end_w_date <- as.POSIXct(HOLA_data_use$precip_early_end_w_int)


HOLA_data_use$precip_late_start_w_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$precip_late_start_w), format="%Y %j") 
HOLA_data_use$precip_late_start_w_date <- as.POSIXct(HOLA_data_use$precip_late_start_w_int)

HOLA_data_use$precip_late_end_w_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$precip_late_end_w), format="%Y %j") 
HOLA_data_use$precip_late_end_w_date <- as.POSIXct(HOLA_data_use$precip_late_end_w_int)


### Mass
HOLA_data_use$precip_start_m_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$precip_start_m), format="%Y %j") 
HOLA_data_use$precip_start_m_date <- as.POSIXct(HOLA_data_use$precip_start_m_int)

HOLA_data_use$precip_end_m_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$precip_end_m), format="%Y %j") 
HOLA_data_use$precip_end_m_date <- as.POSIXct(HOLA_data_use$precip_end_m_int)



###############
### Storms
###############


### Wing length
HOLA_data_use$storms_start_w_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$storms_start_w), format="%Y %j") 
HOLA_data_use$storms_start_w_date <- as.POSIXct(HOLA_data_use$storms_start_w_int)

HOLA_data_use$storms_end_w_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$storms_end_w), format="%Y %j") 
HOLA_data_use$storms_end_w_date <- as.POSIXct(HOLA_data_use$storms_end_w_int)


### Tarsus length
HOLA_data_use$storms_start_t_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$storms_start_t), format="%Y %j") 
HOLA_data_use$storms_start_t_date <- as.POSIXct(HOLA_data_use$storms_start_t_int)

HOLA_data_use$storms_end_t_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$storms_end_t), format="%Y %j") 
HOLA_data_use$storms_end_t_date <- as.POSIXct(HOLA_data_use$storms_end_t_int)


### Overal time period

HOLA_data_use$overall_start_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$overall_start), format="%Y %j") 
HOLA_data_use$overall_start_date <- as.POSIXct(HOLA_data_use$overall_start_int)

HOLA_data_use$overall_end_int <- strptime(paste(HOLA_data_use$year, HOLA_data_use$overall_end), format="%Y %j") 
HOLA_data_use$overall_end_date <- as.POSIXct(HOLA_data_use$overall_end_int)




#############################################
### Create loop to extract weather variables for each nestling and save as a csv
#############################################

### Empty lists

avg_temp_w_list <- list()
avg_temp_t_list <- list()
avg_temp_m_list <- list()

min_temp_w_list <- list()
min_temp_t_list <- list()
min_temp_m_list <- list()

max_temp_w_list <- list()
max_temp_t_list <- list()

sub_0_m_list <- list()

sub_5_short_w_list <- list()
sub_5_long_w_list <- list()
sub_5_short_t_list <- list()
sub_5_long_t_list <- list()
sub_5_m_list <- list()

sub_10_w_list <- list()
sub_10_m_list <- list()

precip_early_w_list <- list()
precip_late_w_list <- list()
precip_m_list <- list()

storms_w_list <- list()
storms_t_list <- list()

overall_list <- list()


for (i in 1:361) {
        
        avg_temp_w_list[[i]] <- HOLA_weather %>% 
                filter(date_new >= HOLA_data_use$avg_temp_start_w_date[i] 
                       & date_new <= HOLA_data_use$avg_temp_end_w_date[i]) %>% 
                
                summarise(avg_temp_w = mean(temp_avg_day))
        
        avg_temp_t_list[[i]] <- HOLA_weather %>% 
                filter(date_new >= HOLA_data_use$avg_temp_start_t_date[i] 
                       & date_new <= HOLA_data_use$avg_temp_end_t_date[i]) %>% 
                
                summarise(avg_temp_t = mean(temp_avg_day))
        
        avg_temp_m_list[[i]] <- HOLA_weather %>% 
                filter(date_new >= HOLA_data_use$avg_temp_start_m_date[i] 
                       & date_new <= HOLA_data_use$avg_temp_end_m_date[i]) %>% 
          
                summarise(avg_temp_m = mean(temp_avg_day))
        
        min_temp_w_list[[i]] <- HOLA_weather %>% 
          filter(date_new >= HOLA_data_use$min_temp_start_w_date[i] 
                 & date_new <= HOLA_data_use$min_temp_end_w_date[i]) %>% 
          
          summarise(min_temp_w = min(temp_avg_day))
        
        min_temp_t_list[[i]] <- HOLA_weather %>% 
          filter(date_new >= HOLA_data_use$min_temp_start_t_date[i] 
                 & date_new <= HOLA_data_use$min_temp_end_t_date[i]) %>% 
          
          summarise(min_temp_t = min(temp_avg_day))
        
        min_temp_m_list[[i]] <- HOLA_weather %>% 
          filter(date_new >= HOLA_data_use$min_temp_start_m_date[i] 
                 & date_new <= HOLA_data_use$min_temp_end_m_date[i]) %>% 
          
          summarise(min_temp_m = min(temp_avg_day))
        
        max_temp_w_list[[i]] <- HOLA_weather %>% 
          filter(date_new >= HOLA_data_use$max_temp_start_w_date[i] 
                 & date_new <= HOLA_data_use$max_temp_end_w_date[i]) %>% 
          
          summarise(max_temp_w = max(temp_avg_day))
        
        max_temp_t_list[[i]] <- HOLA_weather %>% 
          filter(date_new >= HOLA_data_use$max_temp_start_t_date[i] 
                 & date_new <= HOLA_data_use$max_temp_end_t_date[i]) %>% 
          
          summarise(max_temp_t = max(temp_avg_day))
        
        
        sub_0_m_list[[i]] <- HOLA_weather %>% 
                filter(date_new >= HOLA_data_use$var_sub0_start_m_date[i] 
                     & date_new <= HOLA_data_use$var_sub0_end_m_date[i]) %>% 
                
                summarise(var_sub0_m = var(temp_sub_zero_hrs))
        
        
        sub_5_short_w_list[[i]] <- HOLA_weather %>% 
                filter(date_new >= HOLA_data_use$sub_5_short_start_w_date[i] 
                       & date_new <= HOLA_data_use$sub_5_short_end_w_date[i]) %>% 
                
                summarise(sub_5_short_w = mean(temp_sub_5_hrs))
        
        sub_5_long_w_list[[i]] <- HOLA_weather %>% 
                filter(date_new >= HOLA_data_use$sub_5_long_start_w_date[i] 
                       & date_new <= HOLA_data_use$sub_5_long_end_w_date[i]) %>% 
          
          summarise(sub_5_long_w = mean(temp_sub_5_hrs))
        
        sub_5_short_t_list[[i]] <- HOLA_weather %>% 
                filter(date_new >= HOLA_data_use$sub_5_short_start_t_date[i] 
                       & date_new <= HOLA_data_use$sub_5_short_end_t_date[i]) %>% 
                
                summarise(sub_5_short_t = mean(temp_sub_5_hrs))
        
        sub_5_long_t_list[[i]] <- HOLA_weather %>% 
                filter(date_new >= HOLA_data_use$sub_5_long_start_t_date[i] 
                      & date_new <= HOLA_data_use$sub_5_long_end_t_date[i]) %>% 
          
          summarise(sub_5_long_t = mean(temp_sub_5_hrs))
        
        sub_5_m_list[[i]] <- HOLA_weather %>% 
                filter(date_new >= HOLA_data_use$sub_5_start_m_date[i] 
                       & date_new <= HOLA_data_use$sub_5_end_m_date[i]) %>% 
                
                summarise(sub_5_m = mean(temp_sub_5_hrs))
        
        sub_10_w_list[[i]] <- HOLA_weather %>% 
                filter(date_new >= HOLA_data_use$sub_10_start_w_date[i] 
                       & date_new <= HOLA_data_use$sub_10_end_w_date[i]) %>% 
                
                summarise(sub_10_w = mean(temp_sub_10_hrs))
      
        sub_10_m_list[[i]] <- HOLA_weather %>% 
                filter(date_new >= HOLA_data_use$sub_10_start_m_date[i] 
                       & date_new <= HOLA_data_use$sub_10_end_m_date[i]) %>% 
                
                summarise(sub_10_m = mean(temp_sub_10_hrs))
        
        
        precip_early_w_list[[i]] <- HOLA_weather %>% 
                filter(date_new >= HOLA_data_use$precip_early_start_w_date[i] 
                       & date_new <= HOLA_data_use$precip_early_end_w_date[i]) %>% 
                
                summarise(total_early_precip_w = sum(precip_days))
        
        precip_late_w_list[[i]] <- HOLA_weather %>% 
                filter(date_new >= HOLA_data_use$precip_late_start_w_date[i] 
                       & date_new <= HOLA_data_use$precip_late_end_w_date[i]) %>% 
                
                summarise(total_late_precip_w = sum(precip_days))
        
        precip_m_list[[i]] <- HOLA_weather %>% 
                filter(date_new >= HOLA_data_use$precip_start_m_date[i] 
                       & date_new <= HOLA_data_use$precip_end_m_date[i]) %>% 
          
                summarise(total_precip_m = sum(precip_days))
        
        storms_w_list[[i]] <- HOLA_weather %>% 
                filter(date_new >= HOLA_data_use$storms_start_w_date[i] 
                       & date_new <= HOLA_data_use$storms_end_w_date[i]) %>% 
                
                summarise(total_storms_w = sum(storms))
        
        storms_t_list[[i]] <- HOLA_weather %>% 
                filter(date_new >= HOLA_data_use$storms_start_t_date[i] 
                       & date_new <= HOLA_data_use$storms_end_t_date[i]) %>% 
                
                summarise(total_storms_t = sum(storms))
        
        overall_list[[i]] <- HOLA_weather %>% 
                filter(date_new >= HOLA_data_use$overall_start_date[i] 
                       & date_new <= HOLA_data_use$overall_end_date[i]) %>% 
          
                summarise(avg_temp = mean(temp_avg_day),
                          min_temp = min(temp_avg_day),
                          max_temp = max(temp_avg_day),
                          avg_5 = mean(temp_sub_5_hrs),
                          avg_10 = mean(temp_sub_10_hrs),
                          total_precip = sum(precip_days),
                          total_storms = sum(storms))
        
        
        ### Convert lists to rows per nestling ID
        
        avg_temp_w_data <- bind_rows(avg_temp_w_list)
        avg_temp_t_data <- bind_rows(avg_temp_t_list)
        avg_temp_m_data <- bind_rows(avg_temp_m_list)
        
        min_temp_w_data <- bind_rows(min_temp_w_list)
        min_temp_t_data <- bind_rows(min_temp_t_list)
        min_temp_m_data <- bind_rows(min_temp_m_list)
        
        max_temp_w_data <- bind_rows(max_temp_w_list)
        max_temp_t_data <- bind_rows(max_temp_t_list)
        
        sub_0_m_data <- bind_rows(sub_0_m_list)
        
        sub_5_short_w_data <- bind_rows(sub_5_short_w_list)
        sub_5_long_w_data <- bind_rows(sub_5_long_w_list)
        sub_5_short_t_data <- bind_rows(sub_5_short_t_list)
        sub_5_long_t_data <- bind_rows(sub_5_long_t_list)
        sub_5_m_data <- bind_rows(sub_5_m_list)
        
        sub_10_w_data <- bind_rows(sub_10_w_list)
        sub_10_m_data <- bind_rows(sub_10_m_list)
        
        precip_early_w_data <- bind_rows(precip_early_w_list)
        precip_late_w_data <- bind_rows(precip_late_w_list)
        precip_m_data <- bind_rows(precip_m_list)
        
        storms_w_data <- bind_rows(storms_w_list)
        storms_t_data <- bind_rows(storms_t_list)
        
        overall_data <- bind_rows(overall_list)
        
        
        
        ### Convert each variable to a column in a larger dataset
        HOLA_windows <- bind_cols(avg_temp_w_data, avg_temp_t_data, avg_temp_m_data,
                                  min_temp_w_data, min_temp_t_data, min_temp_m_data,
                                  max_temp_w_data, max_temp_t_data, sub_0_m_data,
                                  sub_5_short_w_data, sub_5_long_w_data, 
                                  sub_5_short_t_data, sub_5_long_t_data, sub_5_m_data, 
                                  sub_10_w_data, sub_10_m_data, precip_early_w_data, 
                                  precip_late_w_data, precip_m_data, storms_w_data, 
                                  storms_t_data, overall_data)
        
}

head(HOLA_windows)


### Combine with nestling measurement data

head(HOLA_data_use)

HOLA_data_new <- HOLA_data_use[, c(1:12,16)] # Extract only columns of interest from original data

HOLA_final <- bind_cols(HOLA_data_new, HOLA_windows)


### Save data

write.csv(HOLA_final, "HOLA_weather_summarized_final.csv")



############################################################################################
### (2) Dark-eyed Junco sliding window #####################################################
############################################################################################


### Convert measurement date to proper format

DEJU_data$measured <- strptime(paste(DEJU_data$year, DEJU_data_use$measurement_date), format="%Y %j") 
DEJU_data$measured_date_formatted <- format(DEJU_data$measured, format= "%d/%m/%Y")


### Select only nestlings with size measurements

DEJU_data_use <- subset(DEJU_data, nestling_size_data == 1)


### Convert nest ID and year into a factor

DEJU_data_use$nest_id <- factor(DEJU_data_use$nest_id)

DEJU_data_use$year_f <- factor(DEJU_data_use$year)


### Calculate number of nestlings and nests

length(DEJU_data_use$species) # 120 nestlings

length(table(DEJU_data_use$nest_id)) # 35 nests


############
### Integrate weather data 
############

### Reformat

DEJU_weather$date.new <- strptime(paste(DEJU_weather$year, DEJU_weather$j_date), format="%Y %j") 
DEJU_weather$date_use <- format(DEJU_weather$date.new, format= "%d/%m/%Y")

### Remove old date column

head(DEJU_weather)
DEJU_weather <- DEJU_weather[,-11]



#######################################
#### DEJU sliding window #######
#######################################


###############################
### (A) wing length
###############################

### Note: Average clutch size = 4, average incubation = 13, so first egg laid 22 days prior to
### measurement

### Test if year should be included to control for among year differences

DEJU_wo_year <- lmer(wing ~ brood_size + measurement_age + DFE + (1|nest_id), data=DEJU_data_use)

DEJU_year <- lmer(wing ~ brood_size + measurement_age + DFE + year_f + (1|nest_id), data=DEJU_data_use)

anova(DEJU_wo_year, DEJU_year, test = "F") # P = 0.23

### Result: no evidence that year improves model. Do not include.

############
### Sliding windows
############

# Temperature variables
DEJU_wing_sw_temp <- slidingwin(baseline = lmer(wing ~ brood_size + measurement_age + DFE + (1|nest_id), data=DEJU_data_use),  ### this is the base model to which the program adds the climate variables
                       xvar = list(TempMean = DEJU_weather$temp_avg_day,
                                   SubZero = DEJU_weather$temp_sub_zero_hrs,
                                   SubFive = DEJU_weather$temp_sub_5_hrs, 
                                   SubTen  = DEJU_weather$temp_sub_10_hrs), 
                       type = "relative", 
                       range= c(30,0), 
                       stat = c("mean", "min", "max", "var"),
                       cmissing = "method1",
                       func = "lin",
                       cinterval = "day",
                       cdate = DEJU_weather$date_use, bdate = DEJU_data_use$measured_date_formatted) 

# Temperature results
DEJU_wing_sw_temp$combos

head(DEJU_wing_sw_temp[[1]]$Dataset, 20)  # Mean avg temp - (7-0) (+)
head(DEJU_wing_sw_temp[[2]]$Dataset, 20)  # Mean sub zero (not sig different from null)
head(DEJU_wing_sw_temp[[3]]$Dataset, 20)  # Mean sub 5 hrs 5-1 (-)
head(DEJU_wing_sw_temp[[4]]$Dataset, 20)  # Mean sub 10 hrs 7-0 (-)

head(DEJU_wing_sw_temp[[5]]$Dataset, 20) # Min avg temp (none better than the null)
head(DEJU_wing_sw_temp[[9]]$Dataset, 20) # Max avg temp 6-0 (0.93)

head(DEJU_wing_sw_temp[[13]]$Dataset, 20) # Avg temp var (19-16) (-0.50)
head(DEJU_wing_sw_temp[[14]]$Dataset, 20) # sub zero var - no better than null
head(DEJU_wing_sw_temp[[15]]$Dataset, 20) # sub 5 var - no better than null
head(DEJU_wing_sw_temp[[16]]$Dataset, 20) # sub 10 var (6-0) (0.16) 
# Note - variance effect sizes very small


# Precipitation variables
DEJU_wing_sw_precip <- slidingwin(baseline = lmer(wing ~ brood_size + measurement_age + DFE + (1|nest_id), data=DEJU_data_use),  ### this is the base model to which the program adds the climate variables
                                xvar = list(PrecipDays = DEJU_weather$precip_days,
                                            Storms = DEJU_weather$storms), 
                                type = "relative", 
                                range= c(30,0), 
                                stat = c("mean"),
                                cmissing = "method1",
                                func = "lin",
                                cinterval = "day",
                                cdate = DEJU_weather$date_use, bdate = DEJU_data_use$measured_date_formatted) 

# Precipitation results
DEJU_wing_sw_precip$combos

head(DEJU_wing_sw[[1]]$Dataset, 20)  # Precip days 8-0 (-)
head(DEJU_wing_sw[[2]]$Dataset, 20)  # Storms 6-0 (-) 28-21 (-)



##################
### Randomization for DEJU wing
##################

# Mean temperature variables
DEJU_wing_rand_temp_mean <- randwin(repeats = 100, 
                xvar = list(TempMean = DEJU_weather$temp_avg_day,
                            #SubZero = DEJU_weather$temp_sub_zero_hrs, # No different from null above
                            SubFive = DEJU_weather$temp_sub_5_hrs, 
                            SubTen  = DEJU_weather$temp_sub_10_hrs),
                cdate = DEJU_weather$date_use, 
                bdate = DEJU_data_use$measured_date_formatted, 
                baseline = lmer(wing ~ brood_size + measurement_age + DFE + (1|nest_id),
                                data=DEJU_data_use), 
                range = c(30, 0), 
                type = "relative", stat = c("mean"), 
                func = "lin", cinterval = "day")


# Min and max avg temperature
DEJU_wing_rand_temp_min_max <- randwin(repeats = 100, 
                                    xvar = list(TempMean = DEJU_weather$temp_avg_day),
                                    cdate = DEJU_weather$date_use, 
                                    bdate = DEJU_data_use$measured_date_formatted, 
                                    baseline = lmer(wing ~ brood_size + measurement_age + DFE + (1|nest_id),
                                                    data=DEJU_data_use), 
                                    range = c(30, 0), 
                                    type = "relative", stat = c("min", "max"), 
                                    func = "lin", cinterval = "day")


# Temperature variance
DEJU_wing_rand_temp_var <- randwin(repeats = 100, 
                                       xvar = list(TempMean = DEJU_weather$temp_avg_day),
                                       cdate = DEJU_weather$date_use, 
                                       bdate = DEJU_data_use$measured_date_formatted, 
                                       baseline = lmer(wing ~ brood_size + measurement_age + DFE + (1|nest_id),
                                                       data=DEJU_data_use), 
                                       range = c(30, 0), 
                                       type = "relative", stat = c("var"), 
                                       func = "lin", cinterval = "day")


# Precipitation variables
DEJU_wing_rand_precip <- randwin(repeats = 100, 
                                   xvar = list(PrecipDays = DEJU_weather$precip_days,
                                               Storms = DEJU_weather$storms),
                                   cdate = DEJU_weather$date_use, 
                                   bdate = DEJU_data_use$measured_date_formatted, 
                                   baseline = lmer(wing ~ brood_size + measurement_age + DFE + (1|nest_id),
                                                   data=DEJU_data_use), 
                                   range = c(30, 0), 
                                   type = "relative", stat = c("mean"), 
                                   func = "lin", cinterval = "day")


############
### Compare original and randomization
############

length(DEJU_data_use$species) # sample size to use = 120

# Mean avg temp
pvalue(dataset = DEJU_wing_sw_temp[[1]]$Dataset, datasetrand = DEJU_wing_rand_temp_mean[[1]], 
       metric="AIC", sample.size=120) # 0.03

# Mean sub 5 hrs
pvalue(dataset = DEJU_wing_sw_temp[[3]]$Dataset, datasetrand = DEJU_wing_rand_temp_mean[[2]], 
       metric="AIC", sample.size=120) # 0.05

# Mean sub 10 hrs
pvalue(dataset = DEJU_wing_sw_temp[[4]]$Dataset, datasetrand = DEJU_wing_rand_temp_mean[[3]], 
       metric="AIC", sample.size=120)  # <0.001

# Min avg temp
pvalue(dataset = DEJU_wing_sw_temp[[5]]$Dataset, datasetrand = DEJU_wing_rand_temp_min_max[[1]], 
       metric="AIC", sample.size=120) # 0.19

# Max avg temp
pvalue(dataset = DEJU_wing_sw_temp[[9]]$Dataset, datasetrand = DEJU_wing_rand_temp_min_max[[2]], 
       metric="AIC", sample.size=120) # < 0.001

# Var avg temp
pvalue(dataset = DEJU_wing_sw_temp[[13]]$Dataset, datasetrand = DEJU_wing_rand_temp_var[[1]], 
       metric="AIC", sample.size=120) # 0.31

# Mean precip days
pvalue(dataset = DEJU_wing_sw_precip[[1]]$Dataset, datasetrand = DEJU_wing_rand_precip[[1]], 
       metric="AIC", sample.size=120) # 0.1

# Mean storm events
pvalue(dataset = DEJU_wing_sw_precip[[2]]$Dataset, datasetrand = DEJU_wing_rand_precip[[1]], 
       metric="AIC", sample.size=120)  # 0.1

### Conclusion:
# Retain mean avg temp, mean sub 5 hrs, mean sub 10 hrs, and max avg temp.



###############################
### (B) Tarsus length 
###############################

### Test if year should be included to control for among year differences

DEJU_tarsus_wo_year <- lmer(tarsus ~ brood_size + measurement_age + DFE + (1|nest_id), data=DEJU_data_use)

DEJU_tarsus_year <- lmer(tarsus ~ brood_size + measurement_age + DFE + year_f + (1|nest_id), data=DEJU_data_use)

anova(DEJU_tarsus_wo_year, DEJU_tarsus_year, test = "F") # P = 0.07, AIC -1.24

### Result: weak evidence that year improves model. Do not include.

#######################
### Sliding windows
#######################

# Temperature variables
DEJU_tarsus_sw_temp <- slidingwin(baseline = lmer(tarsus ~ brood_size + measurement_age + DFE + (1|nest_id), data=DEJU_data_use),  ### this is the base model to which the program adds the climate variables
                          xvar = list(TempMean=DEJU_weather$temp_avg_day, 
                                      SubZero = DEJU_weather$temp_sub_zero_hrs,
                                      SubFive = DEJU_weather$temp_sub_5_hrs, 
                                      SubTen  = DEJU_weather$temp_sub_10_hrs),
                          type = "relative", 
                          range= c(30,0), 
                          stat = c("mean", "min", "max", "var"), 
                          cmissing = "method1",
                          func = "lin",
                          cinterval = "day",
                          cdate = DEJU_weather$date_use, bdate = DEJU_data_use$measured_date_formatted) 

# Temperature results
DEJU_tarsus_sw_temp$combos

head(DEJU_tarsus_sw_temp[[1]]$Dataset, 50) # Mean avg temp 7-0 (+) 30-19 (-)
head(DEJU_tarsus_sw_temp[[2]]$Dataset, 50) # Mean sub zero hrs - no different from null
head(DEJU_tarsus_sw_temp[[3]]$Dataset, 50) # Mean sub 5 hrs  6-0 (-)
head(DEJU_tarsus_sw_temp[[4]]$Dataset, 50) # Mean sub 10 hrs 7-0 (-) 30-17 (+)

head(DEJU_tarsus_sw_temp[[5]]$Dataset, 20) # Min avg temp  28-18 (-0.28)
head(DEJU_tarsus_sw_temp[[9]]$Dataset, 20) # Max avg temp  6-0 (0.28)

head(DEJU_tarsus_sw_temp[[13]]$Dataset, 20) # Variance avg temp (no different from null)
head(DEJU_tarsus_sw_temp[[14]]$Dataset, 20) # Variance sub zero (no better than null)
head(DEJU_tarsus_sw_temp[[15]]$Dataset, 20) # Variance sub 5  (no better than null)
head(DEJU_tarsus_sw_temp[[16]]$Dataset, 20) # Variance sub 10 (30-27) (-0.07) 


# Precipitation variables
DEJU_tarsus_sw_precip <- slidingwin(baseline = lmer(tarsus ~ brood_size + measurement_age + DFE + (1|nest_id), data=DEJU_data_use),  ### this is the base model to which the program adds the climate variables
                                  xvar = list(PrecipDays = DEJU_weather$precip_days,
                                              Storms = DEJU_weather$storms),
                                  type = "relative", 
                                  range= c(30,0), 
                                  stat = c("mean"), 
                                  cmissing = "method1",
                                  func = "lin",
                                  cinterval = "day",
                                  cdate = DEJU_weather$date_use, bdate = DEJU_data_use$measured_date_formatted) 

# Precipitation results
DEJU_tarsus_sw_precip$combos

head(DEJU_tarsus_sw_precip[[1]]$Dataset, 20) # Precip days 6-0 (-) 27-18 (+)
head(DEJU_tarsus_sw_precip[[2]]$Dataset, 20) # Storms 5-0 (-) 28-18 (-)


#################
### Randomization for DEJU tarsus
#################

# Mean temperature variables
DEJU_tarsus_rand_temp_mean <- randwin(repeats = 100, 
                        xvar = list(TempMean=DEJU_weather$temp_avg_day,
                                    SubFive = DEJU_weather$temp_sub_5_hrs, 
                                    SubTen  = DEJU_weather$temp_sub_10_hrs),
                cdate = DEJU_weather$date_use, 
                bdate = DEJU_data_use$measured_date_formatted, 
                baseline = lmer(tarsus ~ brood_size + measurement_age + DFE + (1|nest_id),
                                data=DEJU_data_use), 
                range = c(30, 0), 
                type = "relative", stat = c("mean"), 
                func = "lin", cinterval = "day")


# Min and max avg temperature
DEJU_tarsus_rand_temp_min_max <- randwin(repeats = 100, 
                                      xvar = list(TempMean=DEJU_weather$temp_avg_day),
                                      cdate = DEJU_weather$date_use, 
                                      bdate = DEJU_data_use$measured_date_formatted, 
                                      baseline = lmer(tarsus ~ brood_size + measurement_age + DFE + (1|nest_id),
                                                      data=DEJU_data_use), 
                                      range = c(30, 0), 
                                      type = "relative", stat = c("min", "max"), 
                                      func = "lin", cinterval = "day")

# Precipitation variables
DEJU_tarsus_rand_precip <- randwin(repeats = 100, 
                                         xvar = list(PrecipDays = DEJU_weather$precip_days,
                                                     Storms = DEJU_weather$storms),
                                         cdate = DEJU_weather$date_use, 
                                         bdate = DEJU_data_use$measured_date_formatted, 
                                         baseline = lmer(tarsus ~ brood_size + measurement_age + DFE + (1|nest_id),
                                                         data=DEJU_data_use), 
                                         range = c(30, 0), 
                                         type = "relative", stat = c("mean"), 
                                         func = "lin", cinterval = "day")


################
### Compare randomization to sliding window results
################

# Mean avg temp
pvalue(dataset = DEJU_tarsus_sw_temp[[1]]$Dataset, datasetrand = DEJU_tarsus_rand_temp_mean[[1]], 
       metric="AIC", sample.size=120) # 0.05

# Mean sub 5 hrs
pvalue(dataset = DEJU_tarsus_sw_temp[[3]]$Dataset, datasetrand = DEJU_tarsus_rand_temp_mean[[2]], 
       metric="AIC", sample.size=120) # 0.01

# Mean sub 10 hrs
pvalue(dataset = DEJU_tarsus_sw_temp[[4]]$Dataset, datasetrand = DEJU_tarsus_rand_temp_mean[[3]], 
       metric="AIC", sample.size=120)  # 0.01


# Min avg temp
pvalue(dataset = DEJU_tarsus_sw_temp[[5]]$Dataset, datasetrand = DEJU_tarsus_rand_temp_min_max[[1]], 
       metric="AIC", sample.size=120) # 0.11

# Max avg temp
pvalue(dataset = DEJU_tarsus_sw_temp[[9]]$Dataset, datasetrand = DEJU_tarsus_rand_temp_min_max[[2]], 
       metric="AIC", sample.size=120) # 0.01

# Mean precip days
pvalue(dataset = DEJU_tarsus_sw_precip[[1]]$Dataset, datasetrand = DEJU_tarsus_rand_precip[[1]], 
       metric="AIC", sample.size=120)  # 0.3

# Mean storm days
pvalue(dataset = DEJU_tarsus_sw_precip[[2]]$Dataset, datasetrand = DEJU_tarsus_rand_precip[[2]], 
       metric="AIC", sample.size=120)  # 0.11

##################
### Conclusion:
# Retain all except min avg temp, mean precip days, and mean storms which are not supported.



###############################
### (C) Mass
###############################

### Test if year should be included to control for among year differences

DEJU_mass_wo_year <- lmer(mass ~ brood_size + measurement_age + DFE + (1|nest_id), data=DEJU_data_use)

DEJU_mass_year <- lmer(mass ~ brood_size + measurement_age + DFE + year_f + (1|nest_id), data=DEJU_data_use)

anova(DEJU_mass_wo_year, DEJU_mass_year, test = "F") # P = 0.42

### Result: no evidence that year improves model. Do not include.

#########################
### Sliding windows
#########################

# Temperature variables
DEJU_mass_sw_temp <- slidingwin(baseline = lmer(mass ~ brood_size + measurement_age + DFE + (1|nest_id), data=DEJU_data_use),  ### this is the base model to which the program adds the climate variables
                          xvar = list(TempMean=DEJU_weather$temp_avg_day, 
                                      SubZero = DEJU_weather$temp_sub_zero_hrs,
                                      SubFive = DEJU_weather$temp_sub_5_hrs, 
                                      SubTen  = DEJU_weather$temp_sub_10_hrs), 
                          type = "relative", 
                          range= c(30,0), 
                          stat = c("mean","min","max","var"),
                          cmissing = "method1",
                          func = "lin",
                          cinterval = "day",
                          cdate = DEJU_weather$date_use, bdate = DEJU_data_use$measured_date_formatted) 
# Temperature results
DEJU_mass_sw_temp$combos

head(DEJU_mass_sw_temp[[1]]$Dataset, 20) # Mean avg temp (none different than the null)
head(DEJU_mass_sw_temp[[2]]$Dataset, 20) # Mean sub zero hrs (none different than the null)
head(DEJU_mass_sw_temp[[3]]$Dataset, 20) # Mean sub 5 hrs (none different than the null)
head(DEJU_mass_sw_temp[[4]]$Dataset, 20) # Mean sub 10 hrs (none different than the null)

head(DEJU_mass_sw_temp[[5]]$Dataset, 20) # Min avg temp (no better than null)
head(DEJU_mass_sw_temp[[9]]$Dataset, 20) # Max avg temp  6-0 (0.25)

head(DEJU_mass_sw_temp[[13]]$Dataset, 20) # Variance avg temp (no different from null)
head(DEJU_mass_sw_temp[[14]]$Dataset, 20) # Variance sub zero (no better than null)
head(DEJU_mass_sw_temp[[15]]$Dataset, 20) # Variance sub 5 (no better than null)
head(DEJU_mass_sw_temp[[16]]$Dataset, 20) # Variance sub 10 (no better than null) 


# Precipitation variables
DEJU_mass_sw_precip <- slidingwin(baseline = lmer(mass ~ brood_size + measurement_age + DFE + (1|nest_id), data=DEJU_data_use),  ### this is the base model to which the program adds the climate variables
                                xvar = list(PrecipDays = DEJU_weather$precip_days,
                                            Storms = DEJU_weather$storms), 
                                type = "relative", 
                                range= c(30,0), 
                                stat = c("mean"),
                                cmissing = "method1",
                                func = "lin",
                                cinterval = "day",
                                cdate = DEJU_weather$date_use, bdate = DEJU_data_use$measured_date_formatted) 
# Precipitation results
DEJU_mass_sw_precip$combos

head(DEJU_mass_sw_precip[[1]]$Dataset, 50) # Precip days - 27-7 (+) and 6-0 (-)
head(DEJU_mass_sw_precip[[2]]$Dataset, 50) # Storms 28-18 (-) 22-6 (+) 5-0 (-)


#################
### Randomization for DEJU mass
#################

# Max avg temperature
DEJU_mass_rand_temp_max <- randwin(repeats = 100, 
                xvar = list(TempMean=DEJU_weather$temp_avg_day),
                cdate = DEJU_weather$date_use, 
                bdate = DEJU_data_use$measured_date_formatted, 
                baseline = lmer(mass ~ brood_size + measurement_age + DFE + (1|nest_id),
                                data=DEJU_data_use), 
                range = c(30, 0), 
                type = "relative", stat = c("max"), 
                func = "lin", cinterval = "day")

# Precipitation variables
DEJU_mass_rand_precip <- randwin(repeats = 100, 
                                   xvar = list(PrecipDays = DEJU_weather$precip_days,
                                               Storms = DEJU_weather$storms),
                                   cdate = DEJU_weather$date_use, 
                                   bdate = DEJU_data_use$measured_date_formatted, 
                                   baseline = lmer(mass ~ brood_size + measurement_age + DFE + (1|nest_id),
                                                   data=DEJU_data_use), 
                                   range = c(30, 0), 
                                   type = "relative", stat = c("mean"), 
                                   func = "lin", cinterval = "day")

#################
### Compare between randomizations and sliding window
#################

# Max avg temp
pvalue(dataset = DEJU_mass_sw_temp[[9]]$Dataset, datasetrand = DEJU_mass_rand_temp_max[[1]], 
       metric="AIC", sample.size=120) # <0.001

# Mean precip days
pvalue(dataset = DEJU_mass_sw_precip[[1]]$Dataset, datasetrand = DEJU_mass_rand_precip[[1]], 
       metric="AIC", sample.size=120) # 0.29

# Mean storm events
pvalue(dataset = DEJU_mass_sw_precip[[2]]$Dataset, datasetrand = DEJU_mass_rand_precip[[2]], 
       metric="AIC", sample.size=120) # 0.06

#############################
### Conclusion:
### Retain only max avg temp and storm events.


###########################################################################################
### DEJU: Extract weather variables from time periods identified by sliding window ########
###########################################################################################


###############
### Avg temp
###############

### DEJU wing: Best window for avg temp (7-0)

DEJU_data_use$avg_temp_start_w <- DEJU_data_use$measurement_date-7
DEJU_data_use$avg_temp_end_w <- DEJU_data_use$measurement_date-0

### DEJU wing: Best window for max temp (6-0)

DEJU_data_use$max_temp_start_w <- DEJU_data_use$measurement_date-6
DEJU_data_use$max_temp_end_w <- DEJU_data_use$measurement_date-0


### DEJU tarsus: Best window for avg temp (7-0 & 30-19, late,early)

DEJU_data_use$avg_temp_late_start_t <- DEJU_data_use$measurement_date-7
DEJU_data_use$avg_temp_late_end_t <- DEJU_data_use$measurement_date-0

DEJU_data_use$avg_temp_early_start_t <- DEJU_data_use$measurement_date-30
DEJU_data_use$avg_temp_early_end_t <- DEJU_data_use$measurement_date-19


### DEJU tarsus: Best window for min temp (28-18)

DEJU_data_use$min_temp_start_t <- DEJU_data_use$measurement_date-28
DEJU_data_use$min_temp_end_t <- DEJU_data_use$measurement_date-18


### DEJU tarsus: Best window for max temp (6-0)

DEJU_data_use$max_temp_start_t <- DEJU_data_use$measurement_date-6
DEJU_data_use$max_temp_end_t <- DEJU_data_use$measurement_date-0


### DEJU mass: Best window for max temp (6-0)

DEJU_data_use$max_temp_start_m <- DEJU_data_use$measurement_date-6
DEJU_data_use$max_temp_end_m <- DEJU_data_use$measurement_date-0



###############
### Sub 5 hrs
###############

### DEJU wing: Best window for < 5 hrs (5-1)

DEJU_data_use$sub_5_start_w <- DEJU_data_use$measurement_date-5
DEJU_data_use$sub_5_end_w <- DEJU_data_use$measurement_date-1


### DEJU tarsus: Best window for < 5 hrs (6-0)

DEJU_data_use$sub_5_start_t <- DEJU_data_use$measurement_date-6
DEJU_data_use$sub_5_end_t <- DEJU_data_use$measurement_date-0


###############
### Sub 10 hrs
###############

### DEJU wing: Best window for sub 10 hrs (7-0)

DEJU_data_use$sub_10_start_w <- DEJU_data_use$measurement_date-7
DEJU_data_use$sub_10_end_w <- DEJU_data_use$measurement_date-0


### DEJU tarsus: Best window for sub 10 hrs (7-0 & 30-17; late, early)

DEJU_data_use$sub_10_early_start_t <- DEJU_data_use$measurement_date-30
DEJU_data_use$sub_10_early_end_t <- DEJU_data_use$measurement_date-17

DEJU_data_use$sub_10_late_start_t <- DEJU_data_use$measurement_date-7
DEJU_data_use$sub_10_late_end_t <- DEJU_data_use$measurement_date-0


###############
### Storms
###############

### DEJU mass: Best window for storms (28-18 & 5-0; early/late)

DEJU_data_use$storms_early_start_m <- DEJU_data_use$measurement_date-28
DEJU_data_use$storms_early_end_m <- DEJU_data_use$measurement_date-18

DEJU_data_use$storms_middle_start_m <- DEJU_data_use$measurement_date-22
DEJU_data_use$storms_middle_end_m <- DEJU_data_use$measurement_date-6

DEJU_data_use$storms_late_start_m <- DEJU_data_use$measurement_date-5
DEJU_data_use$storms_late_end_m <- DEJU_data_use$measurement_date-0


### Overall for summary table

DEJU_data_use$overall_start <- DEJU_data_use$measurement_date-30
DEJU_data_use$overall_end <- DEJU_data_use$measurement_date-0



###########################################
### Convert julian dates to proper posix format
###########################################

### Convert weather dates to POSIXct or else the following code will not work 
### (i.e., it won't work with same date format as sliding window)

DEJU_weather$date.int <- strptime(paste(DEJU_weather$year, DEJU_weather$j_date), format="%Y %j") 
DEJU_weather$date_new <- as.POSIXct(DEJU_weather$date.int)

head(DEJU_weather)

### Remove date.int or else loop below will not work

DEJU_weather <- DEJU_weather[,-12]


###############
### Avg temp
###############


### DEJU wing length
DEJU_data_use$avg_temp_start_w_int <- strptime(paste(DEJU_data_use$year, DEJU_data_use$avg_temp_start_w), format="%Y %j") 
DEJU_data_use$avg_temp_start_w_date <- as.POSIXct(DEJU_data_use$avg_temp_start_w_int)

DEJU_data_use$avg_temp_end_w_int <- strptime(paste(DEJU_data_use$year, DEJU_data_use$avg_temp_end_w), format="%Y %j") 
DEJU_data_use$avg_temp_end_w_date <- as.POSIXct(DEJU_data_use$avg_temp_end_w_int)


### DEJU tarsus length
DEJU_data_use$avg_temp_early_start_t_int <- strptime(paste(DEJU_data_use$year, DEJU_data_use$avg_temp_early_start_t), format="%Y %j") 
DEJU_data_use$avg_temp_early_start_t_date <- as.POSIXct(DEJU_data_use$avg_temp_early_start_t_int)

DEJU_data_use$avg_temp_early_end_t_int <- strptime(paste(DEJU_data_use$year, DEJU_data_use$avg_temp_early_end_t), format="%Y %j") 
DEJU_data_use$avg_temp_early_end_t_date <- as.POSIXct(DEJU_data_use$avg_temp_early_end_t_int)


DEJU_data_use$avg_temp_late_start_t_int <- strptime(paste(DEJU_data_use$year, DEJU_data_use$avg_temp_late_start_t), format="%Y %j") 
DEJU_data_use$avg_temp_late_start_t_date <- as.POSIXct(DEJU_data_use$avg_temp_late_start_t_int)

DEJU_data_use$avg_temp_late_end_t_int <- strptime(paste(DEJU_data_use$year, DEJU_data_use$avg_temp_late_end_t), format="%Y %j") 
DEJU_data_use$avg_temp_late_end_t_date <- as.POSIXct(DEJU_data_use$avg_temp_late_end_t_int)


###############
### Min temp
###############


### DEJU tarsus length
DEJU_data_use$min_temp_start_t_int <- strptime(paste(DEJU_data_use$year, DEJU_data_use$min_temp_start_t), format="%Y %j") 
DEJU_data_use$min_temp_start_t_date <- as.POSIXct(DEJU_data_use$min_temp_start_t_int)

DEJU_data_use$min_temp_end_t_int <- strptime(paste(DEJU_data_use$year, DEJU_data_use$min_temp_end_t), format="%Y %j") 
DEJU_data_use$min_temp_end_t_date <- as.POSIXct(DEJU_data_use$min_temp_end_t_int)


###############
### Max temp
###############

### DEJU wing length
DEJU_data_use$max_temp_start_w_int <- strptime(paste(DEJU_data_use$year, DEJU_data_use$max_temp_start_w), format="%Y %j") 
DEJU_data_use$max_temp_start_w_date <- as.POSIXct(DEJU_data_use$max_temp_start_w_int)

DEJU_data_use$max_temp_end_w_int <- strptime(paste(DEJU_data_use$year, DEJU_data_use$max_temp_end_w), format="%Y %j") 
DEJU_data_use$max_temp_end_w_date <- as.POSIXct(DEJU_data_use$max_temp_end_w_int)


### DEJU tarsus length
DEJU_data_use$max_temp_start_t_int <- strptime(paste(DEJU_data_use$year, DEJU_data_use$max_temp_start_t), format="%Y %j") 
DEJU_data_use$max_temp_start_t_date <- as.POSIXct(DEJU_data_use$max_temp_start_t_int)

DEJU_data_use$max_temp_end_t_int <- strptime(paste(DEJU_data_use$year, DEJU_data_use$max_temp_end_t), format="%Y %j") 
DEJU_data_use$max_temp_end_t_date <- as.POSIXct(DEJU_data_use$max_temp_end_t_int)


### DEJU mass length
DEJU_data_use$max_temp_start_m_int <- strptime(paste(DEJU_data_use$year, DEJU_data_use$max_temp_start_m), format="%Y %j") 
DEJU_data_use$max_temp_start_m_date <- as.POSIXct(DEJU_data_use$max_temp_start_m_int)

DEJU_data_use$max_temp_end_m_int <- strptime(paste(DEJU_data_use$year, DEJU_data_use$max_temp_end_m), format="%Y %j") 
DEJU_data_use$max_temp_end_m_date <- as.POSIXct(DEJU_data_use$max_temp_end_m_int)


###############
### Sub 5 hrs
###############

### DEJU wing length
DEJU_data_use$sub_5_start_w_int <- strptime(paste(DEJU_data_use$year, DEJU_data_use$sub_5_start_w), format="%Y %j") 
DEJU_data_use$sub_5_start_w_date <- as.POSIXct(DEJU_data_use$sub_5_start_w_int)

DEJU_data_use$sub_5_end_w_int <- strptime(paste(DEJU_data_use$year, DEJU_data_use$sub_5_end_w), format="%Y %j") 
DEJU_data_use$sub_5_end_w_date <- as.POSIXct(DEJU_data_use$sub_5_end_w_int)


### DEJU tarsus
DEJU_data_use$sub_5_start_t_int <- strptime(paste(DEJU_data_use$year, DEJU_data_use$sub_5_start_t), format="%Y %j") 
DEJU_data_use$sub_5_start_t_date <- as.POSIXct(DEJU_data_use$sub_5_start_t_int)

DEJU_data_use$sub_5_end_t_int <- strptime(paste(DEJU_data_use$year, DEJU_data_use$sub_5_end_t), format="%Y %j") 
DEJU_data_use$sub_5_end_t_date <- as.POSIXct(DEJU_data_use$sub_5_end_t_int)


###############
### Sub 10 hrs
###############

### DEJU wing length
DEJU_data_use$sub_10_start_w_int <- strptime(paste(DEJU_data_use$year, DEJU_data_use$sub_10_start_w), format="%Y %j") 
DEJU_data_use$sub_10_start_w_date <- as.POSIXct(DEJU_data_use$sub_10_start_w_int)

DEJU_data_use$sub_10_end_w_int <- strptime(paste(DEJU_data_use$year, DEJU_data_use$sub_10_end_w), format="%Y %j") 
DEJU_data_use$sub_10_end_w_date <- as.POSIXct(DEJU_data_use$sub_10_end_w_int)


### DEJU tarsus length
DEJU_data_use$sub_10_early_start_t_int <- strptime(paste(DEJU_data_use$year, DEJU_data_use$sub_10_early_start_t), format="%Y %j") 
DEJU_data_use$sub_10_early_start_t_date <- as.POSIXct(DEJU_data_use$sub_10_early_start_t_int)

DEJU_data_use$sub_10_early_end_t_int <- strptime(paste(DEJU_data_use$year, DEJU_data_use$sub_10_early_end_t), format="%Y %j") 
DEJU_data_use$sub_10_early_end_t_date <- as.POSIXct(DEJU_data_use$sub_10_early_end_t_int)


DEJU_data_use$sub_10_late_start_t_int <- strptime(paste(DEJU_data_use$year, DEJU_data_use$sub_10_late_start_t), format="%Y %j") 
DEJU_data_use$sub_10_late_start_t_date <- as.POSIXct(DEJU_data_use$sub_10_late_start_t_int)

DEJU_data_use$sub_10_late_end_t_int <- strptime(paste(DEJU_data_use$year, DEJU_data_use$sub_10_late_end_t), format="%Y %j") 
DEJU_data_use$sub_10_late_end_t_date <- as.POSIXct(DEJU_data_use$sub_10_late_end_t_int)



###############
### Storms
###############


### DEJU mass
DEJU_data_use$storms_early_start_m_int <- strptime(paste(DEJU_data_use$year, DEJU_data_use$storms_early_start_m), format="%Y %j") 
DEJU_data_use$storms_early_start_m_date <- as.POSIXct(DEJU_data_use$storms_early_start_m_int)

DEJU_data_use$storms_early_end_m_int <- strptime(paste(DEJU_data_use$year, DEJU_data_use$storms_early_end_m), format="%Y %j") 
DEJU_data_use$storms_early_end_m_date <- as.POSIXct(DEJU_data_use$storms_early_end_m_int)

DEJU_data_use$storms_middle_start_m_int <- strptime(paste(DEJU_data_use$year, DEJU_data_use$storms_middle_start_m), format="%Y %j") 
DEJU_data_use$storms_middle_start_m_date <- as.POSIXct(DEJU_data_use$storms_middle_start_m_int)

DEJU_data_use$storms_middle_end_m_int <- strptime(paste(DEJU_data_use$year, DEJU_data_use$storms_middle_end_m), format="%Y %j") 
DEJU_data_use$storms_middle_end_m_date <- as.POSIXct(DEJU_data_use$storms_middle_end_m_int)

DEJU_data_use$storms_late_start_m_int <- strptime(paste(DEJU_data_use$year, DEJU_data_use$storms_late_start_m), format="%Y %j") 
DEJU_data_use$storms_late_start_m_date <- as.POSIXct(DEJU_data_use$storms_late_start_m_int)

DEJU_data_use$storms_late_end_m_int <- strptime(paste(DEJU_data_use$year, DEJU_data_use$storms_late_end_m), format="%Y %j") 
DEJU_data_use$storms_late_end_m_date <- as.POSIXct(DEJU_data_use$storms_late_end_m_int)


### Overall

DEJU_data_use$overall_start_int <- strptime(paste(DEJU_data_use$year, DEJU_data_use$overall_start), format="%Y %j") 
DEJU_data_use$overall_start_date <- as.POSIXct(DEJU_data_use$overall_start_int)

DEJU_data_use$overall_end_int <- strptime(paste(DEJU_data_use$year, DEJU_data_use$overall_end), format="%Y %j") 
DEJU_data_use$overall_end_date <- as.POSIXct(DEJU_data_use$overall_end_int)



#############################################
### Create loop to extract weather variables for each nestling and save as a csv
#############################################

### Empty lists

DEJU_avg_temp_w_list <- list()
DEJU_avg_temp_early_t_list <- list()
DEJU_avg_temp_late_t_list <- list()

DEJU_min_temp_t_list <- list()

DEJU_max_temp_w_list <- list()
DEJU_max_temp_t_list <- list()
DEJU_max_temp_m_list <- list()

DEJU_sub_5_w_list <- list()
DEJU_sub_5_t_list <- list()

DEJU_sub_10_w_list <- list()
DEJU_early_sub_10_t_list <- list()
DEJU_late_sub_10_t_list <- list()

DEJU_early_storms_m_list <- list()

DEJU_late_storms_m_list <- list()

DEJU_middle_storms_m_list <- list()

DEJU_overall_list <- list()


### Loop

for (i in 1:120) {
        
        DEJU_avg_temp_w_list[[i]] <- DEJU_weather %>% 
                filter(date_new >= DEJU_data_use$avg_temp_start_w_date[i] 
                       & date_new <= DEJU_data_use$avg_temp_end_w_date[i]) %>% 
                
                summarise(avg_temp_w = mean(temp_avg_day))
        
        DEJU_avg_temp_early_t_list[[i]] <- DEJU_weather %>% 
                filter(date_new >= DEJU_data_use$avg_temp_early_start_t_date[i] 
                       & date_new <= DEJU_data_use$avg_temp_early_end_t_date[i]) %>% 
                
                summarise(avg_temp_early_t = mean(temp_avg_day))
        
        DEJU_avg_temp_late_t_list[[i]] <- DEJU_weather %>% 
          filter(date_new >= DEJU_data_use$avg_temp_late_start_t_date[i] 
                 & date_new <= DEJU_data_use$avg_temp_late_end_t_date[i]) %>% 
          
          summarise(avg_temp_late_t = mean(temp_avg_day))
        
        
        DEJU_min_temp_t_list[[i]] <- DEJU_weather %>% 
          filter(date_new >= DEJU_data_use$min_temp_start_t_date[i] 
                 & date_new <= DEJU_data_use$min_temp_end_t_date[i]) %>% 
          
          summarise(min_temp_t = min(temp_avg_day))
        
        DEJU_max_temp_w_list[[i]] <- DEJU_weather %>% 
          filter(date_new >= DEJU_data_use$max_temp_start_w_date[i] 
                 & date_new <= DEJU_data_use$max_temp_end_w_date[i]) %>% 
          
          summarise(max_temp_w = max(temp_avg_day))
        
        DEJU_max_temp_t_list[[i]] <- DEJU_weather %>% 
          filter(date_new >= DEJU_data_use$max_temp_start_t_date[i] 
                 & date_new <= DEJU_data_use$max_temp_end_t_date[i]) %>% 
          
          summarise(max_temp_t = max(temp_avg_day))
        
        DEJU_max_temp_m_list[[i]] <- DEJU_weather %>% 
          filter(date_new >= DEJU_data_use$max_temp_start_m_date[i] 
                 & date_new <= DEJU_data_use$max_temp_end_m_date[i]) %>% 
          
          summarise(max_temp_m = max(temp_avg_day))
        
        
        
        DEJU_sub_5_w_list[[i]] <- DEJU_weather %>% 
                filter(date_new >= DEJU_data_use$sub_5_start_w_date[i] 
                       & date_new <= DEJU_data_use$sub_5_end_w_date[i]) %>% 
                
                summarise(sub_5_w = mean(temp_sub_5_hrs))
        
        DEJU_sub_5_t_list[[i]] <- DEJU_weather %>% 
                filter(date_new >= DEJU_data_use$sub_5_start_t_date[i] 
                       & date_new <= DEJU_data_use$sub_5_end_t_date[i]) %>% 
                
                summarise(sub_5_t = mean(temp_sub_5_hrs))
        
       DEJU_sub_10_w_list[[i]] <- DEJU_weather %>% 
                filter(date_new >= DEJU_data_use$sub_10_start_w_date[i] 
                       & date_new <= DEJU_data_use$sub_10_end_w_date[i]) %>% 
                
                summarise(sub_10_w = mean(temp_sub_10_hrs))
        
       DEJU_early_sub_10_t_list[[i]] <- DEJU_weather %>% 
                filter(date_new >= DEJU_data_use$sub_10_early_start_t_date[i] 
                       & date_new <= DEJU_data_use$sub_10_early_end_t_date[i]) %>% 
                
                summarise(sub_10_early_t = mean(temp_sub_10_hrs))
       
       DEJU_late_sub_10_t_list[[i]] <- DEJU_weather %>% 
               filter(date_new >= DEJU_data_use$sub_10_late_start_t_date[i] 
                      & date_new <= DEJU_data_use$sub_10_late_end_t_date[i]) %>% 
               
               summarise(sub_10_late_t = mean(temp_sub_10_hrs))
        
       
       DEJU_early_storms_m_list[[i]] <- DEJU_weather %>% 
               filter(date_new >= DEJU_data_use$storms_early_start_m_date[i] 
                      & date_new <= DEJU_data_use$storms_early_end_m_date[i]) %>% 
               
               summarise(total_storms_early_m = sum(storms))
       
       DEJU_late_storms_m_list[[i]] <- DEJU_weather %>% 
               filter(date_new >= DEJU_data_use$storms_late_start_m_date[i] 
                      & date_new <= DEJU_data_use$storms_late_end_m_date[i]) %>% 
               
               summarise(total_storms_late_m = sum(storms))
       
       DEJU_middle_storms_m_list[[i]] <- DEJU_weather %>% 
                filter(date_new >= DEJU_data_use$storms_middle_start_m_date[i] 
                      & date_new <= DEJU_data_use$storms_middle_end_m_date[i]) %>% 
         
                summarise(total_storms_middle_m = sum(storms))
       
       DEJU_overall_list[[i]] <- DEJU_weather %>% 
                filter(date_new >= DEJU_data_use$overall_start_date[i] 
                      & date_new <= DEJU_data_use$overall_end_date[i]) %>% 
         
                summarise(avg_temp = mean(temp_avg_day),
                          min_temp = min(temp_avg_day),
                          max_temp = max(temp_avg_day),
                          avg_5 = mean(temp_sub_5_hrs),
                          avg_10 = mean(temp_sub_10_hrs),
                          total_precip = sum(precip_days),
                          total_storms = sum(storms))
       
       
        
        ### Convert lists to rows per nestling ID
        
       
       DEJU_avg_temp_w_data <- bind_rows(DEJU_avg_temp_w_list)
       DEJU_avg_temp_early_t_data <- bind_rows(DEJU_avg_temp_early_t_list)
       DEJU_avg_temp_late_t_data <- bind_rows(DEJU_avg_temp_late_t_list)
       
       DEJU_min_temp_t_data <- bind_rows(DEJU_min_temp_t_list)
       
       DEJU_max_temp_w_data <- bind_rows(DEJU_max_temp_w_list)
       DEJU_max_temp_t_data <- bind_rows(DEJU_max_temp_t_list)
       DEJU_max_temp_m_data <- bind_rows(DEJU_max_temp_m_list)
       
       DEJU_sub_5_w_data <- bind_rows(DEJU_sub_5_w_list)
       DEJU_sub_5_t_data <- bind_rows(DEJU_sub_5_t_list)
       
       DEJU_sub_10_w_data <- bind_rows(DEJU_sub_10_w_list)
       DEJU_early_sub_10_t_data <- bind_rows(DEJU_early_sub_10_t_list)
       DEJU_late_sub_10_t_data <- bind_rows(DEJU_late_sub_10_t_list)
       
       DEJU_early_storms_m_data <- bind_rows(DEJU_early_storms_m_list)
       
       DEJU_late_storms_m_data <- bind_rows(DEJU_late_storms_m_list)
       
       DEJU_middle_storms_m_data <- bind_rows(DEJU_middle_storms_m_list)
       
       DEJU_overall_data <- bind_rows(DEJU_overall_list)
       
       
       
       ### Convert each variable to a column in a larger dataset
       DEJU_windows <- bind_cols(DEJU_avg_temp_w_data, DEJU_avg_temp_early_t_data, DEJU_avg_temp_late_t_data,
                                 DEJU_min_temp_t_data,
                                 DEJU_max_temp_w_data, DEJU_max_temp_t_data, DEJU_max_temp_m_data,
                                 DEJU_sub_5_w_data, DEJU_sub_5_t_data,
                                 DEJU_sub_10_w_data, DEJU_early_sub_10_t_data, DEJU_late_sub_10_t_data,
                                 DEJU_early_storms_m_data,
                                 DEJU_late_storms_m_data,
                                 DEJU_middle_storms_m_data, DEJU_overall_data)
        
}

head(DEJU_windows)


### Combine with original measurement data

head(DEJU_data_use)

DEJU_data_new <- DEJU_data_use[, c(1:12)] # Select only original columns of interest

DEJU_final <- bind_cols(DEJU_data_new, DEJU_windows)

### Save data

write.csv(DEJU_final, "DEJU_weather_summarized_final.csv")




############################################################################################
### (3) Savannah Sparrow sliding window ####################################################
############################################################################################


#### Convert measurement date into proper date format

SAVS_data$measured <- strptime(paste(SAVS_data$year, SAVS_data$measurement_date), format="%Y %j") 
SAVS_data$measured_date_formatted <- format(SAVS_data$measured, format= "%d/%m/%Y")

### Remove interim weather variable

head(SAVS_data)

SAVS_data_use <- SAVS_data[,-14]

### Convert year and nest ID to factor

SAVS_data_use$year_f <- factor(SAVS_data_use$year)

SAVS_data_use$nest_id <- factor(SAVS_data_use$nest_id)


### Select only nestlings with either tarsus or mass measurements.

SAVS_tarsus_data <- subset(SAVS_data_use, tarsus >0)

SAVS_mass_data <- subset(SAVS_data_use, mass >0)


### Number of nestlings and nests

SAVS_mass_data$nest_id <- factor(SAVS_mass_data$nest_id)
SAVS_tarsus_data$nest_id <- factor(SAVS_tarsus_data$nest_id)

length(SAVS_tarsus_data$tarsus) # 96 nestlings
length(table(SAVS_tarsus_data$nest_id)) # 25 nests


#############
### Incorporate weather data 
#############

### Convert date to correct format

SAVS_weather$date.new <- strptime(paste(SAVS_weather$year, SAVS_weather$j_date), format="%Y %j") 
SAVS_weather$date_use <- format(SAVS_weather$date.new, format= "%d/%m/%Y")

head(SAVS_weather)

### Remove interim weather variable

SAVS_weather <- SAVS_weather[,-11]


#######################################
#### SAVS sliding window #######
#######################################

### Note: Average clutch size = 5, average incubation = 13, so first egg laid 23 days prior to
### measurement


###############################
### (A) Tarsus length
###############################

### Test if year should be included to control for among year differences

SAVS_wo_year <- lmer(tarsus ~ brood_size + measurement_age + DFE + (1|nest_id), data=SAVS_tarsus_data)

SAVS_year <- lmer(tarsus ~ brood_size + measurement_age + DFE + year_f + (1|nest_id), data=SAVS_tarsus_data)

anova(SAVS_wo_year, SAVS_year, test = "F") # P = 0.42

### Decision: No evidence that year improves the model. Do not include.


##################################
### Sliding windows
##################################

# Temperature variables
SAVS_tarsus_sw_temp <- slidingwin(baseline = lmer(tarsus ~ brood_size + measurement_age + DFE + (1|nest_id), data=SAVS_tarsus_data),  ### this is the base model to which the program adds the climate variables
                        xvar = list(TempMean = SAVS_weather$temp_avg_day,
                                SubZero = SAVS_weather$temp_sub_zero_hrs,
                                SubFive = SAVS_weather$temp_sub_5_hrs, 
                                SubTen  = SAVS_weather$temp_sub_10_hrs), 
                        type = "relative", 
                        range= c(30,0), 
                        stat = c("mean","min","max","var"), 
                        cmissing = "method1",
                        func = "lin",
                        cinterval = "day",
                        cdate = SAVS_weather$date_use, bdate = SAVS_tarsus_data$measured_date_formatted) 

# Temperature results
SAVS_tarsus_sw_temp$combos

head(SAVS_tarsus_sw_temp[[1]]$Dataset, 20)  # Mean avg temp - 30-25 (+)
head(SAVS_tarsus_sw_temp[[2]]$Dataset, 20)  # Mean sub zero 27-19 (+) & 17-0 (-)
head(SAVS_tarsus_sw_temp[[3]]$Dataset, 20)  # Mean sub 5 hrs (no different from null)
head(SAVS_tarsus_sw_temp[[4]]$Dataset, 20)  # Mean sub 10 hrs (no different from null)

head(SAVS_tarsus_sw_temp[[5]]$Dataset, 20) # Min mean temp (29-19) (1.35)
head(SAVS_tarsus_sw_temp[[9]]$Dataset, 20) # Max mean temp (none better than null)

head(SAVS_tarsus_sw_temp[[13]]$Dataset, 20) # Variance avg temp (none better than null)
head(SAVS_tarsus_sw_temp[[14]]$Dataset, 20) # Variance sub zero 27-19 (+) - randomization too rank deficient to fit
head(SAVS_tarsus_sw_temp[[15]]$Dataset, 20) # Varaince sub 5 (none better than null)
head(SAVS_tarsus_sw_temp[[16]]$Dataset, 20) # Variance sub 10 (none better than null)


# Precipitation variables
SAVS_tarsus_sw_precip <- slidingwin(baseline = lmer(tarsus ~ brood_size + measurement_age + DFE + (1|nest_id), data=SAVS_tarsus_data),  ### this is the base model to which the program adds the climate variables
                                  xvar = list(PrecipDays = SAVS_weather$precip_days,
                                              Storms = SAVS_weather$storms), 
                                  type = "relative", 
                                  range= c(30,0), 
                                  stat = c("mean"), 
                                  cmissing = "method1",
                                  func = "lin",
                                  cinterval = "day",
                                  cdate = SAVS_weather$date_use, bdate = SAVS_tarsus_data$measured_date_formatted) 

# Precipitation results
SAVS_tarsus_sw_precip$combos

head(SAVS_tarsus_sw_precip[[1]]$Dataset, 20)  # Precip days 22-15 (+)
head(SAVS_tarsus_sw_precip[[2]]$Dataset, 20)  # Storms 25-14 (+)


##################
### Randomization for SAVS tarsus
##################

# Mean temperature variables
SAVS_tarsus_rand_temp_mean <- randwin(repeats = 100, 
                     xvar = list(TempMean = SAVS_weather$temp_avg_day,
                                 SubZero = SAVS_weather$temp_sub_zero_hrs),
                     cdate = SAVS_weather$date_use, 
                     bdate = SAVS_tarsus_data$measured_date, 
                     baseline = lmer(tarsus ~ brood_size + measurement_age + DFE + (1|nest_id),
                                     data=SAVS_tarsus_data), 
                     range = c(30, 0), 
                     type = "relative", stat = c("mean"), 
                     func = "lin", cinterval = "day")

# Min avg temp
SAVS_tarsus_rand_temp_min <- randwin(repeats = 100, 
                                      xvar = list(TempMean = SAVS_weather$temp_avg_day),
                                      cdate = SAVS_weather$date_use, 
                                      bdate = SAVS_tarsus_data$measured_date, 
                                      baseline = lmer(tarsus ~ brood_size + measurement_age + DFE + (1|nest_id),
                                                      data=SAVS_tarsus_data), 
                                      range = c(30, 0), 
                                      type = "relative", stat = c("min"), 
                                      func = "lin", cinterval = "day")

# Precipitation variables
SAVS_tarsus_rand_precip <- randwin(repeats = 100, 
                                     xvar = list(PrecipDays = SAVS_weather$precip_days,
                                                 Storms = SAVS_weather$storms),
                                     cdate = SAVS_weather$date_use, 
                                     bdate = SAVS_tarsus_data$measured_date, 
                                     baseline = lmer(tarsus ~ brood_size + measurement_age + DFE + (1|nest_id),
                                                     data=SAVS_tarsus_data), 
                                     range = c(30, 0), 
                                     type = "relative", stat = c("mean"), 
                                     func = "lin", cinterval = "day")

###########################
### Compare original and randomization
###########################

length(SAVS_tarsus_data$species) # sample size to use = 96

# Mean avg temp
pvalue(dataset = SAVS_tarsus_sw_temp[[1]]$Dataset, datasetrand = SAVS_tarsus_rand_temp_mean[[1]], 
       metric="AIC", sample.size=96) #  0.06

# Mean sub 0 hrs
pvalue(dataset = SAVS_tarsus_sw_temp[[2]]$Dataset, datasetrand = SAVS_tarsus_rand_temp_mean[[2]], 
       metric="AIC", sample.size=96) # 0.33

# Min avg temp
pvalue(dataset = SAVS_tarsus_sw_temp[[5]]$Dataset, datasetrand = SAVS_tarsus_rand_temp_min[[1]], 
       metric="AIC", sample.size=96) #  0.09

# Mean precip days
pvalue(dataset = SAVS_tarsus_sw_precip[[1]]$Dataset, datasetrand = SAVS_tarsus_rand_precip[[1]], 
       metric="AIC", sample.size=96)  # 0.12

# Mean storm events
pvalue(dataset = SAVS_tarsus_sw_precip[[2]]$Dataset, datasetrand = SAVS_tarsus_rand_precip[[2]], 
       metric="AIC", sample.size=96) # 0.22


###########################
### Conclusion: 
# Retain mean avg temp and min avg temp.
# There is no support for the rest of the variables.



###############################
### (B) Mass
###############################


## Test if year should be included to control for among year differences

SAVS_wo_year_m <- lmer(mass ~ brood_size + measurement_age + DFE + (1|nest_id), data=SAVS_mass_data)

SAVS_year_m <- lmer(mass ~ brood_size + measurement_age + DFE + year_f + (1|nest_id), data=SAVS_mass_data)

anova(SAVS_wo_year_m, SAVS_year_m, test = "F") # P = 0.48

### Decision: No evidence that year improves the model. Do not include.


###########################
### Sliding windows
###########################

# Temperature variables
SAVS_mass_sw_temp <- slidingwin(baseline = lmer(mass ~ brood_size + measurement_age + DFE + (1|nest_id), data=SAVS_mass_data),  ### this is the base model to which the program adds the climate variables
                          xvar = list(TempMean = SAVS_weather$temp_avg_day,
                                      SubZero = SAVS_weather$temp_sub_zero_hrs,
                                      SubFive = SAVS_weather$temp_sub_5_hrs, 
                                      SubTen  = SAVS_weather$temp_sub_10_hrs), 
                          type = "relative", 
                          range= c(30,0), 
                          stat = c("mean","min","max","var"), 
                          cmissing = "method1",
                          func = "lin",
                          cinterval = "day",
                          cdate = SAVS_weather$date_use, bdate = SAVS_mass_data$measured_date_formatted) 

# Temperature results
SAVS_mass_sw_temp$combos

head(SAVS_mass_sw_temp[[1]]$Dataset, 20)  # Mean avg temp - 3-1 (+), 30-23 (+)
head(SAVS_mass_sw_temp[[2]]$Dataset, 20)  # Mean sub zero 27-19 (+)
head(SAVS_mass_sw_temp[[3]]$Dataset, 20)  # Mean sub 5 hrs - 30-23 (-)
head(SAVS_mass_sw_temp[[4]]$Dataset, 20)  # Mean sub 10 hrs 30-28 (-)

head(SAVS_mass_sw_temp[[5]]$Dataset, 20)  # Min avg temp (29-19) (2.39)
head(SAVS_mass_sw_temp[[9]]$Dataset, 20)  # Max avg temp (3-1 & 29-27) (0.84 & 0.48)

head(SAVS_mass_sw_temp[[13]]$Dataset, 20) # Variance avg temp (23-18) (-)
head(SAVS_mass_sw_temp[[14]]$Dataset, 20) # Variance sub zero (27-19) (+) - Randomization unidentifiable
head(SAVS_mass_sw_temp[[15]]$Dataset, 20) # Variance sub 5 (none different from null)
head(SAVS_mass_sw_temp[[16]]$Dataset, 20) # Variance sub 10 (none different from null)



# Precipitation variables
SAVS_mass_sw_precip <- slidingwin(baseline = lmer(mass ~ brood_size + measurement_age + DFE + (1|nest_id), data=SAVS_mass_data),  ### this is the base model to which the program adds the climate variables
                                xvar = list(PrecipDays = SAVS_weather$precip_days,
                                            Storms = SAVS_weather$storms), 
                                type = "relative", 
                                range= c(30,0), 
                                stat = c("mean"), 
                                cmissing = "method1",
                                func = "lin",
                                cinterval = "day",
                                cdate = SAVS_weather$date_use, bdate = SAVS_mass_data$measured_date_formatted) 

# Precipitation results
SAVS_mass_sw_precip$combos

head(SAVS_mass_sw_precip[[1]]$Dataset, 20)  # Precip days 30-12 (+)
head(SAVS_mass_sw_precip[[2]]$Dataset, 20)  # Storms 30-12 (+)



##################
### Randomization for SAVS mass
##################

# Mean temperature variables
SAVS_mass_rand_temp_mean <- randwin(repeats = 100, 
                       xvar = list(TempMean = SAVS_weather$temp_avg_day,
                                   SubZero = SAVS_weather$temp_sub_zero_hrs,
                                   SubFive = SAVS_weather$temp_sub_5_hrs, 
                                   SubTen  = SAVS_weather$temp_sub_10_hrs),
                       cdate = SAVS_weather$date_use, 
                       bdate = SAVS_mass_data$measured_date_formatted, 
                       baseline = lmer(mass ~ brood_size + measurement_age + DFE + (1|nest_id),
                                       data=SAVS_mass_data), 
                       range = c(30, 0), 
                       type = "relative", stat = c("mean"), 
                       func = "lin", cinterval = "day")


# Temperature min and max
SAVS_mass_rand_temp_min_max <- randwin(repeats = 100, 
                                    xvar = list(TempMean = SAVS_weather$temp_avg_day),
                                    cdate = SAVS_weather$date_use, 
                                    bdate = SAVS_mass_data$measured_date_formatted, 
                                    baseline = lmer(mass ~ brood_size + measurement_age + DFE + (1|nest_id),
                                                    data=SAVS_mass_data), 
                                    range = c(30, 0), 
                                    type = "relative", stat = c("min","max"), 
                                    func = "lin", cinterval = "day")

# Temperature variance
SAVS_mass_rand_temp_var <- randwin(repeats = 100, 
                                       xvar = list(TempMean = SAVS_weather$temp_avg_day),
                                       cdate = SAVS_weather$date_use, 
                                       bdate = SAVS_mass_data$measured_date_formatted, 
                                       baseline = lmer(mass ~ brood_size + measurement_age + DFE + (1|nest_id),
                                                       data=SAVS_mass_data), 
                                       range = c(30, 0), 
                                       type = "relative", stat = c("var"), 
                                       func = "lin", cinterval = "day")

# Precipitation variables
SAVS_mass_rand_precip <- randwin(repeats = 100, 
                                       xvar = list(PrecipDays = SAVS_weather$precip_days,
                                                   Storms = SAVS_weather$storms),
                                       cdate = SAVS_weather$date_use, 
                                       bdate = SAVS_mass_data$measured_date_formatted, 
                                       baseline = lmer(mass ~ brood_size + measurement_age + DFE + (1|nest_id),
                                                       data=SAVS_mass_data), 
                                       range = c(30, 0), 
                                       type = "relative", stat = c("mean"), 
                                       func = "lin", cinterval = "day")

#############################
### Compare original and randomization
#############################

length(SAVS_mass_data$species) # sample size used = 86

# Mean avg temp
pvalue(dataset = SAVS_mass_sw_temp[[1]]$Dataset, datasetrand = SAVS_mass_rand_temp_mean[[1]], 
       metric="AIC", sample.size=86) # 0.07

# Mean sub zero
pvalue(dataset = SAVS_mass_sw_temp[[2]]$Dataset, datasetrand = SAVS_mass_rand_temp_mean[[2]], 
       metric="AIC", sample.size=86) # 0.17  

# Mean sub 5
pvalue(dataset = SAVS_mass_sw_temp[[3]]$Dataset, datasetrand = SAVS_mass_rand_temp_mean[[3]], 
       metric="AIC", sample.size=86)  # 0.09

# Mean sub 10
pvalue(dataset = SAVS_mass_sw_temp[[4]]$Dataset, datasetrand = SAVS_mass_rand_temp_mean[[4]], 
       metric="AIC", sample.size=86) # 0.18

# Min avg temp
pvalue(dataset = SAVS_mass_sw_temp[[5]]$Dataset, datasetrand = SAVS_mass_rand_temp_min_max[[1]], 
       metric="AIC", sample.size=86) #  0.03

# Max avg temp
pvalue(dataset = SAVS_mass_sw_temp[[9]]$Dataset, datasetrand = SAVS_mass_rand_temp_min_max[[2]], 
       metric="AIC", sample.size=86) #  0.09

# Var avg temp
pvalue(dataset = SAVS_mass_sw_temp[[13]]$Dataset, datasetrand = SAVS_mass_rand_temp_var[[1]], 
       metric="AIC", sample.size=86) #  0.13

# Mean precip days
pvalue(dataset = SAVS_mass_sw_precip[[1]]$Dataset, datasetrand = SAVS_mass_rand_precip[[1]], 
       metric="AIC", sample.size=86) # <0.001

# Mean storm events
pvalue(dataset = SAVS_mass_sw_precip[[2]]$Dataset, datasetrand = SAVS_mass_rand_precip[[2]], 
       metric="AIC", sample.size=86) # 0.12


##########################
### Conclusion: 
### Retain mean avg temp, mean sub 5, min avg temp, max avg temp, and meanprecip days.



###########################################################################################
### SAVS: Extract weather variables from time periods identified by sliding window ########
###########################################################################################


###############
### Avg temp
###############

### SAVS tarsus: Best window for avg temp (30-25)

SAVS_tarsus_data$avg_temp_start_t <- SAVS_tarsus_data$measurement_date - 30
SAVS_tarsus_data$avg_temp_end_t <- SAVS_tarsus_data$measurement_date - 25

### SAVS tarsus: Best window for min temp (29-19)

SAVS_tarsus_data$min_temp_start_t <- SAVS_tarsus_data$measurement_date - 29
SAVS_tarsus_data$min_temp_end_t <- SAVS_tarsus_data$measurement_date - 19

### SAVS mass: Best window for avg temp (3-1 & 30-23; late, early)

SAVS_mass_data$avg_temp_late_start_m <- SAVS_mass_data$measurement_date - 3
SAVS_mass_data$avg_temp_late_end_m <- SAVS_mass_data$measurement_date - 1

SAVS_mass_data$avg_temp_early_start_m <- SAVS_mass_data$measurement_date - 30
SAVS_mass_data$avg_temp_early_end_m <- SAVS_mass_data$measurement_date - 23


### SAVS mass: Best window for min temp (29-19)

SAVS_mass_data$min_temp_start_m <- SAVS_mass_data$measurement_date - 29
SAVS_mass_data$min_temp_end_m <- SAVS_mass_data$measurement_date - 19


### SAVS mass: Best window for max temp (3-1 & 29-27)

SAVS_mass_data$max_temp_early_start_m <- SAVS_mass_data$measurement_date - 29
SAVS_mass_data$max_temp_early_end_m <- SAVS_mass_data$measurement_date - 27

SAVS_mass_data$max_temp_late_start_m <- SAVS_mass_data$measurement_date - 3
SAVS_mass_data$max_temp_late_end_m <- SAVS_mass_data$measurement_date - 1



###############
### Sub 5 hrs
###############

### SAVS mass: Best window for sub 5 hrs (30-23)

SAVS_mass_data$sub_5_start_m <- SAVS_mass_data$measurement_date-30
SAVS_mass_data$sub_5_end_m <- SAVS_mass_data$measurement_date-23


###############
### Precip days
###############


### SAVS mass: Best window for precip days (22-15)

SAVS_tarsus_data$precip_start_t <- SAVS_tarsus_data$measurement_date-22
SAVS_tarsus_data$precip_end_t <- SAVS_tarsus_data$measurement_date-15


### SAVS mass: Best window for precip days (30-12)

SAVS_mass_data$precip_start_m <- SAVS_mass_data$measurement_date-30
SAVS_mass_data$precip_end_m <- SAVS_mass_data$measurement_date-12


##############
### Overall for summary table
##############

SAVS_tarsus_data$overall_start <- SAVS_tarsus_data$measurement_date-30
SAVS_tarsus_data$overall_end <- SAVS_tarsus_data$measurement_date-0


###########################################
### Convert julian dates to proper posix format
###########################################

### Convert weather dates to POSIXct or else the following code will not work 
### (i.e., it won't work with same date format as sliding window)

SAVS_weather$date.int <- strptime(paste(SAVS_weather$year, SAVS_weather$j_date), format="%Y %j") 
SAVS_weather$date_new <- as.POSIXct(SAVS_weather$date.int)

head(SAVS_weather)

### Remove interim weather variable

SAVS_weather <- SAVS_weather[,-c(12)]


###############
### Avg temp
###############


### SAVS tarsus length
SAVS_tarsus_data$avg_temp_start_t_int <- strptime(paste(SAVS_tarsus_data$year, SAVS_tarsus_data$avg_temp_start_t), format="%Y %j") 
SAVS_tarsus_data$avg_temp_start_t_date <- as.POSIXct(SAVS_tarsus_data$avg_temp_start_t_int)

SAVS_tarsus_data$avg_temp_end_t_int <- strptime(paste(SAVS_tarsus_data$year, SAVS_tarsus_data$avg_temp_end_t), format="%Y %j") 
SAVS_tarsus_data$avg_temp_end_t_date <- as.POSIXct(SAVS_tarsus_data$avg_temp_end_t_int)


### SAVS mass
SAVS_mass_data$avg_temp_early_start_m_int <- strptime(paste(SAVS_mass_data$year, SAVS_mass_data$avg_temp_early_start_m), format="%Y %j") 
SAVS_mass_data$avg_temp_early_start_m_date <- as.POSIXct(SAVS_mass_data$avg_temp_early_start_m_int)

SAVS_mass_data$avg_temp_early_end_m_int <- strptime(paste(SAVS_mass_data$year, SAVS_mass_data$avg_temp_early_end_m), format="%Y %j") 
SAVS_mass_data$avg_temp_early_end_m_date <- as.POSIXct(SAVS_mass_data$avg_temp_early_end_m_int)


SAVS_mass_data$avg_temp_late_start_m_int <- strptime(paste(SAVS_mass_data$year, SAVS_mass_data$avg_temp_late_start_m), format="%Y %j") 
SAVS_mass_data$avg_temp_late_start_m_date <- as.POSIXct(SAVS_mass_data$avg_temp_late_start_m_int)

SAVS_mass_data$avg_temp_late_end_m_int <- strptime(paste(SAVS_mass_data$year, SAVS_mass_data$avg_temp_late_end_m), format="%Y %j") 
SAVS_mass_data$avg_temp_late_end_m_date <- as.POSIXct(SAVS_mass_data$avg_temp_late_end_m_int)



###############
### Min temp
###############


### SAVS tarsus length
SAVS_tarsus_data$min_temp_start_t_int <- strptime(paste(SAVS_tarsus_data$year, SAVS_tarsus_data$min_temp_start_t), format="%Y %j") 
SAVS_tarsus_data$min_temp_start_t_date <- as.POSIXct(SAVS_tarsus_data$min_temp_start_t_int)

SAVS_tarsus_data$min_temp_end_t_int <- strptime(paste(SAVS_tarsus_data$year, SAVS_tarsus_data$min_temp_end_t), format="%Y %j") 
SAVS_tarsus_data$min_temp_end_t_date <- as.POSIXct(SAVS_tarsus_data$min_temp_end_t_int)


### SAVS mass

SAVS_mass_data$min_temp_start_m_int <- strptime(paste(SAVS_mass_data$year, SAVS_mass_data$min_temp_start_m), format="%Y %j") 
SAVS_mass_data$min_temp_start_m_date <- as.POSIXct(SAVS_mass_data$min_temp_start_m_int)

SAVS_mass_data$min_temp_end_m_int <- strptime(paste(SAVS_mass_data$year, SAVS_mass_data$min_temp_end_m), format="%Y %j") 
SAVS_mass_data$min_temp_end_m_date <- as.POSIXct(SAVS_mass_data$min_temp_end_m_int)




###############
### Max temp
###############


SAVS_mass_data$max_temp_early_start_m_int <- strptime(paste(SAVS_mass_data$year, SAVS_mass_data$max_temp_early_start_m), format="%Y %j") 
SAVS_mass_data$max_temp_early_start_m_date <- as.POSIXct(SAVS_mass_data$max_temp_early_start_m_int)

SAVS_mass_data$max_temp_early_end_m_int <- strptime(paste(SAVS_mass_data$year, SAVS_mass_data$max_temp_early_end_m), format="%Y %j") 
SAVS_mass_data$max_temp_early_end_m_date <- as.POSIXct(SAVS_mass_data$max_temp_early_end_m_int)


SAVS_mass_data$max_temp_late_start_m_int <- strptime(paste(SAVS_mass_data$year, SAVS_mass_data$max_temp_late_start_m), format="%Y %j") 
SAVS_mass_data$max_temp_late_start_m_date <- as.POSIXct(SAVS_mass_data$max_temp_late_start_m_int)

SAVS_mass_data$max_temp_late_end_m_int <- strptime(paste(SAVS_mass_data$year, SAVS_mass_data$max_temp_late_end_m), format="%Y %j") 
SAVS_mass_data$max_temp_late_end_m_date <- as.POSIXct(SAVS_mass_data$max_temp_late_end_m_int)



###############
### Sub 5 hrs
###############

### SAVS mass
SAVS_mass_data$sub_5_start_m_int <- strptime(paste(SAVS_mass_data$year, SAVS_mass_data$sub_5_start_m), format="%Y %j") 
SAVS_mass_data$sub_5_start_m_date <- as.POSIXct(SAVS_mass_data$sub_5_start_m_int)

SAVS_mass_data$sub_5_end_m_int <- strptime(paste(SAVS_mass_data$year, SAVS_mass_data$sub_5_end_m), format="%Y %j") 
SAVS_mass_data$sub_5_end_m_date <- as.POSIXct(SAVS_mass_data$sub_5_end_m_int)


###############
### Precip days
###############

### SAVS tarsus
SAVS_tarsus_data$precip_start_t_int <- strptime(paste(SAVS_tarsus_data$year, SAVS_tarsus_data$precip_start_t), format="%Y %j") 
SAVS_tarsus_data$precip_start_t_date <- as.POSIXct(SAVS_tarsus_data$precip_start_t_int)

SAVS_tarsus_data$precip_end_t_int <- strptime(paste(SAVS_tarsus_data$year, SAVS_tarsus_data$precip_end_t), format="%Y %j") 
SAVS_tarsus_data$precip_end_t_date <- as.POSIXct(SAVS_tarsus_data$precip_end_t_int)



### SAVS mass
SAVS_mass_data$precip_start_m_int <- strptime(paste(SAVS_mass_data$year, SAVS_mass_data$precip_start_m), format="%Y %j") 
SAVS_mass_data$precip_start_m_date <- as.POSIXct(SAVS_mass_data$precip_start_m_int)

SAVS_mass_data$precip_end_m_int <- strptime(paste(SAVS_mass_data$year, SAVS_mass_data$precip_end_m), format="%Y %j") 
SAVS_mass_data$precip_end_m_date <- as.POSIXct(SAVS_mass_data$precip_end_m_int)


##############
### Overall
##############

SAVS_tarsus_data$overall_start_int <- strptime(paste(SAVS_tarsus_data$year, SAVS_tarsus_data$overall_start), format="%Y %j") 
SAVS_tarsus_data$overall_start_date <- as.POSIXct(SAVS_tarsus_data$overall_start_int)

SAVS_tarsus_data$overall_end_int <- strptime(paste(SAVS_tarsus_data$year, SAVS_tarsus_data$overall_end), format="%Y %j") 
SAVS_tarsus_data$overall_end_date <- as.POSIXct(SAVS_tarsus_data$overall_end_int)



#############################################
### Create loop to extract weather variables for each nestling and save as a csv
#############################################

### Empty lists

SAVS_avg_temp_t_list <- list()
SAVS_min_temp_t_list <- list()

SAVS_precip_t_list <- list()

SAVS_overall_list <- list()

##########
### Loop - need to create two (tarsus and mass), because different sample sizes
##########


### tarsus

for (i in 1:96) {
  
  SAVS_avg_temp_t_list[[i]] <- SAVS_weather %>% 
            filter(date_new >= SAVS_tarsus_data$avg_temp_start_t_date[i] 
                  & date_new <= SAVS_tarsus_data$avg_temp_end_t_date[i]) %>% 
    
            summarise(avg_temp_t = mean(temp_avg_day))
  
  
  SAVS_min_temp_t_list[[i]] <- SAVS_weather %>% 
    filter(date_new >= SAVS_tarsus_data$min_temp_start_t_date[i] 
           & date_new <= SAVS_tarsus_data$min_temp_end_t_date[i]) %>% 
    
    summarise(min_temp_t = min(temp_avg_day))
  
  SAVS_precip_t_list[[i]] <- SAVS_weather %>% 
            filter(date_new >= SAVS_tarsus_data$precip_start_t_date[i] 
                  & date_new <= SAVS_tarsus_data$precip_end_t_date[i]) %>% 
    
            summarise(total_precip_t = sum(precip_days))
  
  SAVS_overall_list[[i]] <- SAVS_weather %>% 
            filter(date_new >= SAVS_tarsus_data$overall_start_date[i] 
                  & date_new <= SAVS_tarsus_data$overall_end_date[i]) %>% 
    
            summarise(avg_temp = mean(temp_avg_day),
                      min_temp = min(temp_avg_day),
                      max_temp = max(temp_avg_day),
                      avg_5 = mean(temp_sub_5_hrs),
                      avg_10 = mean(temp_sub_10_hrs),
                      total_precip = sum(precip_days),
                      total_storms = sum(storms))
  
}



### Convert list to row per nestling ID
  
SAVS_avg_temp_t_data <- bind_rows(SAVS_avg_temp_t_list)
SAVS_min_temp_t_data <- bind_rows(SAVS_min_temp_t_list)
SAVS_precip_t_data <- bind_rows(SAVS_precip_t_list)
SAVS_overall_data <- bind_rows(SAVS_overall_list)
  
### Create a nestling ID column
  
nestling_id_t <- as.data.frame(SAVS_tarsus_data$nestling_id)
colnames(nestling_id_t) <- "nestling_id"
  
  
### Convert each variable to a column in a larger dataset
SAVS_tarsus_windows <- bind_cols(nestling_id_t, SAVS_avg_temp_t_data, 
                                   SAVS_min_temp_t_data,
                                   SAVS_precip_t_data,
                                   SAVS_overall_data)
head(SAVS_tarsus_windows)



##################  
### Mass
##################

SAVS_avg_temp_early_m_list <- list()
SAVS_avg_temp_late_m_list <- list()
SAVS_min_temp_m_list <- list()
SAVS_max_temp_early_m_list <- list()
SAVS_max_temp_late_m_list <- list()
SAVS_sub_5_m_list <- list()
SAVS_precip_m_list <- list()


for (i in 1:86) {
  
  SAVS_avg_temp_early_m_list[[i]] <- SAVS_weather %>% 
    filter(date_new >= SAVS_mass_data$avg_temp_early_start_m_date[i] 
           & date_new <= SAVS_mass_data$avg_temp_early_end_m_date[i]) %>% 
    
    summarise(avg_temp_early_m = mean(temp_avg_day))
  
  SAVS_avg_temp_late_m_list[[i]] <- SAVS_weather %>% 
    filter(date_new >= SAVS_mass_data$avg_temp_late_start_m_date[i] 
           & date_new <= SAVS_mass_data$avg_temp_late_end_m_date[i]) %>% 
    
    summarise(avg_temp_late_m = mean(temp_avg_day))
  
  
  SAVS_min_temp_m_list[[i]] <- SAVS_weather %>% 
    filter(date_new >= SAVS_mass_data$min_temp_start_m_date[i] 
           & date_new <= SAVS_mass_data$min_temp_end_m_date[i]) %>% 
    
    summarise(min_temp_m = min(temp_avg_day))
  
  
  SAVS_max_temp_early_m_list[[i]] <- SAVS_weather %>% 
    filter(date_new >= SAVS_mass_data$max_temp_early_start_m_date[i] 
           & date_new <= SAVS_mass_data$max_temp_early_end_m_date[i]) %>% 
    
    summarise(max_temp_early_m = max(temp_avg_day))
  
  SAVS_max_temp_late_m_list[[i]] <- SAVS_weather %>% 
    filter(date_new >= SAVS_mass_data$max_temp_late_start_m_date[i] 
           & date_new <= SAVS_mass_data$max_temp_late_end_m_date[i]) %>% 
    
    summarise(max_temp_late_m = max(temp_avg_day))
  
  
  SAVS_sub_5_m_list[[i]] <- SAVS_weather %>% 
    filter(date_new >= SAVS_mass_data$sub_5_start_m_date[i] 
           & date_new <= SAVS_mass_data$sub_5_end_m_date[i]) %>% 
    
    summarise(sub_5_m = mean(temp_sub_5_hrs))
  
  
  SAVS_precip_m_list[[i]] <- SAVS_weather %>% 
    filter(date_new >= SAVS_mass_data$precip_start_m_date[i] 
           & date_new <= SAVS_mass_data$precip_end_m_date[i]) %>% 
    
    summarise(total_precip_m = sum(precip_days))
  
  
  ### Convert list to row per nestling ID
  
  SAVS_avg_temp_early_m_data <- bind_rows(SAVS_avg_temp_early_m_list)
  SAVS_avg_temp_late_m_data <- bind_rows(SAVS_avg_temp_late_m_list)
  SAVS_min_temp_m_data <- bind_rows(SAVS_min_temp_m_list)
  SAVS_max_temp_early_m_data <- bind_rows(SAVS_max_temp_early_m_list)
  SAVS_max_temp_late_m_data <- bind_rows(SAVS_max_temp_late_m_list)
  SAVS_sub_5_m_data <- bind_rows(SAVS_sub_5_m_list)
  SAVS_precip_m_data <- bind_rows(SAVS_precip_m_list)
  
}


### Create vector of nestling ID 

nestling_id_m <- as.data.frame(SAVS_mass_data$nestling_id)
  
colnames(nestling_id_m) <- "nestling_id"
  
### Convert each variable to a column in a larger dataset
SAVS_mass_windows <- bind_cols(nestling_id_m, 
                                 SAVS_avg_temp_early_m_data, SAVS_avg_temp_late_m_data,
                                 SAVS_min_temp_m_data,
                                 SAVS_max_temp_early_m_data, SAVS_max_temp_late_m_data,
                                 SAVS_sub_5_m_data, SAVS_precip_m_data)
  
head(SAVS_mass_windows)


#########
### Combine two SAVS outputs into one data frame
#########


SAVS_windows <- merge(SAVS_tarsus_windows, SAVS_mass_windows, by = "nestling_id", all.x=TRUE)

head(SAVS_windows)

### Combine original measurement data

head(SAVS_data)

SAVS_new <- SAVS_data[, c(1:4,6:9,11,12)]


SAVS_final <- merge(SAVS_new, SAVS_windows, by = "nestling_id", all.x=TRUE)

### Save data

write.csv(SAVS_final, "SAVS_weather_summarized_final.csv")



############################################################################################
########################## End of sliding window code ######################################
############################################################################################