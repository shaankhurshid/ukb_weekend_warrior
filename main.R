### Script to perform main weekend warrior activity analyses
### Assumes availability of certain data as flat files (please see README)
### Assumes availability of specified fields as flat files (please see README)

# Dependencies
library(data.table)
library(survival)
library(prodlim)
library(ggplot2)
library(reshape2)
library(plyr)
library(Cairo)
library(stringr)
library(Epi)
source(file=paste0(getwd(),'/functions/functions.R'))

# Load acceleration files
acceleration <- fread('/accel/ukbb_accel2.csv') # 7089
load(file='/accel/accel_raw_cat.RData') # 7089
withdrawals <- fread(file='/withdrawals/w7089_20230425.csv')

################################################################ Step 1: Basic variables
# Date cleanup
for (j in (c('start_date','end_date'))){set(acceleration,j=j,value=strtrim(acceleration[[j]],width=10))}
for (j in (c('start_date','end_date'))){set(acceleration,j=j,value=as.Date(acceleration[[j]],format='%Y-%m-%d'))}

# Join to get secondary variables, including acceleration epoch availability
setkey(acceleration); setkey(dt_all,sek_id)

# Total time
acceleration[dt_all,':='(mvpa=i.mvpa,sed=i.sed,total=i.total)]

# Sedentary fraction
acceleration[,':='(sed_prop = (sed/total))]

# Standardized sedentary proportion
acceleration[,':='(sed_std = (sed_prop - mean(sed_prop))/sd(sed_prop))]

# Coefficient of variation
acceleration[,':='(cov_10 = (accel_sd/value)*10)]

# Regressed energy expenditure
acceleration[,':='(ee = ifelse((-10.58 + 1.1176*value + 2.9418*sqrt(value) - 0.00059277*value**2) > 0,
                               (-10.58 + 1.1176*value + 2.9418*sqrt(value) - 0.00059277*value**2)/71.225,0))]

##################### Step 2: Derived acceleration variables
# Load derived accel file 
accel_deriv <- fread(file='/weekend_warrior/derived_accel/ukb672918.csv',
                     select=c('eid','40045-0.0','40049-0.0','40037-0.0','40041-0.0','40033-0.0'))

# Friendly renames
setnames(accel_deriv,c('eid','40045-0.0','40049-0.0','40037-0.0','40041-0.0','40033-0.0'),
         c('sample_id','mvpa_daily','mvpa_overall','mvpa_weekday','mvpa_weekend','mvpa_daily_hour'))

# Remove blank MVPA
accel_deriv <- accel_deriv[!is.na(mvpa_overall)]

# Remove blank MVPA daily
accel_deriv2 <- accel_deriv[mvpa_daily != '']

### Link and join daily MVPA
# Extract daily MVPA
mvpa_daily <- data.table(do.call(rbind,sapply(accel_deriv2$mvpa_daily,strsplit,split=',')))
mvpa_daily <- cbind(accel_deriv2$sample_id,mvpa_daily)
names(mvpa_daily) <- c('sample_id','mvpa1','mvpa2','mvpa3','mvpa4','mvpa5','mvpa6','mvpa7')
for (j in names(mvpa_daily)){set(mvpa_daily,j=j,value=as.numeric(mvpa_daily[[j]]))}

# Link
linker <- fread('/accel/7089_17488_linker.csv')

# Join
setkey(linker,app17488); setkey(mvpa_daily,sample_id)
mvpa_daily[linker,sek_id := i.app7089]

setkey(mvpa_daily,sek_id); setkey(acceleration,sample_id)
acceleration[mvpa_daily,':='(mvpa1=i.mvpa1,mvpa2=i.mvpa2,mvpa3=i.mvpa3,mvpa4=i.mvpa4,
                             mvpa5=i.mvpa5,mvpa6=i.mvpa6,mvpa7=i.mvpa7)]

### Link and join overall MVPA
setkey(linker,app17488); setkey(accel_deriv,sample_id)
accel_deriv[linker,sek_id := i.app7089]
setkey(accel_deriv,sek_id); setkey(acceleration,sample_id)
acceleration[accel_deriv,':='(mvpa_overall_average = i.mvpa_overall)]

# Wear time
wear_time1 <- fread(file='/weekend_warrior/monday_wear_duration.csv')
wear_time2 <- fread(file='/weekend_warrior/tuesday_wear_duration.csv')
wear_time3 <- fread(file='/weekend_warrior/wednesday_wear_duration.csv')
wear_time4 <- fread(file='/weekend_warrior/thursday_wear_duration.csv')
wear_time5 <- fread(file='/weekend_warrior/friday_wear_duration.csv')
wear_time6 <- fread(file='/weekend_warrior/saturday_wear_duration.csv')
wear_time7 <- fread(file='/weekend_warrior/sunday_wear_duration.csv')

setkey(wear_time1,sample_id); setkey(wear_time2,sample_id); setkey(wear_time3,sample_id); setkey(acceleration,sample_id)
setkey(wear_time4,sample_id); setkey(wear_time5,sample_id); setkey(wear_time6,sample_id); setkey(wear_time7,sample_id)

acceleration[wear_time1,wear_time1 := i.value]
acceleration[wear_time2,wear_time2 := i.value]
acceleration[wear_time3,wear_time3 := i.value]
acceleration[wear_time4,wear_time4 := i.value]
acceleration[wear_time5,wear_time5 := i.value]
acceleration[wear_time6,wear_time6 := i.value]
acceleration[wear_time7,wear_time7 := i.value]

##################### Step 3: Accelerometer QC and Sample Selection
# Recommended QC for accelerometer data
## Remove revoked consent
acceleration <- acceleration[!(sample_id %in% withdrawals$V1)] # N=103691 - 30 = 103661
## Could not calibrate
no_calibrate <- acceleration[calibration==1] #N=103661 - 11 = 103650
## Inadequate wear time
no_wear_time <- no_calibrate[wear_time==1] #N=103650 - 6985 = 96665
## Remove unrealistically high mean accel
acceleration <- no_wear_time[value <= 100] #N=96665 - 13 = 96652
## MVPA not calculated
acceleration <- acceleration[!c(is.na(mvpa1) | is.na(mvpa2) | is.na(mvpa3) |
                                  is.na(mvpa4) | is.na(mvpa5) | is.na(mvpa6) | is.na(mvpa7))] #96652 - 1159 = 95493 

##################### Step 4: Covariates
### BMI
bmi <- fread(file='/phenotypes/bmi_instance01_202109.csv')

## Smoking
tob0 <- fread('/phenotypes/tobacco_instance_0.csv')
tob1 <- fread('/phenotypes/tobacco_instance1.csv')
setkey(tob0,sample_id); setkey(tob1,sample_id)
tob <- tob0[tob1,':='(value1 = i.value)]
tob[,value_unified := ifelse(!is.na(value1),value1,value)]
tob[,tobacco_accel_selfreport := ifelse(value_unified==2,"Current",
                                        ifelse(value_unified==1,"Former","Never"))]

## TDI
tdi <- fread(file='/phenotypes/townsend_0.csv')

## EtOH
etoh <- fread(file='/phenotypes/etoh_01.csv')

## BP
sbp <- fread('/phenotypes/sbp_combined_instance01.csv')
dbp <- fread('/phenotypes/dbp_combined_instance01.csv')

## anti-HTN
bpmed <- fread('/phenotypes/bpmed_combined_instance01.csv')

## Education
edu <- fread('/phenotypes/education_01.csv')

## Diet
diet <- fread('/phenotypes/diet_instance01.csv')

## Employment status
employment0 <- fread(file='/phenotypes/employment_status_0.csv')
employment0_corrected <- fread(file='/phenotypes/employment_status_corrected.csv')
setnames(employment0_corrected,'20119-0.0','employment_corrected')
employment0_corrected <- employment0_corrected[!is.na(employment_corrected)]

setkey(linker,app17488); setkey(employment0_corrected,eid)
employment0_corrected[linker,app7089 := i.app7089]
setkey(employment0_corrected,app7089); setkey(employment0,sample_id)
employment0[employment0_corrected,value_corrected := i.employment_corrected]
employment0[,value := ifelse(!is.na(value_corrected),value_corrected,value)]

employment1 <- fread(file='/phenotypes/employment_status_1.csv')
employment0 <- employment0[value > 0]
employment1 <- employment1[value > 0]
employment0 <- employment0[,.SD[which.min(value)],by='sample_id']
employment1 <- employment1[,.SD[which.min(value)],by='sample_id']
setkey(employment0,sample_id); setkey(employment1,sample_id)
employment0[employment1,value1 := i.value]
employment0[,value_unified := ifelse(!is.na(value1),value1,value)]

employment0[,employment_status := ifelse(value_unified==1,"Employed","Unemployed/Retired")]

# Self-reported health
self_health0 <- fread(file='/phenotypes/self_health_0.csv')
self_health1 <- fread(file='/phenotypes/self_health_1.csv')
setkey(self_health0,sample_id); setkey(self_health1,sample_id)
self_health0[self_health1,':='(value1 = i.value)]
self_health0[,value_unified := ifelse(!is.na(value1),value1,value)]
self_health0[,self_health := ifelse(value_unified==1,"Excellent",
                                    ifelse(value_unified==2,"Good",
                                           ifelse(value_unified==3,"Fair",
                                                  ifelse(value_unified==4,"Poor",NA))))]

## Race
race <- fread('/phenotypes/race_0.csv')
race[,race_category := ifelse(value %in% c(1,1001,1002,1003),'white',
                              ifelse(value %in% c(2,2001,2002,2003,2004),'mixed',
                                     ifelse(value %in% c(3,3001,3002,3003,3004,5),'asian',
                                            ifelse(value %in% c(4,4001,4002,4003),'black',
                                                   ifelse(value %in% 6,'other',NA)))))]
race[,race_category_adjust := ifelse(race_category=='mixed','other',race_category)]

# Merges
setkey(bmi,sample_id); setkey(acceleration,sample_id); setkey(bpmed,sample_id); setkey(tob,sample_id)
setkey(tdi,sample_id); setkey(dbp,sample_id); setkey(sbp,sample_id); setkey(etoh,sample_id)
setkey(diet,sample_id); setkey(edu,sample_id); setkey(race,sample_id); setkey(self_health0,sample_id);
setkey(employment0,sample_id)

acceleration[bmi, bmi := i.bmi]
acceleration[bpmed, bpmed := i.prev_bpmed]
acceleration[tob, tob := i.tobacco_accel_selfreport]
acceleration[tdi, tdi := i.value]
acceleration[dbp, dbp := i.dbp_combined]
acceleration[sbp, sbp := i.sbp_combined]
acceleration[diet, diet := i.diet_quality]
acceleration[edu, qual_ea := i.qual_ea]
acceleration[etoh, ':='(etoh_grams = i.etoh_grams, etoh_status = i.f_1558)]
acceleration[race, ':='(race_category = i.race_category, race_category_adjust = i.race_category_adjust)]
acceleration[self_health0, self_health := i.self_health]
acceleration[employment0, employment_status := i.employment_status]

no_qual <- acceleration[!is.na(qual_ea)] # 95493 - 831 = 94662
no_diet <- no_qual[!is.na(diet)] # 94662 - 51 = 94611
no_tdi <- no_diet[!is.na(tdi)] # 94611 - 107 = 94504
no_tob <- no_tdi[!is.na(tob)] # 94504 - 2 = 94502
no_race <- no_tob[!is.na(race_category_adjust)] # 94502 - 272 = 94230
no_self_health <- no_race[!is.na(self_health)] # 94230 - 168 = 94062
no_employment <- no_self_health[!is.na(employment_status)] # 94062 - 533 = 93529

# Assume median alcohol intake at listed frequency when not specified
no_employment <- no_employment[,etoh_grams := ifelse(etoh_grams > 2000, NA, etoh_grams)]

monthly_mean <- round(median(no_employment[c(!is.na(etoh_grams) & etoh_status=='Monthly')]$etoh_grams),0)
weekly_mean <- round(median(no_employment[c(!is.na(etoh_grams) & etoh_status=='Weekly')]$etoh_grams),0)

no_employment[c(!is.na(etoh_status) & etoh_status=='Weekly' & is.na(etoh_grams)),
              etoh_grams := (weekly_mean)]
no_employment[c(!is.na(etoh_status) & etoh_status=='Monthly' & is.na(etoh_grams)),
              etoh_grams := (monthly_mean)]

# Remove individuals with missing covars
no_etoh <- no_employment[!is.na(etoh_grams)] # 93530 - 27 = 93502

# Complete acceleration data only
no_etoh[,':='(complete_time = ifelse(total >= 604800,1,0))]
complete_time <- no_etoh[complete_time==1] # 93503 - 3925 = 89577

# MVPA daily totals
complete_time[,":="(mvpa1=mvpa1*1440,mvpa2=mvpa2*1440,mvpa3=mvpa3*1440,
                    mvpa4=mvpa4*1440,mvpa5=mvpa5*1440,
                    mvpa6=mvpa6*1440,mvpa7=mvpa7*1440)]

# Total wear time
complete_time[,":="(total_wear_time = wear_time1 + wear_time2 + wear_time3 + wear_time4 +
                      wear_time5 + wear_time6 + wear_time7)]

# Total MVPA
complete_time[,mvpa_daily_total_per_average := (mvpa_overall_average)*1440]
complete_time[,mvpa_daily_total := mvpa1+mvpa2+mvpa3+mvpa4+mvpa5+mvpa6+mvpa7] 

##################### Step 4: Categorize into WW
## Define categories (any two combination of days >= 50% of total MVPA)
complete_time[,':='(ww_pattern = ifelse(c(((mvpa1+mvpa2) >= (mvpa_daily_total/2)) | ((mvpa1+mvpa3) >= (mvpa_daily_total/2)) | ((mvpa1+mvpa4) >= (mvpa_daily_total/2)) |
                                            ((mvpa1+mvpa5) >= (mvpa_daily_total/2)) | ((mvpa1+mvpa6) >= (mvpa_daily_total/2)) |
                                            ((mvpa1+mvpa7) >= (mvpa_daily_total/2)) | ((mvpa2+mvpa3) >= (mvpa_daily_total/2)) |
                                            ((mvpa2+mvpa4) >= (mvpa_daily_total/2)) | ((mvpa2+mvpa5) >= (mvpa_daily_total/2)) |
                                            ((mvpa2+mvpa6) >= (mvpa_daily_total/2)) | ((mvpa2+mvpa7) >= (mvpa_daily_total/2)) | 
                                            ((mvpa3+mvpa4) >= (mvpa_daily_total/2)) | ((mvpa3+mvpa5) >= (mvpa_daily_total/2)) |
                                            ((mvpa3+mvpa6) >= (mvpa_daily_total/2)) | ((mvpa3+mvpa7) >= (mvpa_daily_total/2)) |
                                            ((mvpa4+mvpa5) >= (mvpa_daily_total/2)) | ((mvpa4+mvpa6) >= (mvpa_daily_total/2)) |
                                            ((mvpa4+mvpa7) >= (mvpa_daily_total/2)) | ((mvpa5+mvpa6) >= (mvpa_daily_total/2)) |
                                            ((mvpa5+mvpa7) >= (mvpa_daily_total/2)) | ((mvpa6+mvpa7) >= (mvpa_daily_total/2))),1,0))]

complete_time[,':='(who_acc = ifelse(mvpa_daily_total >= 150,1,0))]
complete_time[,':='(who_acc_ww = ifelse(c(who_acc==1 & ww_pattern==1),1,0),
                    who_acc_nww = ifelse(c(who_acc==1 & ww_pattern==0),1,0))]
complete_time[,':='(activity_group = ifelse(who_acc==0,'Inactive',
                                            ifelse(who_acc_ww==1,'Active - WW','Active - Regular')))]

# Sample-based thresholds
## Median (230.4)
complete_time[,':='(above_median = ifelse(mvpa_daily_total >= quantile(mvpa_daily_total,0.5),1,0))]
complete_time[,':='(median_ww = ifelse(c(above_median==1 & ww_pattern==1),1,0),
                    median_nww = ifelse(c(above_median==1 & ww_pattern==0),1,0))]
complete_time[,':='(activity_group_median = ifelse(above_median==0,'Inactive',
                                                   ifelse(median_ww==1,'Active - WW','Active - Regular')))]

## 25th percentile (115.2)
complete_time[,':='(above_quarter = ifelse(mvpa_daily_total >= quantile(mvpa_daily_total,0.25),1,0))]
complete_time[,':='(quarter_ww = ifelse(c(above_quarter==1 & ww_pattern==1),1,0),
                    quarter_nww = ifelse(c(above_quarter==1 & ww_pattern==0),1,0))]
complete_time[,':='(activity_group_quarter = ifelse(above_quarter==0,'Inactive',
                                                    ifelse(quarter_ww==1,'Active - WW','Active - Regular')))]

## 75th percentile (403.2)
complete_time[,':='(above_sf = ifelse(mvpa_daily_total >= quantile(mvpa_daily_total,0.75),1,0))]
complete_time[,':='(sf_ww = ifelse(c(above_sf==1 & ww_pattern==1),1,0),
                    sf_nww = ifelse(c(above_sf==1 & ww_pattern==0),1,0))]
complete_time[,':='(activity_group_sf = ifelse(above_sf==0,'Inactive',
                                               ifelse(sf_ww==1,'Active - WW','Active - Regular')))]

# MVPA deciles
complete_time[,":="(mvpa_decile = classifier(mvpa_daily_total,ncuts=10))]

##################### Step 5: Prepare outcome variables
# Loads
af <- fread(file='/phenotypes/af_202109.csv')
mi <- fread(file='/phenotypes/mi_icd_202109.tsv')
hf <- fread(file='/phenotypes/hf_inclusive_202109.tsv')
stroke <- fread(file='/phenotypes/stroke_icd_202109.tsv')
msk <- fread(file='/phenotypes/msk_202109.tsv')
censor <- fread(file='/phenotypes/censor_202109.csv')

setkey(censor,sample_id); setkey(complete_time,sample_id);

## Death handling
complete_time[censor, ':='(death_date = i.death_date, death_censor_date = i.death_censor_date)]

complete_time[,died := ifelse(!is.na(death_date),1,0)]
complete_time[,time_to_death := (ifelse(!is.na(death_date),(as.numeric(death_date) - as.numeric(end_date)),
                                        (as.numeric(death_censor_date) - as.numeric(end_date)))/365.25)]
complete_time[,dummy := 1]

setDF(complete_time)
death_ci <- survivor(data=complete_time,risk_data='dummy',event='died',time='time_to_death',breakpoint=5)
setDT(complete_time)

# Remove small number of people with negative follow up time (4)
complete_time <- complete_time[!(time_to_death < 0)] # 89577 - 4 = 89573

# Censor at death or last follow-up for non-events
af[,censor_date := as.Date(ifelse(has_disease==0 & has_died == 1,
                                  pmin(death_censor_date,censor_date),censor_date),origin='1970-01-01')]
mi[,censor_date := as.Date(ifelse(has_disease==0 & has_died == 1,
                                  pmin(death_censor_date,censor_date),censor_date),origin='1970-01-01')]
hf[,censor_date := as.Date(ifelse(has_disease==0 & has_died == 1,
                                  pmin(death_censor_date,censor_date),censor_date),origin='1970-01-01')]
stroke[,censor_date := as.Date(ifelse(has_disease==0 & has_died == 1,
                                      pmin(death_censor_date,censor_date),censor_date),origin='1970-01-01')]
msk[,censor_date := as.Date(ifelse(has_disease==0 & has_died == 1,
                                   pmin(death_censor_date,censor_date),censor_date),origin='1970-01-01')]

# Generate phenotypes
## Clean up and prejoin
setkey(af,sample_id); setkey(mi,sample_id); setkey(hf,sample_id)
setkey(stroke,sample_id); setkey(msk,sample_id); setkey(complete_time,sample_id);

## Joins
complete_time[af,":="(has_af = i.has_disease, prevalent_af = i.prevalent_disease, 
                      incident_af = i.incident_disease, af_censor_date = i.censor_date)]
complete_time[mi,":="(has_mi = i.has_disease, prevalent_mi = i.prevalent_disease, 
                      incident_mi = i.incident_disease, mi_censor_date = i.censor_date)]
complete_time[hf,":="(has_hf = i.has_disease, prevalent_hf = i.prevalent_disease, 
                      incident_hf = i.incident_disease, hf_censor_date = i.censor_date)]
complete_time[stroke,":="(has_stroke = i.has_disease, prevalent_stroke = i.prevalent_disease, 
                          incident_stroke = i.incident_disease, stroke_censor_date = i.censor_date)]
complete_time[msk,":="(has_msk = i.has_disease, prevalent_msk = i.prevalent_disease, 
                       incident_msk = i.incident_disease, msk_censor_date = i.censor_date)]

# Site-specific censoring dates
# Format dates
dates <- c('af_censor_date','mi_censor_date','hf_censor_date',
           'stroke_censor_date','msk_censor_date')
for (j in dates){set(complete_time,j=j,value=as.Date(complete_time[[j]],format='%Y-%m-%d'))}

# Load center categories
center <- fread(file='/Volumes/medpop_afib/skhurshid/phenotypes/center0.csv')
center_lookup <- fread(file='/Volumes/medpop_afib/skhurshid/phenotypes/enrollment_correspondences.csv')

# Add center value to dataset
setkey(complete_time,sample_id); setkey(center,sample_id)
complete_time[center,':='(center_code = i.value)]

setkey(complete_time,center_code); setkey(center_lookup,Code)
complete_time[center_lookup,':='(center_location = i.Region)]

# Modify censor dates based on location
complete_time[,':='(af_censor_date = as.Date(ifelse(c(has_af==1 | center_location=='England'),af_censor_date,
                                                    ifelse(center_location=='Scotland',pmin(af_censor_date,as.Date('2021-03-31',format='%Y-%m-%d')),
                                                           pmin(af_censor_date,as.Date('2018-02-28',format='%Y-%m-%d')))),origin='1970-01-01'))]
complete_time[,':='(mi_censor_date = as.Date(ifelse(c(has_mi==1 | center_location=='England'),mi_censor_date,
                                                    ifelse(center_location=='Scotland',pmin(mi_censor_date,as.Date('2021-03-31',format='%Y-%m-%d')),
                                                           pmin(mi_censor_date,as.Date('2018-02-28',format='%Y-%m-%d')))),origin='1970-01-01'))]
complete_time[,':='(hf_censor_date = as.Date(ifelse(c(has_hf==1 | center_location=='England'),hf_censor_date,
                                                    ifelse(center_location=='Scotland',pmin(hf_censor_date,as.Date('2021-03-31',format='%Y-%m-%d')),
                                                           pmin(hf_censor_date,as.Date('2018-02-28',format='%Y-%m-%d')))),origin='1970-01-01'))]
complete_time[,':='(stroke_censor_date = as.Date(ifelse(c(has_stroke==1 | center_location=='England'),stroke_censor_date,
                                                        ifelse(center_location=='Scotland',pmin(stroke_censor_date,as.Date('2021-03-31',format='%Y-%m-%d')),
                                                               pmin(stroke_censor_date,as.Date('2018-02-28',format='%Y-%m-%d')))),origin='1970-01-01'))]
complete_time[,':='(msk_censor_date = as.Date(ifelse(c(has_msk==1 | center_location=='England'),msk_censor_date,
                                                     ifelse(center_location=='Scotland',pmin(msk_censor_date,as.Date('2021-03-31',format='%Y-%m-%d')),
                                                            pmin(msk_censor_date,as.Date('2018-02-28',format='%Y-%m-%d')))),origin='1970-01-01'))]

# Define prev and incd with respect to accel
# Prev
complete_time[,":="(prev_af_accel = ifelse(c(!is.na(has_af) & (has_af==1) 
                                             & (af_censor_date <= end_date)),1,0),
                    prev_mi_accel = ifelse(c(!is.na(has_mi) & (has_mi==1) 
                                             & (mi_censor_date <= end_date)),1,0),
                    prev_hf_accel = ifelse(c(!is.na(has_hf) & (has_hf==1) 
                                             & (hf_censor_date <= end_date)),1,0),
                    prev_stroke_accel = ifelse(c(!is.na(has_stroke) & (has_stroke==1) 
                                                 & (stroke_censor_date <= end_date)),1,0),
                    prev_msk_accel = ifelse(c(!is.na(has_msk) & (has_msk==1) 
                                              & (msk_censor_date <= end_date)),1,0))]
# Incd
complete_time[,":="(incd_af_accel = ifelse(c(!is.na(has_af) & (has_af==1) 
                                             & (af_censor_date > end_date)),1,0),
                    incd_mi_accel = ifelse(c(!is.na(has_mi) & (has_mi==1) 
                                             & (mi_censor_date > end_date)),1,0),
                    incd_hf_accel = ifelse(c(!is.na(has_hf) & (has_hf==1) 
                                             & (hf_censor_date > end_date)),1,0),
                    incd_stroke_accel = ifelse(c(!is.na(has_stroke) & (has_stroke==1) 
                                                 & (stroke_censor_date > end_date)),1,0),
                    incd_msk_accel = ifelse(c(!is.na(has_msk) & (has_msk==1) 
                                              & (msk_censor_date > end_date)),1,0))]

# Time
complete_time[,":="(time_accel_to_af = as.numeric(af_censor_date - end_date)/365.25,
                    time_accel_to_mi = as.numeric(mi_censor_date - end_date)/365.25,
                    time_accel_to_hf = as.numeric(hf_censor_date - end_date)/365.25,
                    time_accel_to_stroke = as.numeric(stroke_censor_date - end_date)/365.25,
                    time_accel_to_msk = as.numeric(msk_censor_date - end_date)/365.25)]

# Age
setkey(complete_time,sample_id); setkey(censor,sample_id)
complete_time[censor,birth_date := i.birthdate]
complete_time[,age_accel := (as.numeric(end_date) - as.numeric(birth_date))/365.25]

# Sex
complete_time[censor,sex := i.sex]

# Save out intermediate
#write.csv(complete_time,file='/weekend_warrior/complete_time_052223.csv',row.names=F)

##################### Step 6: Outcomes (Main analysis)
## Relevels
complete_time[,activity_group := factor(activity_group,levels=c('Inactive','Active - WW','Active - Regular'))]
complete_time[,activity_group_active_ref := factor(activity_group,levels=c('Active - WW','Active - Regular','Inactive'))]
complete_time[,activity_group_median := factor(activity_group_median,levels=c('Inactive','Active - WW','Active - Regular'))]

##################### AF
af_set <- complete_time[c((prev_af_accel == 0) & (time_accel_to_af > 0))]

### Models
mod_af <- coxph(Surv(time_accel_to_af,incd_af_accel) ~ activity_group + age_accel + sex + race_category_adjust +
                  tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=af_set)
mod_af_far <- float(mod_af)

mod_af_median <- coxph(Surv(time_accel_to_af,incd_af_accel) ~ activity_group_median + age_accel + sex + race_category_adjust +
                         tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=af_set)
mod_af_far_median <- float(mod_af_median)

### KM
# Graphical variables
af_set[,activity_group_graphical := factor(activity_group,levels=c('Inactive','Active - Regular','Active - WW'))]

# AF KM
prodlim_af <- prodlim(Hist(time_accel_to_af,incd_af_accel)~activity_group_graphical,data=af_set)

CairoPDF(file='/weekend_warrior/km_all_af.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
plot(prodlim_af,"cuminc",ylim=c(0,0.04),xlim=c(0,5),
     axis2.at=seq(0,0.04,0.01),axis2.las=2,lwd=1.4,background=F,
     axis1.at=c(0,1,2,3,4,5),axis1.labels=as.character(0:5),
     atrisk.times=c(0,1,2,3,4,5),col=c("#d95f02",'darkgray','#1b9e77'),atrisk.col='black',confint=TRUE,
     legend.x=0,legend.y=0.04,axis1.cex.axis=2.5,axis2.cex.axis=2.5,axis1.padj=0.5,
     legend.cex=2.2,legend.legend=c("Inactive","Active-Regular","Active-Weekend Warrior"),
     atrisk.title=("                     "),atrisk.pos=-0.5,atrisk.line=c(1.2,2.8,4.4),
     atrisk.cex=1.8,atrisk.interspace=1.4,xlab='',ylab='')
mtext("Cumulative risk (%)",side=2,line=-1.2,at=0.02,cex=2.5)
mtext("Years",side=1, line=-1,cex=2.5)
mtext('Stratum',side=1, line=-0.3,cex=1.8,at=-1.0)
dev.off()

### KM MEDIAN
# Graphical variables
af_set[,activity_group_median_graphical := factor(activity_group_median,levels=c('Inactive','Active - Regular','Active - WW'))]

# AF KM
prodlim_af <- prodlim(Hist(time_accel_to_af,incd_af_accel)~activity_group_median_graphical,data=af_set)

CairoPDF(file='/weekend_warrior/km_all_af_median.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
plot(prodlim_af,"cuminc",ylim=c(0,0.04),xlim=c(0,5),
     axis2.at=seq(0,0.04,0.01),axis2.las=2,lwd=1.4,background=F,
     axis1.at=c(0,1,2,3,4,5),axis1.labels=as.character(0:5),
     atrisk.times=c(0,1,2,3,4,5),col=c("#d95f02",'darkgray','#1b9e77'),atrisk.col='black',confint=TRUE,
     legend.x=0,legend.y=0.04,axis1.cex.axis=2.5,axis2.cex.axis=2.5,axis1.padj=0.5,
     legend.cex=2.2,legend.legend=c("Inactive","Active-Regular","Active-Weekend Warrior"),
     atrisk.title=("                     "),atrisk.pos=-0.5,atrisk.line=c(1.2,2.8,4.4),
     atrisk.cex=1.8,atrisk.interspace=1.8,xlab='',ylab='')
mtext("Cumulative risk (%)",side=2,line=-1.2,at=0.02,cex=2.5)
mtext("Years",side=1, line=-1,cex=2.5)
mtext('Stratum',side=1, line=-0.3,cex=1.8,at=-1.0)
dev.off()

##################### MI
mi_set <- complete_time[c((prev_mi_accel == 0) & (time_accel_to_mi > 0))]

### Models
mod_mi <- coxph(Surv(time_accel_to_mi,incd_mi_accel) ~ activity_group + age_accel + sex + race_category_adjust +
                  tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=mi_set)
mod_mi_far <- float(mod_mi)

mod_mi_median <- coxph(Surv(time_accel_to_mi,incd_mi_accel) ~ activity_group_median + age_accel + sex + race_category_adjust +
                         tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=mi_set)
mod_mi_far_median <- float(mod_mi_median)

# Graphical variables
mi_set[,activity_group_graphical := factor(activity_group,levels=c('Inactive','Active - Regular','Active - WW'))]

### KM
prodlim_mi <- prodlim(Hist(time_accel_to_mi,incd_mi_accel)~activity_group_graphical,data=mi_set)

CairoPDF(file='/weekend_warrior/km_all_mi.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
plot(prodlim_mi,"cuminc",ylim=c(0,0.04),xlim=c(0,5),
     axis2.at=seq(0,0.04,0.01),axis2.las=2,lwd=1.4,background=F,
     axis1.at=c(0,1,2,3,4,5),axis1.labels=as.character(0:5),
     atrisk.times=c(0,1,2,3,4,5),col=c("#d95f02",'darkgray','#1b9e77'),atrisk.col='black',confint=TRUE,
     legend.x=0,legend.y=0.04,axis1.cex.axis=2.5,axis2.cex.axis=2.5,axis1.padj=0.5,
     legend.cex=2.2,legend.legend=c("Inactive","Active-Regular","Active-Weekend Warrior"),
     atrisk.title=("                     "),atrisk.pos=-0.5,atrisk.line=c(1.2,2.8,4.4),
     atrisk.cex=1.8,atrisk.interspace=1.4,xlab='',ylab='')
mtext("Cumulative risk (%)",side=2,line=-1.2,at=0.02,cex=2.5)
mtext("Years",side=1, line=-1,cex=2.5)
mtext('Stratum',side=1, line=-0.3,cex=1.8,at=-1.0)
dev.off()

# Graphical variables
mi_set[,activity_group_graphical_median := factor(activity_group_median,levels=c('Inactive','Active - Regular','Active - WW'))]

### KM
prodlim_mi <- prodlim(Hist(time_accel_to_mi,incd_mi_accel)~activity_group_graphical_median,data=mi_set)

CairoPDF(file='/weekend_warrior/km_all_mi_median.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
plot(prodlim_mi,"cuminc",ylim=c(0,0.04),xlim=c(0,5),
     axis2.at=seq(0,0.04,0.01),axis2.las=2,lwd=1.4,background=F,
     axis1.at=c(0,1,2,3,4,5),axis1.labels=as.character(0:5),
     atrisk.times=c(0,1,2,3,4,5),col=c("#d95f02",'darkgray','#1b9e77'),atrisk.col='black',confint=TRUE,
     legend.x=0,legend.y=0.04,axis1.cex.axis=2.5,axis2.cex.axis=2.5,axis1.padj=0.5,
     legend.cex=2.2,legend.legend=c("Inactive","Active-Regular","Active-Weekend Warrior"),
     atrisk.title=("                     "),atrisk.pos=-0.5,atrisk.line=c(1.2,2.8,4.4),
     atrisk.cex=1.8,atrisk.interspace=1.4,xlab='',ylab='')
mtext("Cumulative risk (%)",side=2,line=-1.2,at=0.02,cex=2.5)
mtext("Years",side=1, line=-1,cex=2.5)
mtext('Stratum',side=1, line=-0.3,cex=1.8,at=-1.0)
dev.off()

##################### HF
hf_set <- complete_time[c((prev_hf_accel == 0) & (time_accel_to_hf > 0))]

### Models
mod_hf <- coxph(Surv(time_accel_to_hf,incd_hf_accel) ~ activity_group + age_accel + sex + race_category_adjust +
                  tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=hf_set)
mod_hf_far <- float(mod_hf)

mod_hf_median <- coxph(Surv(time_accel_to_hf,incd_hf_accel) ~ activity_group_median + age_accel + sex + race_category_adjust +
                         tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=hf_set)
mod_hf_far_median <- float(mod_hf_median)

# Graphical variables
hf_set[,activity_group_graphical := factor(activity_group,levels=c('Inactive','Active - Regular','Active - WW'))]

### KM
prodlim_hf <- prodlim(Hist(time_accel_to_hf,incd_hf_accel)~activity_group_graphical,data=hf_set)

CairoPDF(file='/weekend_warrior/km_all_hf.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
plot(prodlim_hf,"cuminc",ylim=c(0,0.04),xlim=c(0,5),
     axis2.at=seq(0,0.04,0.01),axis2.las=2,lwd=1.4,background=F,
     axis1.at=c(0,1,2,3,4,5),axis1.labels=as.character(0:5),
     atrisk.times=c(0,1,2,3,4,5),col=c("#d95f02",'darkgray','#1b9e77'),atrisk.col='black',confint=TRUE,
     legend.x=0,legend.y=0.04,axis1.cex.axis=2.5,axis2.cex.axis=2.5,axis1.padj=0.5,
     legend.cex=2.2,legend.legend=c("Inactive","Active-Regular","Active-Weekend Warrior"),
     atrisk.title=("                     "),atrisk.pos=-0.5,atrisk.line=c(1.2,2.8,4.4),
     atrisk.cex=1.8,atrisk.interspace=1.4,xlab='',ylab='')
mtext("Cumulative risk (%)",side=2,line=-1.2,at=0.02,cex=2.5)
mtext("Years",side=1, line=-1,cex=2.5)
mtext('Stratum',side=1, line=-0.3,cex=1.8,at=-1.0)
dev.off()

# Graphical variables
hf_set[,activity_group_graphical_median := factor(activity_group_median,levels=c('Inactive','Active - Regular','Active - WW'))]

### KM
prodlim_hf <- prodlim(Hist(time_accel_to_hf,incd_hf_accel)~activity_group_graphical_median,data=mi_set)

CairoPDF(file='/weekend_warrior/km_all_hf_median.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
plot(prodlim_hf,"cuminc",ylim=c(0,0.04),xlim=c(0,5),
     axis2.at=seq(0,0.04,0.01),axis2.las=2,lwd=1.4,background=F,
     axis1.at=c(0,1,2,3,4,5),axis1.labels=as.character(0:5),
     atrisk.times=c(0,1,2,3,4,5),col=c("#d95f02",'darkgray','#1b9e77'),atrisk.col='black',confint=TRUE,
     legend.x=0,legend.y=0.04,axis1.cex.axis=2.5,axis2.cex.axis=2.5,axis1.padj=0.5,
     legend.cex=2.2,legend.legend=c("Inactive","Active-Regular","Active-Weekend Warrior"),
     atrisk.title=("                     "),atrisk.pos=-0.5,atrisk.line=c(1.2,2.8,4.4),
     atrisk.cex=1.8,atrisk.interspace=1.4,xlab='',ylab='')
mtext("Cumulative risk (%)",side=2,line=-1.2,at=0.02,cex=2.5)
mtext("Years",side=1, line=-1,cex=2.5)
mtext('Stratum',side=1, line=-0.3,cex=1.8,at=-1.0)
dev.off()

##################### Stroke
stroke_set <- complete_time[c((prev_stroke_accel == 0) & (time_accel_to_stroke > 0))]

### Models
mod_stroke <- coxph(Surv(time_accel_to_stroke,incd_stroke_accel) ~ activity_group + age_accel + sex + race_category_adjust +
                      tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=stroke_set)
mod_stroke_far <- float(mod_stroke)

mod_stroke_median <- coxph(Surv(time_accel_to_stroke,incd_stroke_accel) ~ activity_group_median + age_accel + sex + race_category_adjust +
                             tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=stroke_set)
mod_stroke_far_median <- float(mod_stroke_median)

stroke_set[,activity_group_graphical := factor(activity_group,levels=c('Inactive','Active - Regular','Active - WW'))]

### KM
prodlim_stroke <- prodlim(Hist(time_accel_to_stroke,incd_stroke_accel)~activity_group_graphical,data=stroke_set)

CairoPDF(file='/weekend_warrior/km_all_stroke.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
plot(prodlim_stroke,"cuminc",ylim=c(0,0.04),xlim=c(0,5),
     axis2.at=seq(0,0.04,0.01),axis2.las=2,lwd=1.4,background=F,
     axis1.at=c(0,1,2,3,4,5),axis1.labels=as.character(0:5),
     atrisk.times=c(0,1,2,3,4,5),col=c("#d95f02",'darkgray','#1b9e77'),atrisk.col='black',confint=TRUE,
     legend.x=0,legend.y=0.04,axis1.cex.axis=2.5,axis2.cex.axis=2.5,axis1.padj=0.5,
     legend.cex=2.2,legend.legend=c("Inactive","Active-Regular","Active-Weekend Warrior"),
     atrisk.title=("                     "),atrisk.pos=-0.5,atrisk.line=c(1.2,2.8,4.4),
     atrisk.cex=1.8,atrisk.interspace=1.4,xlab='',ylab='')
mtext("Cumulative risk (%)",side=2,line=-1.2,at=0.02,cex=2.5)
mtext("Years",side=1, line=-1,cex=2.5)
mtext('Stratum',side=1, line=-0.3,cex=1.8,at=-1.0)
dev.off()

stroke_set[,activity_group_graphical_median := factor(activity_group_median,levels=c('Inactive','Active - Regular','Active - WW'))]

### KM
prodlim_stroke <- prodlim(Hist(time_accel_to_stroke,incd_stroke_accel)~activity_group_graphical_median,data=stroke_set)

CairoPDF(file='/weekend_warrior/km_all_stroke_median.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
plot(prodlim_stroke,"cuminc",ylim=c(0,0.04),xlim=c(0,5),
     axis2.at=seq(0,0.04,0.01),axis2.las=2,lwd=1.4,background=F,
     axis1.at=c(0,1,2,3,4,5),axis1.labels=as.character(0:5),
     atrisk.times=c(0,1,2,3,4,5),col=c("#d95f02",'darkgray','#1b9e77'),atrisk.col='black',confint=TRUE,
     legend.x=0,legend.y=0.04,axis1.cex.axis=2.5,axis2.cex.axis=2.5,axis1.padj=0.5,
     legend.cex=2.2,legend.legend=c("Inactive","Active-Regular","Active-Weekend Warrior"),
     atrisk.title=("                     "),atrisk.pos=-0.5,atrisk.line=c(1.2,2.8,4.4),
     atrisk.cex=1.8,atrisk.interspace=1.4,xlab='',ylab='')
mtext("Cumulative risk (%)",side=2,line=-1.2,at=0.02,cex=2.5)
mtext("Years",side=1, line=-1,cex=2.5)
mtext('Stratum',side=1, line=-0.3,cex=1.8,at=-1.0)
dev.off()

##################### Step 7: Stratified cox model
# Classify deciles within each subgroup
af_set[,":="(mvpa_decile = factor(classifier(mvpa_daily_total,ncuts=10)))]
mi_set[,":="(mvpa_decile = factor(classifier(mvpa_daily_total,ncuts=10)))]
hf_set[,":="(mvpa_decile = factor(classifier(mvpa_daily_total,ncuts=10)))]
stroke_set[,":="(mvpa_decile = factor(classifier(mvpa_daily_total,ncuts=10)))]

### Overall models
mod_af <- coxph(Surv(time_accel_to_af,incd_af_accel) ~ ww_pattern + age_accel + sex + race_category_adjust +
                  tob + tdi + etoh_grams + diet + qual_ea + employment_status + self_health + factor(strata(mvpa_decile)),data=af_set)
mod_mi <- coxph(Surv(time_accel_to_mi,incd_mi_accel) ~ ww_pattern + age_accel + sex + race_category_adjust +
                  tob + tdi + etoh_grams + diet + qual_ea + employment_status + self_health + factor(strata(mvpa_decile)),data=mi_set)
mod_hf <- coxph(Surv(time_accel_to_hf,incd_hf_accel) ~ ww_pattern + age_accel + sex + race_category_adjust +
                  tob + tdi + etoh_grams + diet + qual_ea + employment_status + self_health + factor(strata(mvpa_decile)),data=hf_set)
mod_stroke <- coxph(Surv(time_accel_to_stroke,incd_stroke_accel) ~ ww_pattern + age_accel + sex + race_category_adjust +
                      tob + tdi + etoh_grams + diet + qual_ea + employment_status + self_health + factor(strata(mvpa_decile)),data=stroke_set)

##################### Step 8: 75% sensi analysis (WW defined as >= 75% of total MVPA over 1-2 days)
complete_time[,':='(ww_pattern75 = ifelse(c((mvpa1+mvpa2 >= mvpa_daily_total*0.75) |
                                              (mvpa1+mvpa3 >= mvpa_daily_total*0.75) | (mvpa1+mvpa4 >= mvpa_daily_total*0.75) |
                                              (mvpa1+mvpa5 >= mvpa_daily_total*0.75) | (mvpa1+mvpa6 >= mvpa_daily_total*0.75) |
                                              (mvpa1+mvpa7 >= mvpa_daily_total*0.75) | (mvpa2+mvpa3 >= mvpa_daily_total*0.75) |
                                              (mvpa2+mvpa4 >= mvpa_daily_total*0.75) | (mvpa2+mvpa5 >= mvpa_daily_total*0.75) |
                                              (mvpa2+mvpa6 >= mvpa_daily_total*0.75) | (mvpa2+mvpa7 >= mvpa_daily_total*0.75) | 
                                              (mvpa3+mvpa4 >= mvpa_daily_total*0.75) | (mvpa3+mvpa5 >= mvpa_daily_total*0.75) |
                                              (mvpa3+mvpa6 >= mvpa_daily_total*0.75) | (mvpa3+mvpa7 >= mvpa_daily_total*0.75) |
                                              (mvpa4+mvpa5 >= mvpa_daily_total*0.75) | (mvpa4+mvpa6 >= mvpa_daily_total*0.75) |
                                              (mvpa4+mvpa7 >= mvpa_daily_total*0.75) | (mvpa5+mvpa6 >= mvpa_daily_total*0.75) |
                                              (mvpa5+mvpa7 >= mvpa_daily_total*0.75) | (mvpa6+mvpa7 >= mvpa_daily_total*0.75)),1,0))]

complete_time[,':='(who_acc_ww75 = ifelse(c(who_acc==1 & ww_pattern75==1),1,0),
                    who_acc_nww75 = ifelse(c(who_acc==1 & ww_pattern75==0),1,0))]
complete_time[,':='(activity_group75 = ifelse(who_acc==0,'Inactive',
                                              ifelse(who_acc_ww75==1,'Active - WW','Active - Regular')))]

complete_time[,':='(median_ww75 = ifelse(c(above_median==1 & ww_pattern75==1),1,0),
                    median_nww75 = ifelse(c(above_median==1 & ww_pattern75==0),1,0))]
complete_time[,':='(activity_group75_median = ifelse(above_median==0,'Inactive',
                                                     ifelse(median_ww75==1,'Active - WW','Active - Regular')))]

##################### Outcomes
## Relevels
complete_time[,activity_group75 := factor(activity_group75,levels=c('Inactive','Active - WW','Active - Regular'))]
complete_time[,activity_group75_median := factor(activity_group75_median,levels=c('Inactive','Active - WW','Active - Regular'))]

##################### AF
af_set <- complete_time[c((prev_af_accel == 0) & (time_accel_to_af > 0))]

### Models
mod_af <- coxph(Surv(time_accel_to_af,incd_af_accel) ~ activity_group75 + age_accel + sex + race_category_adjust +
                  tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=af_set)
mod_af_far <- float(mod_af)

mod_af_median <- coxph(Surv(time_accel_to_af,incd_af_accel) ~ activity_group75_median + age_accel + sex + race_category_adjust +
                         tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=af_set)
mod_af_far_median <- float(mod_af_median)

##################### MI
mi_set <- complete_time[c((prev_mi_accel == 0) & (time_accel_to_mi > 0))]

### Models
mod_mi <- coxph(Surv(time_accel_to_mi,incd_mi_accel) ~ activity_group75 + age_accel + sex + race_category_adjust +
                  tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=mi_set)
mod_mi_far <- float(mod_mi)

mod_mi_median <- coxph(Surv(time_accel_to_mi,incd_mi_accel) ~ activity_group75_median + age_accel + sex + race_category_adjust +
                         tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=mi_set)
mod_mi_far_median <- float(mod_mi_median)

##################### HF
hf_set <- complete_time[c((prev_hf_accel == 0) & (time_accel_to_hf > 0))]

### Models
mod_hf <- coxph(Surv(time_accel_to_hf,incd_hf_accel) ~ activity_group75 + age_accel + sex + race_category_adjust +
                  tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=hf_set)
mod_hf_far <- float(mod_hf)

mod_hf_median <- coxph(Surv(time_accel_to_hf,incd_hf_accel) ~ activity_group75_median + age_accel + sex + race_category_adjust +
                         tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=hf_set)
mod_hf_far_median <- float(mod_hf_median)

##################### Stroke
stroke_set <- complete_time[c((prev_stroke_accel == 0) & (time_accel_to_stroke > 0))]

### Models
mod_stroke <- coxph(Surv(time_accel_to_stroke,incd_stroke_accel) ~ activity_group75 + age_accel + sex + race_category_adjust +
                      tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=stroke_set)
mod_stroke_far <- float(mod_stroke)

mod_stroke_median <- coxph(Surv(time_accel_to_stroke,incd_stroke_accel) ~ activity_group75_median + age_accel + sex + race_category_adjust +
                             tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=stroke_set)
mod_stroke_far_median <- float(mod_stroke_median)

##################### Step 9: True weekends (defined as days 6 and 7)
## Define categories
complete_time[,':='(weekend_ww = ifelse(c((mvpa6+mvpa7) >= mvpa_daily_total/2),1,0))]
complete_time[,':='(who_acc_weekend_ww = ifelse(c(who_acc==1 & weekend_ww==1),1,0),
                    who_acc_weekend_nww = ifelse(c(who_acc==1 & weekend_ww==0),1,0))]
complete_time[,':='(weekend_activity_group = ifelse(who_acc==0,'Inactive',
                                                    ifelse(who_acc_weekend_ww==1,'Active - WW','Active - Regular')))]

complete_time[,':='(median_weekend_ww = ifelse(c(above_median==1 & weekend_ww==1),1,0),
                    median_weekend_nww = ifelse(c(above_median==1 & weekend_ww==0),1,0))]
complete_time[,':='(weekend_activity_group_median = ifelse(above_median==0,'Inactive',
                                                           ifelse(median_weekend_ww==1,'Active - WW','Active - Regular')))]

##################### Outcomes
## Relevels
complete_time[,weekend_activity_group := factor(weekend_activity_group,levels=c('Inactive','Active - WW','Active - Regular'))]
complete_time[,weekend_activity_group_median := factor(weekend_activity_group_median,levels=c('Inactive','Active - WW','Active - Regular'))]

##################### AF
af_set <- complete_time[c((prev_af_accel == 0) & (time_accel_to_af > 0))]

### Models
mod_af <- coxph(Surv(time_accel_to_af,incd_af_accel) ~ weekend_activity_group + age_accel + sex + race_category_adjust +
                  tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=af_set)
mod_af_far <- float(mod_af)

mod_af_median <- coxph(Surv(time_accel_to_af,incd_af_accel) ~ weekend_activity_group_median + age_accel + sex + race_category_adjust +
                         tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=af_set)
mod_af_far_median <- float(mod_af_median)

##################### MI
mi_set <- complete_time[c((prev_mi_accel == 0) & (time_accel_to_mi > 0))]

### Models
mod_mi <- coxph(Surv(time_accel_to_mi,incd_mi_accel) ~ weekend_activity_group + age_accel + sex + race_category_adjust +
                  tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=mi_set)
mod_mi_far <- float(mod_mi)

mod_mi_median <- coxph(Surv(time_accel_to_mi,incd_mi_accel) ~ weekend_activity_group_median + age_accel + sex + race_category_adjust +
                         tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=mi_set)
mod_mi_far_median <- float(mod_mi_median)

##################### HF
hf_set <- complete_time[c((prev_hf_accel == 0) & (time_accel_to_hf > 0))]

### Models
mod_hf <- coxph(Surv(time_accel_to_hf,incd_hf_accel) ~ weekend_activity_group + age_accel + sex + race_category_adjust +
                  tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=hf_set)
mod_hf_far <- float(mod_hf)

mod_hf_median <- coxph(Surv(time_accel_to_hf,incd_hf_accel) ~ weekend_activity_group_median + age_accel + sex + race_category_adjust +
                         tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=hf_set)
mod_hf_far_median <- float(mod_hf_median)

##################### Stroke
stroke_set <- complete_time[c((prev_stroke_accel == 0) & (time_accel_to_stroke > 0))]

### Models
mod_stroke <- coxph(Surv(time_accel_to_stroke,incd_stroke_accel) ~ weekend_activity_group + age_accel + sex + race_category_adjust +
                      tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=stroke_set)
mod_stroke_far <- float(mod_stroke)

mod_stroke_median <- coxph(Surv(time_accel_to_stroke,incd_stroke_accel) ~ weekend_activity_group_median + age_accel + sex + race_category_adjust +
                             tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=stroke_set)
mod_stroke_far_median <- float(mod_stroke_median)

##################### Step 10: 2 consecutive days
## Define categories
complete_time[,':='(ww_pattern_consec = ifelse(c((mvpa1+mvpa2 >= mvpa_daily_total/2) | (mvpa2+mvpa3 >= mvpa_daily_total/2) | (mvpa3+mvpa4 >= mvpa_daily_total/2) |
                                                   (mvpa4+mvpa5 >= mvpa_daily_total/2) | (mvpa5+mvpa6 >= mvpa_daily_total/2) |
                                                   (mvpa6+mvpa7 >= mvpa_daily_total/2)),1,0))]

complete_time[,':='(who_acc = ifelse(mvpa_daily_total >= 150,1,0))]
complete_time[,':='(who_acc_consec_ww = ifelse(c(who_acc==1 & ww_pattern_consec==1),1,0),
                    who_acc_consec_nww = ifelse(c(who_acc==1 & ww_pattern_consec==0),1,0))]
complete_time[,':='(activity_group_consec = ifelse(who_acc==0,'Inactive',
                                                   ifelse(who_acc_consec_ww==1,'Active - WW','Active - Regular')))]

complete_time[,':='(median_consec_ww = ifelse(c(above_median==1 & ww_pattern_consec==1),1,0),
                    median_consec_nww = ifelse(c(above_median==1 & ww_pattern_consec==0),1,0))]
complete_time[,':='(activity_group_consec_median = ifelse(above_median==0,'Inactive',
                                                          ifelse(median_consec_ww==1,'Active - WW','Active - Regular')))]

## Relevels
complete_time[,activity_group_consec := factor(activity_group_consec,levels=c('Inactive','Active - WW','Active - Regular'))]
complete_time[,activity_group_consec_median := factor(activity_group_consec_median,levels=c('Inactive','Active - WW','Active - Regular'))]

##################### AF
af_set <- complete_time[c((prev_af_accel == 0) & (time_accel_to_af > 0))]

### Models
mod_af <- coxph(Surv(time_accel_to_af,incd_af_accel) ~ activity_group_consec + age_accel + sex + race_category_adjust +
                  tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=af_set)
mod_af_far <- float(mod_af)

mod_af_median <- coxph(Surv(time_accel_to_af,incd_af_accel) ~ activity_group_consec_median + age_accel + sex + race_category_adjust +
                         tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=af_set)
mod_af_far_median <- float(mod_af_median)

##################### MI
mi_set <- complete_time[c((prev_mi_accel == 0) & (time_accel_to_mi > 0))]

### Models
mod_mi <- coxph(Surv(time_accel_to_mi,incd_mi_accel) ~ activity_group_consec + age_accel + sex + race_category_adjust +
                  tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=mi_set)
mod_mi_far <- float(mod_mi)

mod_mi_median <- coxph(Surv(time_accel_to_mi,incd_mi_accel) ~ activity_group_consec_median + age_accel + sex + race_category_adjust +
                         tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=mi_set)
mod_mi_far_median <- float(mod_mi_median)

##################### HF
hf_set <- complete_time[c((prev_hf_accel == 0) & (time_accel_to_hf > 0))]

### Models
mod_hf <- coxph(Surv(time_accel_to_hf,incd_hf_accel) ~ activity_group_consec + age_accel + sex + race_category_adjust +
                  tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=hf_set)
mod_hf_far <- float(mod_hf)

mod_hf_median <- coxph(Surv(time_accel_to_hf,incd_hf_accel) ~ activity_group_consec_median + age_accel + sex + race_category_adjust +
                         tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=hf_set)
mod_hf_far_median <- float(mod_hf_median)

##################### Stroke
stroke_set <- complete_time[c((prev_stroke_accel == 0) & (time_accel_to_stroke > 0))]

### Models
mod_stroke <- coxph(Surv(time_accel_to_stroke,incd_stroke_accel) ~ activity_group_consec + age_accel + sex + race_category_adjust +
                      tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=stroke_set)
mod_stroke_far <- float(mod_stroke)

mod_stroke_median <- coxph(Surv(time_accel_to_stroke,incd_stroke_accel) ~ activity_group_consec_median + age_accel + sex + race_category_adjust +
                             tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=stroke_set)
mod_stroke_far_median <- float(mod_stroke_median)

##################### Step 11: Additional Thresholds
################## 25
## Relevels
complete_time[,activity_group_quarter := factor(activity_group_quarter,levels=c('Inactive','Active - WW','Active - Regular'))]

##################### AF
af_set <- complete_time[c((prev_af_accel == 0) & (time_accel_to_af > 0))]

### Models
mod_af <- coxph(Surv(time_accel_to_af,incd_af_accel) ~ activity_group_quarter + age_accel + sex + race_category_adjust +
                  tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=af_set)
mod_af_far <- float(mod_af)

##################### MI
mi_set <- complete_time[c((prev_mi_accel == 0) & (time_accel_to_mi > 0))]

### Models
mod_mi <- coxph(Surv(time_accel_to_mi,incd_mi_accel) ~ activity_group_quarter + age_accel + sex + race_category_adjust +
                  tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=mi_set)
mod_mi_far <- float(mod_mi)

##################### HF
hf_set <- complete_time[c((prev_hf_accel == 0) & (time_accel_to_hf > 0))]

### Models
mod_hf <- coxph(Surv(time_accel_to_hf,incd_hf_accel) ~ activity_group_quarter + age_accel + sex + race_category_adjust +
                  tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=hf_set)
mod_hf_far <- float(mod_hf)

##################### Stroke
stroke_set <- complete_time[c((prev_stroke_accel == 0) & (time_accel_to_stroke > 0))]

### Models
mod_stroke <- coxph(Surv(time_accel_to_stroke,incd_stroke_accel) ~ activity_group_quarter + age_accel + sex + race_category_adjust +
                      tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=stroke_set)
mod_stroke_far <- float(mod_stroke)

################## 75
## Relevels
complete_time[,activity_group_sf := factor(activity_group_sf,levels=c('Inactive','Active - WW','Active - Regular'))]

##################### AF
af_set <- complete_time[c((prev_af_accel == 0) & (time_accel_to_af > 0))]

### Models
mod_af <- coxph(Surv(time_accel_to_af,incd_af_accel) ~ activity_group_sf + age_accel + sex + race_category_adjust +
                  tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=af_set)
mod_af_far <- float(mod_af)

##################### MI
mi_set <- complete_time[c((prev_mi_accel == 0) & (time_accel_to_mi > 0))]

### Models
mod_mi <- coxph(Surv(time_accel_to_mi,incd_mi_accel) ~ activity_group_sf + age_accel + sex + race_category_adjust +
                  tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=mi_set)
mod_mi_far <- float(mod_mi)

##################### HF
hf_set <- complete_time[c((prev_hf_accel == 0) & (time_accel_to_hf > 0))]

### Models
mod_hf <- coxph(Surv(time_accel_to_hf,incd_hf_accel) ~ activity_group_sf + age_accel + sex + race_category_adjust +
                  tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=hf_set)
mod_hf_far <- float(mod_hf)

##################### Stroke
stroke_set <- complete_time[c((prev_stroke_accel == 0) & (time_accel_to_stroke > 0))]

### Models
mod_stroke <- coxph(Surv(time_accel_to_stroke,incd_stroke_accel) ~ activity_group_sf + age_accel + sex + race_category_adjust +
                      tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=stroke_set)
mod_stroke_far <- float(mod_stroke)

##################### Step 12: 2-year Blanking period
# Create blanked date variable
complete_time[,":="(blanked_end_date = end_date+365.25*2)]

# Define prev and incd with respect to accel
# Prev
complete_time[,":="(prev_af_accel_blanked = ifelse(c(!is.na(has_af) & (has_af==1) 
                                                     & (af_censor_date <= blanked_end_date)),1,0),
                    prev_mi_accel_blanked = ifelse(c(!is.na(has_mi) & (has_mi==1) 
                                                     & (mi_censor_date <= blanked_end_date)),1,0),
                    prev_hf_accel_blanked = ifelse(c(!is.na(has_hf) & (has_hf==1) 
                                                     & (hf_censor_date <= blanked_end_date)),1,0),
                    prev_stroke_accel_blanked = ifelse(c(!is.na(has_stroke) & (has_stroke==1) 
                                                         & (stroke_censor_date <= blanked_end_date)),1,0))]
# Incd
complete_time[,":="(incd_af_accel_blanked  = ifelse(c(!is.na(has_af) & (has_af==1) 
                                                      & (af_censor_date > blanked_end_date)),1,0),
                    incd_mi_accel_blanked  = ifelse(c(!is.na(has_mi) & (has_mi==1) 
                                                      & (mi_censor_date > blanked_end_date)),1,0),
                    incd_hf_accel_blanked  = ifelse(c(!is.na(has_hf) & (has_hf==1) 
                                                      & (hf_censor_date > blanked_end_date)),1,0),
                    incd_stroke_accel_blanked  = ifelse(c(!is.na(has_stroke) & (has_stroke==1) 
                                                          & (stroke_censor_date > blanked_end_date)),1,0))]

# Time
complete_time[,":="(time_accel_to_af_blanked  = as.numeric(af_censor_date - blanked_end_date)/365.25,
                    time_accel_to_mi_blanked  = as.numeric(mi_censor_date - blanked_end_date)/365.25,
                    time_accel_to_hf_blanked  = as.numeric(hf_censor_date - blanked_end_date)/365.25,
                    time_accel_to_stroke_blanked  = as.numeric(stroke_censor_date - blanked_end_date)/365.25)]

##################### AF
af_set <- complete_time[c((prev_af_accel_blanked == 0) & (time_accel_to_af_blanked > 0))]

mod_af <- coxph(Surv(time_accel_to_af_blanked,incd_af_accel_blanked) ~ activity_group + age_accel + sex + race_category_adjust +
                  tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=af_set)
mod_af_far <- float(mod_af)

mod_af_median <- coxph(Surv(time_accel_to_af_blanked,incd_af_accel_blanked) ~ activity_group_median + age_accel + sex + race_category_adjust +
                         tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=af_set)
mod_af_far_median <- float(mod_af_median)

##################### MI
mi_set <- complete_time[c((prev_mi_accel_blanked == 0) & (time_accel_to_mi_blanked > 0))]

mod_mi <- coxph(Surv(time_accel_to_mi_blanked,incd_mi_accel_blanked) ~ activity_group + age_accel + sex + race_category_adjust +
                  tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=mi_set)
mod_mi_far <- float(mod_mi)

mod_mi_median <- coxph(Surv(time_accel_to_mi_blanked,incd_mi_accel_blanked) ~ activity_group_median + age_accel + sex + race_category_adjust +
                         tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=mi_set)
mod_mi_far_median <- float(mod_mi_median)

##################### HF
hf_set <- complete_time[c((prev_hf_accel_blanked == 0) & (time_accel_to_hf_blanked > 0))]

mod_hf <- coxph(Surv(time_accel_to_hf_blanked,incd_hf_accel_blanked) ~ activity_group + age_accel + sex + race_category_adjust +
                  tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=hf_set)
mod_hf_far <- float(mod_hf)

mod_hf_median <- coxph(Surv(time_accel_to_hf_blanked,incd_hf_accel_blanked) ~ activity_group_median + age_accel + sex + race_category_adjust +
                         tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=hf_set)
mod_hf_far_median <- float(mod_hf_median)

##################### Stroke
stroke_set <- complete_time[c((prev_stroke_accel_blanked == 0) & (time_accel_to_stroke_blanked > 0))]

mod_stroke <- coxph(Surv(time_accel_to_stroke_blanked,incd_stroke_accel_blanked) ~ activity_group + age_accel + sex + race_category_adjust +
                      tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=stroke_set)
mod_stroke_far <- float(mod_stroke)

mod_stroke_median <- coxph(Surv(time_accel_to_stroke_blanked,incd_stroke_accel_blanked) ~ activity_group_median + age_accel + sex + race_category_adjust +
                             tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=stroke_set)
mod_stroke_far_median <- float(mod_stroke_median)

##################### Step 13: Additional adjustment
complete_time_add_covars <- complete_time[c(!is.na(sbp) & !is.na(dbp) & !is.na(bmi) &
                                              !is.na(bpmed))]

# DM
dm <- fread(file='/phenotypes/any_dm_202106.csv')

setkey(dm,sample_id); setkey(complete_time_add_covars,sample_id);
setkey(complete_time,sample_id)

complete_time_add_covars[dm,":="(has_dm = i.has_disease, prevalent_dm = i.prevalent_disease, 
                                 incident_dm = i.incident_disease, dm_censor_date = i.censor_date)]
complete_time[dm,":="(has_dm = i.has_disease, prevalent_dm = i.prevalent_disease, 
                      incident_dm = i.incident_disease, dm_censor_date = i.censor_date)]

# Define prev and incd with respect to accel
# Prev
complete_time_add_covars[,":="(prev_dm_accel = ifelse(c(!is.na(has_dm) & (has_dm==1) 
                                                        & (dm_censor_date <= end_date)),1,0))]

complete_time[,":="(prev_dm_accel = ifelse(c(!is.na(has_dm) & (has_dm==1) 
                                             & (dm_censor_date <= end_date)),1,0))]

##################### AF
af_set <- complete_time_add_covars[c((prev_af_accel == 0) & (time_accel_to_af > 0))]

mod_af <- coxph(Surv(time_accel_to_af,incd_af_accel) ~ activity_group + age_accel + sex + race_category_adjust +
                  tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status + bmi + sbp + dbp + bpmed + prev_dm_accel,data=af_set)
mod_af_far <- float(mod_af)

mod_af_median <- coxph(Surv(time_accel_to_af,incd_af_accel) ~ activity_group_median + age_accel + sex + race_category_adjust +
                         tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status + bmi + sbp + dbp + bpmed + prev_dm_accel,data=af_set)
mod_af_far_median <- float(mod_af_median)

##################### MI
mi_set <- complete_time_add_covars[c((prev_mi_accel == 0) & (time_accel_to_mi > 0))]

mod_mi <- coxph(Surv(time_accel_to_mi,incd_mi_accel) ~ activity_group + age_accel + sex + race_category_adjust +
                  tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status + bmi + sbp + dbp + bpmed + prev_dm_accel,data=mi_set)
mod_mi_far <- float(mod_mi)

mod_mi_median <- coxph(Surv(time_accel_to_mi,incd_mi_accel) ~ activity_group_median + age_accel + sex + race_category_adjust +
                         tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status + bmi + sbp + dbp + bpmed + prev_dm_accel,data=mi_set)
mod_mi_far_median <- float(mod_mi_median)

##################### HF
hf_set <- complete_time_add_covars[c((prev_hf_accel == 0) & (time_accel_to_hf > 0))]

mod_hf <- coxph(Surv(time_accel_to_hf,incd_hf_accel) ~ activity_group + age_accel + sex + race_category_adjust +
                  tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status + bmi + sbp + dbp + bpmed + prev_dm_accel,data=hf_set)
mod_hf_far <- float(mod_hf)

mod_hf_median <- coxph(Surv(time_accel_to_hf,incd_hf_accel) ~ activity_group_median + age_accel + sex + race_category_adjust +
                         tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status + bmi + sbp + dbp + bpmed + prev_dm_accel,data=hf_set)
mod_hf_far_median <- float(mod_hf_median) 

##################### Stroke
stroke_set <- complete_time_add_covars[c((prev_stroke_accel == 0) & (time_accel_to_stroke > 0))]

mod_stroke <- coxph(Surv(time_accel_to_stroke,incd_stroke_accel) ~ activity_group + age_accel + sex + race_category_adjust +
                      tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status + bmi + sbp + dbp + bpmed + prev_dm_accel,data=stroke_set)
mod_stroke_far <- float(mod_stroke)

mod_stroke_median <- coxph(Surv(time_accel_to_stroke,incd_stroke_accel) ~ activity_group_median + age_accel + sex + race_category_adjust +
                             tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status + bmi + sbp + dbp + bpmed + prev_dm_accel,data=stroke_set)
mod_stroke_far_median <- float(mod_stroke_median)

##################### Step 14: MSK
msk_set <- complete_time[c((prev_msk_accel == 0) & (time_accel_to_msk > 0))]

mod_msk <- coxph(Surv(time_accel_to_msk,incd_msk_accel) ~ activity_group + age_accel + sex + race_category_adjust +
                   tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=msk_set)
mod_msk_far <- float(mod_msk)

mod_msk_median <- coxph(Surv(time_accel_to_msk,incd_msk_accel) ~ activity_group_median + age_accel + sex + race_category_adjust +
                          tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=msk_set)
mod_msk_far_median <- float(mod_msk_median)

# Graphical variables
msk_set[,activity_group_graphical_median := factor(activity_group_median,levels=c('Inactive','Active - Regular','Active - WW'))]

### KM
prodlim_msk <- prodlim(Hist(time_accel_to_msk,incd_msk_accel)~activity_group_graphical_median,data=msk_set)

CairoPDF(file='/weekend_warrior/km_all_msk_median.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
plot(prodlim_msk,"cuminc",ylim=c(0,0.12),xlim=c(0,5),
     axis2.at=seq(0,0.12,0.02),axis2.las=2,lwd=1.4,background=F,
     axis1.at=c(0,1,2,3,4,5),axis1.labels=as.character(0:5),
     atrisk.times=c(0,1,2,3,4,5),col=c("#d95f02",'darkgray','#1b9e77'),atrisk.col='black',confint=TRUE,
     legend.x=0,legend.y=0.12,axis1.cex.axis=2.5,axis2.cex.axis=2.5,axis1.padj=0.5,
     legend.cex=2.2,legend.legend=c("Inactive","Active-Regular","Active-Weekend Warrior"),
     atrisk.title=("                     "),atrisk.pos=-0.5,atrisk.line=c(1.2,2.8,4.4),
     atrisk.cex=1.8,atrisk.interspace=1.4,xlab='',ylab='')
mtext("Cumulative risk (%)",side=2,line=-1.2,at=0.06,cex=2.5)
mtext("Years",side=1, line=-1,cex=2.5)
mtext('Stratum',side=1, line=-0.3,cex=1.8,at=-1.0)
dev.off()

# Graphical variables
msk_set[,activity_group_graphical := factor(activity_group,levels=c('Inactive','Active - Regular','Active - WW'))]

### KM
prodlim_msk <- prodlim(Hist(time_accel_to_msk,incd_msk_accel)~activity_group_graphical,data=msk_set)

CairoPDF(file='/weekend_warrior/km_all_msk.pdf',height=3,width=3.5,
         pointsize=4)
par(oma=c(3,1,1,1),mar=c(3,1,1,1))
plot(prodlim_msk,"cuminc",ylim=c(0,0.12),xlim=c(0,5),
     axis2.at=seq(0,0.12,0.02),axis2.las=2,lwd=1.4,background=F,
     axis1.at=c(0,1,2,3,4,5),axis1.labels=as.character(0:5),
     atrisk.times=c(0,1,2,3,4,5),col=c("#d95f02",'darkgray','#1b9e77'),atrisk.col='black',confint=TRUE,
     legend.x=0,legend.y=0.12,axis1.cex.axis=2.5,axis2.cex.axis=2.5,axis1.padj=0.5,
     legend.cex=2.2,legend.legend=c("Inactive","Active-Regular","Active-Weekend Warrior"),
     atrisk.title=("                     "),atrisk.pos=-0.5,atrisk.line=c(1.2,2.8,4.4),
     atrisk.cex=1.8,atrisk.interspace=1.4,xlab='',ylab='')
mtext("Cumulative risk (%)",side=2,line=-1.2,at=0.06,cex=2.5)
mtext("Years",side=1, line=-1,cex=2.5)
mtext('Stratum',side=1, line=-0.3,cex=1.8,at=-1.0)
dev.off()

##################### Step 14: Incomplete wear data
all_time <- no_etoh
setkey(censor,sample_id); setkey(all_time,sample_id);

## Death handling
all_time[censor, ':='(death_date = i.death_date, death_censor_date = i.death_censor_date)]

all_time[,died := ifelse(!is.na(death_date),1,0)]
all_time[,time_to_death := (ifelse(!is.na(death_date),(as.numeric(death_date) - as.numeric(end_date)),
                                   (as.numeric(death_censor_date) - as.numeric(end_date)))/365.25)]

# Remove small number of people with negative follow up time (4)
all_time <- all_time[!(time_to_death < 0)]

# Age
setkey(all_time,sample_id); setkey(censor,sample_id)
all_time[censor,birth_date := i.birthdate]
all_time[,age_accel := (as.numeric(end_date) - as.numeric(birth_date))/365.25]

# Sex
all_time[censor,sex := i.sex]

# Join in DM
setkey(dm,sample_id); setkey(all_time,sample_id);

all_time[dm,":="(has_dm = i.has_disease, prevalent_dm = i.prevalent_disease, 
                 incident_dm = i.incident_disease, dm_censor_date = i.censor_date)]

all_time[,":="(prev_dm_accel = ifelse(c(!is.na(has_dm) & (has_dm==1) 
                                        & (dm_censor_date <= end_date)),1,0))]


incomplete_time <- all_time[complete_time==0] # 93498 - 89573 = 3925

# Generate phenotypes
## Clean up and prejoin
setkey(af,sample_id); setkey(mi,sample_id); setkey(hf,sample_id)
setkey(stroke,sample_id); setkey(all_time,sample_id);

## Joins
all_time[af,":="(has_af = i.has_disease, prevalent_af = i.prevalent_disease, 
                 incident_af = i.incident_disease, af_censor_date = i.censor_date)]
all_time[mi,":="(has_mi = i.has_disease, prevalent_mi = i.prevalent_disease, 
                 incident_mi = i.incident_disease, mi_censor_date = i.censor_date)]
all_time[hf,":="(has_hf = i.has_disease, prevalent_hf = i.prevalent_disease, 
                 incident_hf = i.incident_disease, hf_censor_date = i.censor_date)]
all_time[stroke,":="(has_stroke = i.has_disease, prevalent_stroke = i.prevalent_disease, 
                     incident_stroke = i.incident_disease, stroke_censor_date = i.censor_date)]

# Fix censoring dates
setkey(all_time,sample_id); setkey(center,sample_id)
all_time[center,':='(center_code = i.value)]

setkey(all_time,center_code); setkey(center_lookup,Code)
all_time[center_lookup,':='(center_location = i.Region)]

# Now correct censor dates based on location
all_time[,':='(af_censor_date = as.Date(ifelse(c(has_af==1 | center_location=='England'),af_censor_date,
                                               ifelse(center_location=='Scotland',pmin(af_censor_date,as.Date('2021-03-31',format='%Y-%m-%d')),
                                                      pmin(af_censor_date,as.Date('2018-02-28',format='%Y-%m-%d')))),origin='1970-01-01'))]
all_time[,':='(mi_censor_date = as.Date(ifelse(c(has_mi==1 | center_location=='England'),mi_censor_date,
                                               ifelse(center_location=='Scotland',pmin(mi_censor_date,as.Date('2021-03-31',format='%Y-%m-%d')),
                                                      pmin(mi_censor_date,as.Date('2018-02-28',format='%Y-%m-%d')))),origin='1970-01-01'))]
all_time[,':='(hf_censor_date = as.Date(ifelse(c(has_hf==1 | center_location=='England'),hf_censor_date,
                                               ifelse(center_location=='Scotland',pmin(hf_censor_date,as.Date('2021-03-31',format='%Y-%m-%d')),
                                                      pmin(hf_censor_date,as.Date('2018-02-28',format='%Y-%m-%d')))),origin='1970-01-01'))]
all_time[,':='(stroke_censor_date = as.Date(ifelse(c(has_stroke==1 | center_location=='England'),stroke_censor_date,
                                                   ifelse(center_location=='Scotland',pmin(stroke_censor_date,as.Date('2021-03-31',format='%Y-%m-%d')),
                                                          pmin(stroke_censor_date,as.Date('2018-02-28',format='%Y-%m-%d')))),origin='1970-01-01'))]

# Define prev and incd with respect to accel
# Prev
all_time[,":="(prev_af_accel = ifelse(c(!is.na(has_af) & (has_af==1) 
                                        & (af_censor_date <= end_date)),1,0),
               prev_mi_accel = ifelse(c(!is.na(has_mi) & (has_mi==1) 
                                        & (mi_censor_date <= end_date)),1,0),
               prev_hf_accel = ifelse(c(!is.na(has_hf) & (has_hf==1) 
                                        & (hf_censor_date <= end_date)),1,0),
               prev_stroke_accel = ifelse(c(!is.na(has_stroke) & (has_stroke==1) 
                                            & (stroke_censor_date <= end_date)),1,0))]
# Incd
all_time[,":="(incd_af_accel = ifelse(c(!is.na(has_af) & (has_af==1) 
                                        & (af_censor_date > end_date)),1,0),
               incd_mi_accel = ifelse(c(!is.na(has_mi) & (has_mi==1) 
                                        & (mi_censor_date > end_date)),1,0),
               incd_hf_accel = ifelse(c(!is.na(has_hf) & (has_hf==1) 
                                        & (hf_censor_date > end_date)),1,0),
               incd_stroke_accel = ifelse(c(!is.na(has_stroke) & (has_stroke==1) 
                                            & (stroke_censor_date > end_date)),1,0))]

# Time
all_time[,":="(time_accel_to_af = as.numeric(af_censor_date - end_date)/365.25,
               time_accel_to_mi = as.numeric(mi_censor_date - end_date)/365.25,
               time_accel_to_hf = as.numeric(hf_censor_date - end_date)/365.25,
               time_accel_to_stroke = as.numeric(stroke_censor_date - end_date)/365.25)]

#### ACTIVITY CLASSES
setkey(all_time,sample_id); setkey(complete_time,sample_id)
all_time[complete_time,':='(activity_group = i.activity_group, activity_group_median = i.activity_group_median)]
all_time[is.na(activity_group)]$activity_group <- 'Inactive'
all_time[is.na(activity_group_median)]$activity_group_median <- 'Inactive'

## Relevels
all_time[,activity_group := factor(activity_group,levels=c('Inactive','Active - WW','Active - Regular'))]
all_time[,activity_group_median := factor(activity_group_median,levels=c('Inactive','Active - WW','Active - Regular'))]

##################### AF
af_set <- all_time[c((prev_af_accel == 0) & (time_accel_to_af > 0))]

mod_af <- coxph(Surv(time_accel_to_af,incd_af_accel) ~ activity_group + age_accel + sex + race_category_adjust +
                  tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=af_set)
mod_af_far <- float(mod_af)

mod_af_median <- coxph(Surv(time_accel_to_af,incd_af_accel) ~ activity_group_median + age_accel + sex + race_category_adjust +
                         tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=af_set)
mod_af_far_median <- float(mod_af_median)

##################### MI
mi_set <- all_time[c((prev_mi_accel == 0) & (time_accel_to_mi > 0))]

mod_mi <- coxph(Surv(time_accel_to_mi,incd_mi_accel) ~ activity_group + age_accel + sex + race_category_adjust +
                  tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=mi_set)
mod_mi_far <- float(mod_mi)

mod_mi_median <- coxph(Surv(time_accel_to_mi,incd_mi_accel) ~ activity_group_median + age_accel + sex + race_category_adjust +
                         tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=mi_set)
mod_mi_far_median <- float(mod_mi_median)

##################### HF
hf_set <- all_time[c((prev_hf_accel == 0) & (time_accel_to_hf > 0))]

mod_hf <- coxph(Surv(time_accel_to_hf,incd_hf_accel) ~ activity_group + age_accel + sex + race_category_adjust +
                  tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=hf_set)
mod_hf_far <- float(mod_hf)

mod_hf_median <- coxph(Surv(time_accel_to_hf,incd_hf_accel) ~ activity_group_median + age_accel + sex + race_category_adjust +
                         tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=hf_set)
mod_hf_far_median <- float(mod_hf_median) 

##################### Stroke
stroke_set <- all_time[c((prev_stroke_accel == 0) & (time_accel_to_stroke > 0))]

mod_stroke <- coxph(Surv(time_accel_to_stroke,incd_stroke_accel) ~ activity_group + age_accel + sex + race_category_adjust +
                      tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=stroke_set)
mod_stroke_far <- float(mod_stroke)

mod_stroke_median <- coxph(Surv(time_accel_to_stroke,incd_stroke_accel) ~ activity_group_median + age_accel + sex + race_category_adjust +
                             tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=stroke_set)
mod_stroke_far_median <- float(mod_stroke_median)

##################### Step 15: Complete true wear time
complete_wear <- complete_time[c(wear_time1 >= 20 & wear_time2 >= 20 & wear_time3 >= 20 & wear_time4 >= 20
                                 & wear_time5 >= 20 & wear_time6 >= 20 & wear_time7 >= 20)]

##################### AF
af_set <- complete_wear[c((prev_af_accel == 0) & (time_accel_to_af > 0))]

mod_af <- coxph(Surv(time_accel_to_af,incd_af_accel) ~ activity_group + age_accel + sex + race_category_adjust +
                  tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=af_set)
mod_af_far <- float(mod_af)

mod_af_median <- coxph(Surv(time_accel_to_af,incd_af_accel) ~ activity_group_median + age_accel + sex + race_category_adjust +
                         tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=af_set)
mod_af_far_median <- float(mod_af_median)

##################### MI
mi_set <- complete_wear[c((prev_mi_accel == 0) & (time_accel_to_mi > 0))]

mod_mi <- coxph(Surv(time_accel_to_mi,incd_mi_accel) ~ activity_group + age_accel + sex + race_category_adjust +
                  tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=mi_set)
mod_mi_far <- float(mod_mi)

mod_mi_median <- coxph(Surv(time_accel_to_mi,incd_mi_accel) ~ activity_group_median + age_accel + sex + race_category_adjust +
                         tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=mi_set)
mod_mi_far_median <- float(mod_mi_median)

##################### HF
hf_set <- complete_wear[c((prev_hf_accel == 0) & (time_accel_to_hf > 0))]

mod_hf <- coxph(Surv(time_accel_to_hf,incd_hf_accel) ~ activity_group + age_accel + sex + race_category_adjust +
                  tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=hf_set)
mod_hf_far <- float(mod_hf)

mod_hf_median <- coxph(Surv(time_accel_to_hf,incd_hf_accel) ~ activity_group_median + age_accel + sex + race_category_adjust +
                         tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=hf_set)
mod_hf_far_median <- float(mod_hf_median) 

##################### Stroke
stroke_set <- complete_wear[c((prev_stroke_accel == 0) & (time_accel_to_stroke > 0))]

mod_stroke <- coxph(Surv(time_accel_to_stroke,incd_stroke_accel) ~ activity_group + age_accel + sex + race_category_adjust +
                      tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=stroke_set)
mod_stroke_far <- float(mod_stroke)

mod_stroke_median <- coxph(Surv(time_accel_to_stroke,incd_stroke_accel) ~ activity_group_median + age_accel + sex + race_category_adjust +
                             tob + tdi + etoh_grams + diet + qual_ea + self_health + employment_status,data=stroke_set)
mod_stroke_far_median <- float(mod_stroke_median)

### Save out final
write.csv(complete_time,file='complete_time.csv',row.names=F)


