# Depends
library(data.table)
library(plyr)
library(stringr)

# Read in data
diet_combined <- fread(file='/phenotypes/diet_combined.csv')
salt <- fread(file='/phenotypes/added_salt_all_instance.csv')

# Friendly renames
names(diet_combined)[str_detect(names(diet_combined),'\\_\\d+')] <- paste0('field',names(diet_combined)[str_detect(names(diet_combined),'\\_\\d+')])

# Isolate to instance 0 or 1 for accelerometer analysis
diet_combined <- diet_combined[instance %in% c(0,1)]
salt <- salt[instance %in% c(0,1)]

# Instance 1 if it's there
diet_summarized <- diet_combined[,.SD[which.max(instance)],by='sample_id']
salt_summarized <- salt[,.SD[which.max(instance)],by='sample_id']

# Combine
setkey(salt_summarized,sample_id); setkey(diet_summarized,sample_id)
diet_summarized[salt_summarized,':='(sum_salt = i.value)]

# Recode salt
diet_summarized[,':='(sum_salt = ifelse(sum_salt==-3,NA,sum_salt))]

# Recode
## Helper functions
classify_fruitveg <- function(field){
  ifelse(field=='Do not know',NA,
         ifelse(field=='Prefer not to answer',NA,
                ifelse(field=='Less than one',0.5,field)))
}
classify_meat <- function(field){
  ifelse(field %in% c('Prefer not to answer','','Do not know'),NA,
         ifelse(field=='Never',0,
                ifelse(field=='Less than once a week',1,
                       ifelse(field=='Once a week',2,
                              ifelse(field=='2-4 times a week',3,
                                     ifelse(field=='5-6 times a week',4,
                                            ifelse(field=='Once or more daily',5,6)))))))
}

## Apply functions
diet_summarized[,':='(field_1289 = sapply(field_1289,FUN=classify_fruitveg),
                      field_1299 = sapply(field_1299,FUN=classify_fruitveg),
                      field_1309 = sapply(field_1309,FUN=classify_fruitveg),
                      field_1319 = sapply(field_1309,FUN=classify_fruitveg))]
diet_summarized[,':='(field_1369 = sapply(field_1369,FUN=classify_meat),
                      field_1379 = sapply(field_1379,FUN=classify_meat),
                      field_1389 = sapply(field_1389,FUN=classify_meat),
                      field_1349 = sapply(field_1329,FUN=classify_meat))]

## Cut out lean meat (not used)
diet_summarized[,':='(field_1329 = NULL, field_1339 = NULL, field_1359 = NULL)]

# Cast as numeric
for (j in names(diet_summarized)[str_detect(names(diet_summarized),'field')]){
  set(diet_summarized,j=j,value=as.numeric(diet_summarized[[j]]))
}

# Create summary vars
diet_summarized[,':='(sum_fruitveg = apply(.SD,FUN=sum,MARGIN=1,na.rm=T)),
                .SDcols=c('field_1289','field_1299','field_1309','field_1319')]
diet_summarized[,':='(sum_meat = apply(.SD,FUN=sum,MARGIN=1,na.rm=T)),
                .SDcols=c('field_1369','field_1379','field_1389','field_1349')]

# Go back and set to NA if no answers at all
diet_summarized[,':='(sum_fruitveg = ifelse(is.na(field_1289) & is.na(field_1299)
                                            & is.na(field_1309) & is.na(field_1319),NA,sum_fruitveg),
                      sum_meat = ifelse(is.na(field_1369) & is.na(field_1379) 
                                        & is.na(field_1389) & is.na(field_1349),NA,sum_meat))]

# Create summary vars
diet_summarized[,':='(low_fruitveg = ifelse(sum_fruitveg < median((sum_fruitveg),na.rm=T),1,0),
                      high_meat = ifelse(sum_meat > median((sum_meat),na.rm=T),1,0),
                      high_salt = ifelse(sum_salt > median((sum_salt),na.rm=T),1,0))]
diet_summarized[,':='(diet_quality = ifelse(c(is.na(low_fruitveg) | is.na(high_meat) | is.na(high_salt)),NA,
                                            ifelse(low_fruitveg==1 & c(high_meat==1 | high_salt==1),'poor',
                                                   ifelse(low_fruitveg==0 & high_meat==0 & high_salt==0,'good','intermediate'))))]

# Save out
write.csv(diet_summarized,file='/phenotypes/diet_instance01.csv',row.names=F)



