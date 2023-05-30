# Script to unify BP values

# Dependencies
library(data.table)

# Load
sbp <- fread(file='/phenotypes/sbp_combined_all.csv')
dbp <- fread(file='/phenotypes/dbp_combined_all.csv')

# Combine all instances 0-2 and favor most recent
## To wide
sbp <- dcast(sbp,sample_id ~ FieldID + instance,value.var='value',fun.aggregate = mean)

dbp <- dcast(dbp,sample_id ~ FieldID + instance,value.var='value',fun.aggregate = mean)

## If has manual and automatic, take within-instance average
sbp[,sbp1 := apply(.SD,FUN=mean,na.rm=T,MARGIN=1),.SDcols=c('93_1','4080_1')]
sbp[,sbp0 := apply(.SD,FUN=mean,na.rm=T,MARGIN=1),.SDcols=c('93_0','4080_0')]

dbp[,dbp1 := apply(.SD,FUN=mean,na.rm=T,MARGIN=1),.SDcols=c('94_1','4079_1')]
dbp[,dbp0 := apply(.SD,FUN=mean,na.rm=T,MARGIN=1),.SDcols=c('94_0','4079_0')]

## Now, prioritize the most recent value
sbp[,sbp_combined := ifelse(is.na(sbp1),sbp0,sbp1)]
dbp[,dbp_combined := ifelse(is.na(dbp1),dbp0,dbp1)]

sbp <- sbp[,c('sample_id','sbp_combined')]
dbp <- dbp[,c('sample_id','dbp_combined')]

# Save out
write.csv(sbp,file='/phenotypes/sbp_combined_instance01.csv',row.names=F)
write.csv(dbp,file='/phenotypes/dbp_combined_instance01.csv',row.names=F)

