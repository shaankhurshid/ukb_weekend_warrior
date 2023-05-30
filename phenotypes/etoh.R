library(data.table)
library(plyr)

etoh_zero <- fread(file='/phenotypes/etoh_visit0.csv')
etoh_one <- fread(file='/phenotypes/etoh_visit1.csv')

# Remove non-responses
etoh_zero <- etoh_0[!(value %in% c(-1,-3))]
etoh_one <- etoh_1[!(value %in% c(-1,-3))]

# Now combine
etoh_all <- rbind(etoh_zero,etoh_one)
etoh_all <- etoh_all[,.SD[which.max(instance)],by='sample_id']
etoh_all$FieldID <- paste0('f_',etoh_all$FieldID)

# And widen
etoh_wide <- dcast(etoh_all,sample_id ~ FieldID, value.var='value')

# Mark people who did not report a frequency of any type
etoh_wide[,':='(freq_summary = ifelse(f_1558==6,"never",
                                      ifelse(c(!is.na(f_1568) | !is.na(f_1578) | !is.na(f_1588) | !is.na(f_1598) | !is.na(f_1608) | !is.na(f_5364)),'weekly',
                                             ifelse(c(!is.na(f_4407) | !is.na(f_4418) | !is.na(f_4429) | !is.na(f_4440) | !is.na(f_4451) | !is.na(f_4462)),'monthly',NA))))]

# Remove NA for frquency (can't calculate)
etoh_wide_nomissing <- etoh_wide[!is.na(freq_summary)]

# Now safely change NA to zero
for (j in c('f_1568','f_1578','f_1588','f_1598','f_1608','f_4407','f_4418',
            'f_4429','f_4440','f_4451','f_4462')){
  set(etoh_wide_nomissing,i=which(is.na(etoh_wide_nomissing[[j]])),j=j,value=0)
}

# Create etoh grams
etoh_wide_nomissing[freq_summary=='weekly',etoh_grams := f_1568*16 + f_1578*16 + f_1588*16 + f_1598*8 + f_1608*16]
etoh_wide_nomissing[freq_summary=='monthly',etoh_grams := (f_4407*16 + f_4418*16 + f_4440*16 + f_4451*8 + f_4462*16)/4]
etoh_wide_nomissing[freq_summary=='never',etoh_grams := 0]

# Rejoin to etoh wide
setkey(etoh_wide_nomissing,sample_id); setkey(etoh_wide,sample_id)
etoh_wide[etoh_wide_nomissing,etoh_grams := i.etoh_grams]
etoh_wide[,f_1558 := ifelse(f_1558==6,'Never',ifelse(f_1558 %in% c(1,2,3),'Weekly','Monthly'))]

#prep for saving
write.csv(etoh_wide,file="/phenotypes/etoh_01.csv",row.names=F)