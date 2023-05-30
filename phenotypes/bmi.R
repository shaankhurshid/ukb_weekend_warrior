# Script to unify BP values

# Dependencies
library(data.table)

# Load
bmi <- fread(file='/phenotypes/bmi_all_instance.csv')
bmi[,FieldID := paste0('f_',FieldID)]

# Combine all instances 0-1 and favor most recent
## To wide
bmi_wide <- dcast(bmi,sample_id ~ FieldID + instance,value.var='value')

## If has manual and automatic, take average
bmi_wide[,bmi := ifelse(!is.na(f_21001_1),f_21001_1,f_21001_0)]
bmi_wide <- bmi_wide[!is.na(bmi)]

# Save out
write.csv(bmi_wide,file='/phenotypes/bmi_instance01.csv',row.names=F)

