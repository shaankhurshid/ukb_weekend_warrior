# Depends
library(data.table)
library(plyr)
library(stringr)

# Read in data
education_age <- fread(file='/phenotypes/age_completed_education_all_instance.csv')
qual <- fread(file='/phenotypes/qualification_all_instance.csv')

# Narrow to instance 0 or 1
education_age <- education_age[instance==0 | instance==1]
qual <- qual[instance==0 | instance==1]

# Convert qual variables to EA
qual[,':='(qual_ea = ifelse(value==1,20,
                            ifelse(value==2,13,
                                   ifelse(c(value==3 | value==4),10,
                                          ifelse(value==5,0,
                                                 ifelse(value==6,15,
                                                        ifelse(value==-7,7,NA)))))))]

# For each individual, pick the max EA that is not NA
qual <- qual[!is.na(qual_ea),.SD[which.max(qual_ea)],by='sample_id']

# Now join years attained
education_age[,':='(value = ifelse(value < 0,NA,value))]
setkey(education_age,sample_id); setkey(qual,sample_id)
qual[education_age,':='(education_age = i.value)]

# Recalculate
qual[,":="(qual_ea = ifelse(qual_ea==0,education_age-5,qual_ea))]

# Write out
write.csv(qual,file='/phenotypes/education_01.csv',row.names=F)

