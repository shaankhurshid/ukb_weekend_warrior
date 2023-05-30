# Script to unify BP values

# Dependencies
library(data.table)

# Load
bpmed0a <- fread(file='/phenotypes/med_6153_instance0.csv')
bpmed1a <- fread(file='/phenotypes/med_6153_instance1.csv')

bpmed0b <- fread(file='/phenotypes/med_6177_instance0.csv')
bpmed1b <- fread(file='/phenotypes/med_6177_instance1.csv')

# Combine
bpmed_a <- rbind(bpmed0a,bpmed1a)
bpmed_b <- rbind(bpmed0b,bpmed1b)

# For each person, filter to the most recent instance with data
bpmed_a_latest <- bpmed_a[,.SD[which.max(instance)],by='sample_id']
bpmed_b_latest <- bpmed_b[,.SD[which.max(instance)],by='sample_id']
on_bpmed <- unique(c(bpmed_a_latest[value==2]$sample_id,bpmed_b_latest[value==2]$sample_id))

# Now just filter on whether BP med value exists
bp_summary <- data.table(sample_id = unique(c(bpmed_a_latest$sample_id,bpmed_b_latest$sample_id)))
bp_summary[,prev_bpmed := ifelse(sample_id %in% on_bpmed,1,0)]

# Save out
write.csv(bp_summary,file='/phenotypes/bpmed_combined_instance01.csv',row.names=F)