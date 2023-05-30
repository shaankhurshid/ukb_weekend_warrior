# Script for assessing calibration of combined clinical/genetic score
classifier <- function(risk,ncuts){
  cuts <- quantile(risk,probs=seq(0,1,1/ncuts),na.rm=T)
  index <- rep(NA,length(risk))
  for (i in 1:(length(cuts)-1)){
    for (j in 1:length(risk)){
      index[j] <- ifelse(risk[j] >= cuts[i],i,index[j])}}
  return(index)
}

# Function to generate survival estimates per AF risk quantile
survivor <- function(data,risk_data,event,time,breakpoint){
  est <- rep(NA,times=length(unique(data[,risk_data])))
  lower <- rep(NA,times=length(unique(data[,risk_data])))
  upper <- rep(NA,times=length(unique(data[,risk_data])))
  level_name <- rep(NA,times=length(unique(data[,risk_data])))
  for (i in 1:length(unique(data[,risk_data]))){
    subset <- data[data[,risk_data]==unique(data[,risk_data])[order(unique(data[,risk_data]))][i],]
    if (nrow(subset[subset[,time] > breakpoint,]) > 0){
      km <- survfit(Surv(subset[,time],subset[,event]) ~ 1, data=subset)
      time_index <- km$time - breakpoint
      end_time <- which(time_index == max(time_index[time_index <= 0]))
      est[i] <- 1-stepfun(km$time[1:end_time], c(1, km$surv[1:end_time]))(breakpoint)
      upper[i] <- 1-stepfun(km$time[1:end_time], c(1, km$lower[1:end_time]))(breakpoint)
      lower[i] <- 1-stepfun(km$time[1:end_time], c(1, km$upper[1:end_time]))(breakpoint)
      level_name[i] <- as.character(unique(data[,risk_data])[order(unique(data[,risk_data]))][i])
      print(level_name[i])
    }
    else {est[i] <- upper[i] <- lower[i] <- NA
    level_name[i] <- as.character(unique(data[,risk_data])[order(unique(data[,risk_data]))][i])
    print(level_name[i])}
  }
  return(data.frame(level=level_name,est=est,upper=upper,lower=lower))
}