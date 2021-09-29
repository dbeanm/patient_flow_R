###
# generate some dummy performance data to use in differential network analysis example
###

source("patient_flow_functions.R")

transfers <- read.csv("data/journal.pone.0185912.s003.csv", stringsAsFactors = F)

# get approx dates from transfer data and use these dates for performance data
transfers$SendingFromDate <- calculate_start_of_week(transfers$Week, transfers$Year)
sending_dates <- unique(transfers$SendingFromDate)

#for the sake of an example limit the data to DH
site <- "DH"
transfers <- transfers[which(transfers$Site_source == site | transfers$Site_target == site),] 

#this gives 101 values and there are 85 rows of data
perf_t <- seq(0,10,0.1)
y <- sin(perf_t)

y <- y+1 # sin(x) is in range -1 - +1 a 
y <- y * 50 #get from range 0-2 to 0-100
y <- y[0:length(sending_dates)]

perf <- data.frame('percent_under_4h' = y, 'Date' = sort(sending_dates))

write.csv(perf, 'data/dummy_performance_data.csv', row.names = FALSE)
