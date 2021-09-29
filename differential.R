###
# save the best and worst networks
###

source("patient_flow_functions.R")
library(ggplot2)

###############################################################
### Load and prepare dataset
###############################################################
# load csv file containing ward transfers with dates, convert
# to a 3-dimensional adjacency matrix of node x node x time

#convert from year-week to year-month
transfers <- read.csv("data/journal.pone.0185912.s003.csv", stringsAsFactors = F)

#usual setup to get the released data to look like raw data
#(assign an approximate date, repeat each transfer n times, limit to DH)
SITE <- "DH"
transfers <- transfers[which(transfers$Site_source == SITE | transfers$Site_target == SITE),] 
transfers$SendingFromDate <- calculate_start_of_week(transfers$Week, transfers$Year)
n_transfers <- transfers$Transfers
transfers <- transfers[rep(seq_len(nrow(transfers)), n_transfers),]

#performance data generated in generate_performance_data.R
perf <- read.csv('data/dummy_performance_data.csv')
perf$Date <- as.Date(perf$Date)



###############################################################
### parameters
###############################################################
#number of sd away from the mean to be considered significant
DIFF_SD <- 2

#restrict to edges present in x% of days
DAY_THRESH <- 0.01 #0.1 i.e. 10% of all days

#top and bottom x% of days considered extremes of performance
EXTREME <- 0.1 #0.1 = top 10% vs bottom 10%



###############################################################
### Find the best vs worst 
###############################################################
n_extreme <- floor(dim(perf)[1] * EXTREME)
extr_perf <- perf[order(perf$percent_under_4h, decreasing = T),]
class3 <- rep("Mid",dim(extr_perf)[1])
class3[1:n_extreme] <- "High"
class3[(length(class3) - n_extreme):length(class3)] <- "Low"
extreme_dates <- extr_perf$Date[class3 != "Mid"]
top_dates <- extr_perf$Date[class3 == "High"]
low_dates <- extr_perf$Date[class3 == "Low"]

perf <- perf[perf$Date %in% extreme_dates,]
transfers <- transfers[transfers$SendingFromDate %in% extreme_dates,]

###############################################################
### Get networks for best and worst days
###############################################################
#get a matrix of the weights of all edges in all the days we've identified as
#having extreme performance
# edge types are rows, dates are columns
#this is the same data as in the adjacency matrix but flattened to 2 dimensions 
#(source+target, date) rather than 3 (source, target, date)

sending_dates <- unique(transfers$SendingFromDate)

transfers$key <- paste(transfers$Source, transfers$Target, sep="_")

#this is all the unique transfers between pairs of wards that ever happened on our extreme days
all_edge_types <- unique(transfers$key)

#make an empty matrix to store the data then work through all the extreme days
#and get the number of transfers on that day
edge_weight_matrix <- matrix(nrow=length(all_edge_types), ncol=length(sending_dates))
rownames(edge_weight_matrix) <- all_edge_types
colnames(edge_weight_matrix) <- sending_dates
for(j in 1:length(sending_dates)){
  d <- sending_dates[j]
  for(i in c(1:length(all_edge_types))){
    k <- all_edge_types[i]
    edge_weight_matrix[i,j] <- length(which(transfers$key[transfers$SendingFromDate == d] == k))
  }
}

edge_weights <- t(edge_weight_matrix)

###############################################################
### identify differential patient flows between best and worst networks
###############################################################
#we want to analyse normalised flows i.e. transfers through each edge as a proportion of all transfers
#to do that we first work out the totals per group
is_best <- sending_dates %in% top_dates 
sum_best <- colSums(edge_weights[is_best,]) #total transfers through each edge over all best days
sum_worst <- colSums(edge_weights[!is_best,])
total_out <- rbind(sum_best, sum_worst)
rownames(total_out) <- c("best","worst")

#now we can normalise each edge to a proportion of transfers
#divide by total number of transfers in each group (best or worst) to get proportion
normalised_flow <- total_out
normalised_flow[1,] <- normalised_flow[1,]/sum(normalised_flow[1,]) #best
normalised_flow[2,] <- normalised_flow[2,]/sum(normalised_flow[2,]) #worst

#what is the difference in average flow between best and worst days
normalised_flow <- rbind(normalised_flow, normalised_flow[1,] - normalised_flow[2,])
rownames(normalised_flow)[3] <- "diff"
normalised_flow <- t(normalised_flow)

#now we have the differences, the next section will label them as a significant
#increase or decrease

###############################################################
### build the main output dataframe
###############################################################

#put the source and target names back to build the output table
sources <- rep(NA, dim(normalised_flow)[1])
targets <- rep(NA, dim(normalised_flow)[1])
edge_names <- rownames(normalised_flow)
for(n in 1:length(edge_names)){
  x <- unlist(strsplit(edge_names[n],"_")) #at the start we labelled edges as "source_target" so split back
  sources[n] <- x[1]
  targets[n] <- x[2]
}
normalised_flow <- data.frame(normalised_flow)
normalised_flow$source <- sources
normalised_flow$target <- targets

#absolute difference is useful when visualising these networks
#you can scale edge thickness to absolute difference and colour by direction of change (positive, negative)
normalised_flow$abs_diff <- abs(normalised_flow$diff) 

#add labels for change in flow based on significance threshold
#initially assign everything to the unchanged group since this is likely to be most edges
#then later add labels for positive/negative change
normalised_flow$edge_bin <- rep("SAME",dim(normalised_flow)[1])

#use standard deviation to set a threshold for a significant change
m <- mean(normalised_flow$diff)
s <- sd(normalised_flow$diff)
high_thresh <- m + (DIFF_SD*s)
low_thresh <- m - (DIFF_SD*s)

normalised_flow$edge_bin[normalised_flow$diff > high_thresh] <- "POSITIVE"
normalised_flow$edge_bin[normalised_flow$diff < low_thresh] <- "NEGATIVE"

#save the results to file, ready for other analysis/to visualise
fname <- paste("edges_pooled_",SITE,"_bw_all_norm.csv", sep="")
write.table(normalised_flow, fname, sep=",", row.names = F, quote=F)

###############################################################
### other output 
###############################################################

#histogram of fold change in edge weights in best vs worst, blue lines represent cutoff for significance
#expect this to have most edges around 0 i.e. not that different best vs worst on average. 
diff_hist <- ggplot(data = normalised_flow, aes(x=diff)) + geom_histogram(bins = 40) + geom_vline(xintercept = c(low_thresh, high_thresh), colour="blue")
diff_hist <- diff_hist + labs(title=SITE)
print(diff_hist)

#print some summary info to console
cat("for",SITE,"\n")
cat(SITE,"total transfers in best",sum(sum_best),"worst",sum(sum_worst),"\n")
 cat("high thresh")
print(high_thresh)
cat("low thresh")
print(low_thresh)
cat("mean norm diff",mean(normalised_flow$diff),"sd",sd(normalised_flow$diff),"\n")
cat(dim(normalised_flow)[1],"edges total \n")
cat(length(which(normalised_flow$edge_bin != "SAME")),"edge weights above/below threshold",length(which(normalised_flow$edge_bin == "POSITIVE")),"above,",length(which(normalised_flow$edge_bin == "NEGATIVE")),"below \n")


###############################################################
### END
###############################################################

