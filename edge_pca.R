###
# PCA of edge weights
###

source("patient_flow_functions.R")
library(ggplot2)
library(grid)
library(gridExtra)

transfers <- read.csv("data/journal.pone.0185912.s003.csv", stringsAsFactors = F)

#parameters for the analysis
#filter out edges present in less than DAY_THRESH proportion of all days
DAY_THRESH <- 0.01 #0.1 i.e. 10% of all days

#for the sake of an example limit the data to DH
SITE <- "DH"
transfers <- transfers[which(transfers$Site_source == SITE | transfers$Site_target == SITE),] 

# sort out dates then convert to a continuous index by year-month
# the first part is very specific to the released flow data format
#first convert from year and week to an actual date
#take the first day of that week as the date, will be aggregated to monthly
#so specific day of week does not matter anyway
transfers$SendingFromDate <- calculate_start_of_week(transfers$Week, transfers$Year)
sending_dates <- unique(transfers$SendingFromDate)

## work out all edges that exist
transfers$key <- paste(transfers$Source, transfers$Target, sep="_")

all_edge_types <- unique(transfers$key)

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


#normalise weights to a proportion of all transfers that day
day_total_tf <- rowSums(edge_weights)
for(i in 1:dim(edge_weights)[1]){
  edge_weights[i,] <- edge_weights[i,]/day_total_tf[i]
}
##end of normalisation chunk

#filter out any edges that are present in less than THRESH % of days
n_days <- dim(edge_weights)[1]
keep <- rep(T,dim(edge_weights)[2])
for(i in 1:dim(edge_weights)[2]){
  count <- sum(edge_weights[,i] > 0)
  prop_days_present <- count/n_days
  if(prop_days_present < DAY_THRESH){
    keep[i] <- FALSE
  }
}
cat("with threshold",DAY_THRESH,"removing", sum(!keep),"edges before PCA, leaving",sum(keep),"\n")
edge_weights <- edge_weights[,keep]


edge.pca <- prcomp(edge_weights, center = T, scale = T)
plot(edge.pca, type="l")

#save the pca object
save(list=c("edge.pca", "edge_weights"), file=paste("edge.pca_",SITE,"_thresh",DAY_THRESH,".RData",sep=""))
#save the first 4 PC to file
pc_out_max <- min(dim(edge.pca$x)[2], 4) #output up to 4 pc if there are more than 4
pc_output <- data.frame("Date" = sending_dates, edge.pca$x[,1:pc_out_max]) 
write.table(pc_output, paste("pca_components_",SITE,"_thresh",DAY_THRESH,".csv",sep=""),row.names = F, sep=",")

#dataframe for plotting
plotdata <- data.frame(edge.pca$x)


#weekday vs weekend is not possible in the released data but this is just an example to show
#plotting the PCA results with days highlighted by a class. Real weekday/weekend data, A&E performance
#data etc. would all be good candidates to use.
day_of_week <- weekdays(sending_dates)
day_type <- rep("Weekday",length(day_of_week))
#first line is how to assign actual weekends. Since our dates are all mondays, assign a few days to weekend.
#day_type[day_of_week %in% c("Saturday", "Sunday")] <- "Weekend"
nwd <- floor(length(day_type)/7)
day_type[c(1:nwd)] <- "Weekend"

day_type <- factor(day_type, levels=c("Weekday", "Weekend"))
plotdata$day_type <- day_type
save(list=c("plotdata"), file=paste("plotdata_daytype_",SITE,"_thresh",DAY_THRESH,".RData",sep=""))

#creat scatter plots for all pairs of PC up to max_pc with days highlighted by day_type
max_pc <- 4
for(i in 1:(max_pc-1)){
  for(j in (i+1):max_pc){
    tmp <- data.frame("x" = plotdata[,i], "y" = plotdata[,j], "day_type" = plotdata$day_type)
    lab_x <- paste("PC",i,sep="")
    lab_y <- paste("PC",j,sep="")
    pcp <- ggplot(data=tmp, aes(x=x, y=y, colour=day_type)) + geom_point(alpha=0.7) + labs(title=SITE, x=lab_x, y=lab_y)
    print(pcp)
  }  
}

#look for a PC that separates the day type in the first 6 PC
plts <- NULL
for(i in 1:6){
  tmp <- data.frame("Component_x" = plotdata[,i], "day_type" = plotdata$day_type)
  p <- ggplot(data=tmp, aes(x=day_type, y=Component_x, fill=day_type)) + geom_boxplot() + labs(y=paste("PC",i,sep="")) + theme(legend.position="none")
  plts[[i]] <-p
}
grid.arrange(grobs=plts, top = paste("Separation of groups by PC (",SITE,")", sep=""),layout_matrix = matrix(c(1,2,3,4,5,6), ncol=2, byrow=TRUE))

#extract the most important edges for any components of interest
#e.g. to get PCs 1-2 (keep_pc)
keep_pc <- c(1,2)
for(pc in keep_pc){
  pc.rot <- edge.pca$rotation[,pc]
  pc.rot <- pc.rot[order(abs(pc.rot), decreasing = T)]
  #look at the distribution of the weights for the top 6 edges in PC (positive or negative)
  pc.sig <- names(pc.rot[1:6])
  w <- which(colnames(edge_weights) %in% pc.sig)
  edge_weights_pc <- data.frame(edge_weights[,w])
  edge_weights_pc$day_type <- day_type
  
  #boxplots
  plts <- NULL
  for(i in 1:(dim(edge_weights_pc)[2]-1)){ #there's column of factors at the end
    tmp <- data.frame("edge" = edge_weights_pc[,i], "day_type" = edge_weights_pc$day_type)
    tt <- gsub("_"," to ",colnames(edge_weights_pc)[i])
    extra <- max(tmp$edge)*0.2
    ymax <- max(tmp$edge) + extra
    bp <- ggplot(data=tmp, aes(y=edge, x=day_type, fill=day_type)) + geom_boxplot() + labs(title=tt,y="Weight") + coord_cartesian(ylim=c(0,ymax))+theme(legend.position="none")
    #ttests
    a_b <- t.test(tmp$edge[tmp$day_type == "Weekday"],tmp$edge[tmp$day_type == "Weekend"])$p.value
    a_b_s <- p_to_str(a_b)
    
    line_low <- signif(max(tmp$edge) + (extra * 0.1), 3)
    
    l <- "---"
    r <- "---"
    bp <- bp + annotate("text", x = 1.5, y = line_low, label = paste(l,r,sep=""), size = 5)
    bp <- bp + annotate("text", x = 1.5, y = line_low, label = a_b_s, size = 5)
    
    #print(bp)
    plts[[i]] <- bp
  }
  #multiplot(plotlist=plts, cols =2)
  grid.arrange(grobs=plts, top = paste("Top edges for PC",pc," (",SITE,")",sep=""),layout_matrix = matrix(c(1,2,3,4,5,6), ncol=2, byrow=TRUE))
}

