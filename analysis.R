### 
# weight variability
# paper figures 2 and 3
###

library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
source("patient_flow_functions.R")

###############################################################
### Load and prepare dataset
###############################################################
# load csv file containing ward transfers with dates, convert
# to a 3-dimensional adjacency matrix of node x node x time

#convert from year-week to year-month
transfers <- read.csv("data/journal.pone.0185912.s003.csv", stringsAsFactors = F)

#for the sake of an example limit the data to DH
site <- "DH"
transfers <- transfers[which(transfers$Site_source == site | transfers$Site_target == site),] 

# sort out dates then convert to a continuous index by year-month
# the first part is very specific to the released flow data format
#first convert from year and week to an actual date
#take the first day of that week as the date, will be aggregated to monthly
#so specific day of week does not matter anyway
# n.b. function is defined in patient_flow_functions.R
transfers$SendingFromDate <- calculate_start_of_week(transfers$Week, transfers$Year)

#Create factor columns for year and month
transfers$SendingMonth <- months(transfers$SendingFromDate)
month_names <- c("January","February","March","April","May","June","July","August","September","October","November","December")
transfers$SendingMonth <- factor(transfers$SendingMonth, levels=month_names)

transfers$SendingYear <- format(transfers$SendingFromDate, format="%Y")
year_names <- c("2014", "2015", "2016")
transfers$SendingYear <- factor(transfers$SendingYear, levels=year_names)

#this part is important for the adjacency matrix functions to work
#whatever aggregation you want to use, you need to create a numeric index for those timepoints
#convert to 1-based month index throughout all timepoints
transfers$SendingIdx <- as.numeric(transfers$SendingMonth) + ((as.numeric(transfers$SendingYear) - 1) * 12)
months_included <- unique(transfers$SendingIdx)
months_included <- months_included[order(months_included, decreasing = F)]

#this step is again very specific to the released flow data
#the adjacency function does NOT currently take into account the n for each of these edges
#it expects 'raw' transfer data at an individual level. So need to repeat each edge N times. 
n_transfers <- transfers$Transfers
transfers <- transfers[rep(seq_len(nrow(transfers)), n_transfers),]

###############################################################
### Analysis - calculations
###############################################################

### make the adjacency matrix

adj <- adjacency_over_time(transfers, "SendingIdx")


## calculate weight variability centrality over time
# weight_variability_arr() applies the weight_variability_centrality() function,
# which calculates the centrality for one node, to all nodes in a single adjacency
# matrix. apply_over_time() applies any function that works on a whole adjacency 
# matrix to all the timepoints from adjacency_over_time().
in_var <- apply_over_time(adj, weight_variability_arr, 'in')
out_var <- apply_over_time(adj, weight_variability_arr, 'out')

## calculate degree centrality over time in the same way
indegree <- apply_over_time(adj, degree_adj, 'in')
outdegree <- apply_over_time(adj, degree_adj, 'out')


###############################################################
### Analysis - example figures
###############################################################

#these first plots are a mess with my data but could theoretically be interesting
#so this is a quick example to plot over time
#also adds variability results to figure_data$in_var etc

#ward colours for all plots
wards <- unique(c(transfers$Source, transfers$Target))
cols <- colorRampPalette( brewer.pal(length(wards),"Set1"))(length(wards))
month_labels <- rep(month_names, length(year_names))[months_included]
figure_data = list()

#input weight variability over time
plot(range(months_included),range(in_var, na.rm = T), type="n", xaxt="n",xlab=NA,main=site, ylab="Input weight variability")
for(i in 1:dim(in_var)[2]){
  lines(months_included,in_var[,i],col=cols[i])
}
axis(side=1, at=months_included, labels = month_labels, las=2)
figure_data[[site]]$in_var <- melt(in_var)

#output weight variability over time
plot(range(months_included),range(out_var, na.rm = T), type="n", xaxt="n",xlab=NA,main=site, ylab="Output weight variability")
for(i in 1:dim(out_var)[2]){
  lines(months_included,out_var[,i],col=cols[i])
}
axis(side=1, at=months_included, labels = month_labels, las=2)
figure_data[[site]]$out_var <- melt(out_var)


###############################################################
### CORE NETWORK - edges present at all timepoints
###############################################################
#this prints the stats for core network used in the paper to terminal
#and makes a couple of quick plots to inspect core network over time

#first some numbers we need later to normalise to
#get total strength
monthly_total_strength <- apply(adj, 3, sum)
#total degree
total_degree <- indegree + outdegree
#edges per month
monthly_total_edges <- rowSums(total_degree)
#nodes per month
present <- total_degree > 0
n_present <- apply(present, 1, sum)


# what are the core edges - the ones that always exist
exists <- adj > 0 #could increase weight required to be considered "core"
always_exists <- matrix(nrow=dim(exists)[1], ncol = dim(exists)[2])
for(i in 1:dim(exists)[1]){
  for(j in 1:dim(exists)[2]){
    always_exists[i,j] <- all(exists[i,j,])
  }
}
cat(length(which(always_exists)), " edges exist at all timepoints\n")

#save core to file
core_edges <- matrix(nrow=length(which(always_exists)), ncol=2)
colnames(core_edges) <- c("Source","Target")
core_row_position <- 1
for(i in 1:dim(always_exists)[1]){
  for(j in 1:dim(always_exists)[2]){
    if(always_exists[i,j]){
      core_edges[core_row_position,] <- c(as.character(wards[i]), as.character(wards[j]))
      core_row_position <- core_row_position + 1
    }
  }
}
write.table(core_edges, file=paste(site,"_core_edges.csv",sep=""),sep=",",row.names = F)

#how much of the total strength each month do these edges account for
core_monthly_strength <- NULL
for(i in 1:dim(adj)[3]){
  x <- adj[,,i]
  core_monthly_strength <-c(core_monthly_strength, sum(x[always_exists]))
}

core_monthly_strength_prop <- core_monthly_strength/monthly_total_strength
plot(months_included,core_monthly_strength_prop, type="o",xaxt="n",xlab=NA, ylab="core as proportion of total transfers",main=site)
axis(side=1, at=months_included, labels = month_labels, las=2)
cat("mean proportion of transfers through core edges:",mean(core_monthly_strength_prop),"sd",sd(core_monthly_strength_prop),"\n")

hist(core_monthly_strength_prop, main=site, xlab="Proportion of transfers through core edges")

#what proportion of all edges in the network are the "core" edges per month
monthly_prop_core <- length(which(always_exists))/monthly_total_edges
cat("mean proportion of edges in code",mean(monthly_prop_core),"sd",sd(monthly_prop_core),"\n")
plot(months_included,monthly_prop_core, type="o",xaxt="n",xlab=NA, ylab="core as proportion of total edges",main=site)
axis(side=1, at=months_included, labels = month_labels, las=2)

#proportion of monthly nodes that are in the core
core_n_nodes <- length(unique(c(core_edges[,1], core_edges[,2])))
cat("nodes in core:",core_n_nodes,"\n")
monthly_nodes_prop_core <- core_n_nodes/n_present
cat("mean proportion of nodes in core",mean(monthly_nodes_prop_core),"sd",sd(monthly_nodes_prop_core),"\n")

###############################################################
### Degree plot (paper figure 2)
### note that we just have the DH data twice here
###############################################################
#pdf(file="degree difference over time.pdf",width=10,height=10)

#data for the indegree-vs-outdegree line
indeg_minus_outdeg <- indegree - outdegree
figure_data[[site]]$indeg_minus_outdeg = melt(data.frame(indeg_minus_outdeg))

#data for indegree-vs-outdegree scatter
max_poss_degree <- length(wards)-1 #max either in or out, total is 2*this
min_poss_degree <- 0 #a node could have no edges for some months
deg_range <- c(min_poss_degree:max_poss_degree)
intensity_matrix <- matrix(data = 0, nrow=length(deg_range), ncol = length(deg_range))
#indegree values in columns, outdegree values in rows
for(i in 1:dim(indegree)[1]){
  for(j in 1:dim(indegree)[2]){
    ind <- indegree[i,j]
    outd <- outdegree[i,j]
    #+1 to value because arrays are 1-indexed but value can be zero
    intensity_matrix[outd+1, ind+1] <- intensity_matrix[outd+1, ind+1] + 1
  }
}

ind_vec <- as.vector(indegree)
out_vec <- as.vector(outdegree)
intensity <- rep(0,length(ind_vec))
for(i in 1:length(ind_vec)){
  intensity[i] <- intensity_matrix[out_vec[i] + 1, ind_vec[i] + 1]
}
df <- data.frame('indegree' = ind_vec, 'outdegree' = out_vec, 'frequency' = intensity)
#remove point 0,0 - this corresponds do nodes that didn't have any edges that month which isn't relevant
df <- df[-which(df$indegree == 0 & df$outdegree == 0),]
figure_data[[site]]$in_vs_out <- df

#dates
mn <- paste(month_names, rep(year_names, each=12))[months_included]
figure_data$DH$indeg_minus_outdeg$Date = factor(mn,levels = mn)
figure_data$PRUH <- figure_data$DH #for demo just copy over PRUH data

#DH colour scale
wards <- unique(figure_data$DH$indeg_minus_outdeg$variable)
first <- rep(0,length(wards))
for(i in 1:length(wards)){
  w <- wards[i]
  first[i] <- figure_data$DH$indeg_minus_outdeg$value[which(figure_data$DH$indeg_minus_outdeg$variable == w & figure_data$DH$indeg_minus_outdeg$Date == mn[1])]
}
names(first) <- wards
first <- first[order(abs(first), decreasing = T)]
myColors_dh <- rep("#000000",length(wards))
myColors_dh[1:3] <- brewer.pal(3,"Dark2")
names(myColors_dh) <- names(first)
highlight_dh <- names(myColors_dh)[1:4]
labels_dh <- highlight_dh
labels_dh[4] <- "Other wards"
labels_dh[which(labels_dh == "KCH.EmergencyDept")] <- "A&E"
labels_dh[which(labels_dh == "ClinicalDecisionUnit.A.E..DH")] <- "ClinicalDecisionUnit"
dh_col_scale <- scale_colour_manual(name = "Ward",breaks=highlight_dh, labels=labels_dh, values = myColors_dh)
pruh_col_scale <- dh_col_scale #the pruh data is just a copy of dh for this demo

#titles
pruh_title <- expression(paste(bold("a"),italic(" PRUH")))
pruh_title_sc <- expression(paste(bold("c"),italic(" PRUH")))
dh_title <- expression(paste(bold("b"),italic(" DH")))
dh_title_sc <- expression(paste(bold("d"),italic(" DH")))

#pruh over time
pruh_inout <- ggplot(data=figure_data$PRUH$indeg_minus_outdeg, aes(x=Date, y=value, group=variable, colour=variable)) + geom_line()
pruh_inout <- pruh_inout + theme(axis.title.y=element_text(size = 8),axis.text.x=element_text(angle=45, hjust=1), plot.title=element_text(hjust=0)) + labs(y="In-degree - out-degree", title=pruh_title)
pruh_inout <- pruh_inout + pruh_col_scale# + guides(colour=FALSE)

#dh over time
dh_inout <- ggplot(data=figure_data$DH$indeg_minus_outdeg, aes(x=Date, y=value, group=variable, colour=variable)) + geom_line()
dh_inout <- dh_inout + theme(axis.title.y=element_text(size = 8),axis.text.x=element_text(angle=45, hjust=1), plot.title=element_text(hjust=0)) + labs(y="In-degree - out-degree", title=dh_title)
dh_inout <- dh_inout + dh_col_scale #+ guides(colour=FALSE)

#pruh in vs out
pruh_scatter <-ggplot(figure_data$PRUH$in_vs_out, aes(x=indegree, y=outdegree)) + geom_point()
pruh_scatter <- pruh_scatter + labs(title=pruh_title_sc) + guides(color=F) + theme(plot.title=element_text(hjust=0))

#dh in vs out
dh_scatter <-ggplot(figure_data$DH$in_vs_out, aes(x=indegree, y=outdegree)) + geom_point()
dh_scatter <- dh_scatter + labs(title=dh_title_sc) + guides(color=F)+ theme(plot.title=element_text(hjust=0))

grid.arrange(grobs=list(pruh_inout, dh_inout,pruh_scatter,dh_scatter), layout_matrix=matrix(c(1,1,2,2,3,4),ncol=2, byrow = T))
#dev.off()

###############################################################
### input and output weigth variability histograms (paper figure 3)
### note that we just have the DH data twice here
###############################################################
#in
figure_data$PRUH$in_var$Site <- "PRUH"
figure_data$PRUH$in_var$Type <- "Input"
figure_data$DH$in_var$Site <- "DH"
figure_data$DH$in_var$Type <- "Input"
#out
figure_data$PRUH$out_var$Site <- "PRUH"
figure_data$PRUH$out_var$Type <- "Output"
figure_data$DH$out_var$Site <- "DH"
figure_data$DH$out_var$Type <- "Output"
#combine
all_vars <- rbind(figure_data$PRUH$in_var, figure_data$PRUH$out_var, figure_data$DH$out_var, figure_data$DH$in_var)
all_vars <- all_vars[-which(is.na(all_vars$value)),]
var_plot <- ggplot(data=all_vars, aes(x=value)) + geom_histogram() + facet_grid(Type ~ Site, scales = "free_y") + labs(x="Weight variability score")
#pdf(file="weight variability histograms.pdf", width=5,height=5)
print(var_plot)
#dev.off()