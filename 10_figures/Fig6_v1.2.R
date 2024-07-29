# Author: M Jager
# Date: December 2023
# Version: 1.2
# R version 4.2.2




#### START UP ####
outdir = "/path/to/NanoRCS/output/output_figures/"
indir_data = "/path/to/NanoRCS/10_figures/source_data/"

#2 Libraries
library(ggplot2) #Version: ‘3.4.0’
library(reshape2) #Version: ‘1.4.4’




#### Get and combine the data #####

#1 Get the timeline table of all events
timeline <- read.table(file = paste(indir_data,"Fig6.txt",sep=""), sep = "\t",header = TRUE)

#2 Get the NanoRCS TF output
nano <- read.table(file = paste(indir_data,"Fig6_2.txt",sep=""), sep = "\t",header = TRUE)

#3 Create empty columns for NanoRCS measurements in timeline
timeline$NanoRCS_SNV <- NA
timeline$NanoRCS_SNV_min <- NA
timeline$NanoRCS_SNV_max <- NA
timeline$NanoRCS_CNV <- NA
timeline$NanoRCS_fragment <- NA

#4 Add the NanoRCS data to the timeline
for(rownr in 1:nrow(nano)) {
  sample <- nano[rownr,]$Sample
  sample <- paste("B",as.numeric(tail(unlist(strsplit(sample,split="-B")),n=1))-1,sep="")
  timeline[which(timeline$label == sample),]$NanoRCS_SNV <- nano[rownr,]$NanoRCS_SNV
  timeline[which(timeline$label == sample),]$NanoRCS_SNV_min <- nano[rownr,]$NanoRCS_SNV_min
  timeline[which(timeline$label == sample),]$NanoRCS_SNV_max <- nano[rownr,]$NanoRCS_SNV_max
  timeline[which(timeline$label == sample),]$NanoRCS_CNV <- nano[rownr,]$NanoRCS_CNV
  timeline[which(timeline$label == sample),]$NanoRCS_fragment <- nano[rownr,]$NanoRCS_fragment
  remove(sample)
}
remove(rownr)

#5 Remove NanoRCS table
remove(nano)




#### Prepare a dataframe for plotting ####

#1 Get lines from timeline with blood measurements
plotdf <- timeline[which(timeline$type == "blood"),]

#2 Convert to a nice table for ggplot
plotdf <- data.frame(plotdf$day,
                     plotdf$ddPCR,
                     plotdf$NanoRCS_SNV,
                     plotdf$NanoRCS_CNV,
                     plotdf$NanoRCS_fragment
)
colnames(plotdf)<- c("day","ddPCR","NanoRCS_SNV","NanoRCS_CNV","NanoRCS_fragment")
plotdf <-melt(plotdf, id = c('day'))
colnames(plotdf) <- c("day","Technique","value")
plotdf <- plotdf[!is.na(plotdf$value),]

#3 For ddPCR: add which timepoints are not covered in NanoRCS
plotdf$Technique_covered <- paste(plotdf$Technique,"only",sep="_")
days <- plotdf[which(plotdf$Technique == "NanoRCS_SNV"),]$day
for(day in days) {
  plotdf[which(plotdf$day == day),]$Technique_covered <- as.character(plotdf[which(plotdf$day == day),]$Technique)
}
remove(day, days)

#4 Add error bar data for NanoRCS_SNV
plotdf$min = NA
plotdf$max = NA
for(line in 1:nrow(plotdf)) {
  if(plotdf[line,]$Technique == "NanoRCS_SNV") {
    plotdf[line,]$min <- timeline[which(timeline$day == plotdf[line,]$day),]$NanoRCS_SNV_min[!is.na(timeline[which(timeline$day == plotdf[line,]$day),]$NanoRCS_SNV_min)]
    plotdf[line,]$max <- timeline[which(timeline$day == plotdf[line,]$day),]$NanoRCS_SNV_max[!is.na(timeline[which(timeline$day == plotdf[line,]$day),]$NanoRCS_SNV_max)]
  }
}
remove(line)




#### Create a dataframe with disease and treatment status (for plotting) ####

#1 Create a temporary table with all disease status data
tempdf <- timeline[which(timeline$type == "disease"),]

#2 Transfer this 'disease status data' to a table with 'rectangles' to plot the disease status in the final plot
rects <- data.frame(start = c(tempdf$day -5), 
                    end = c(tempdf$day + 5),
                    label = tempdf$label)
rects[which(rects$start < 0),]$end <- rects[which(rects$start < 0),]$end+12
rects[which(rects$start < 0),]$start <- 0

#3 Create a temporary table with all treatment data
tempdf <- timeline[which(timeline$type == "chemotherapy" | timeline$type == "surgery"),]
tempdf <- tempdf[which(tempdf$day >= 0),]
tempdf <- tempdf[which(tempdf$date <= "2020-03-03"),]

#4 Add this 'treatment data' to the rect table
rects <- rbind(rects,
               data.frame(start =c(tempdf$day[-nrow(tempdf)]), 
                 end = c(tempdf$day[-nrow(tempdf)]),
                 label = tempdf$label[-nrow(tempdf)]))

#5 Add colour that will be used for plotting
rects$group = NA
for(line in 1:nrow(rects)) {
  if(length(grep("S",rects[line,]$label))>0) {rects[line,]$group <- "Stripe"}
  if(length(grep("C",rects[line,]$label))>0) {rects[line,]$group <- "Dot"}
  if(rects[line,]$label == "PD") {rects[line,]$group <- "Red"}
  if(rects[line,]$label == "SD") {rects[line,]$group <- "Orange"}
  if(rects[line,]$label == "PR") {rects[line,]$group <- "Green"}
  if(rects[line,]$label == "DD") {rects[line,]$group <- "Black"}
}
remove(line)

#6 Take time frame relevant for NanoRCS study
rects <- rects[which(rects$end < 605),]

#7 Cleanup
remove(tempdf)




#### Plot ####
#1 Plot
plot <- ggplot(plotdf,aes(x=day, y= value, group = Technique)) +
  ggtitle("Timeline GCT02") +
  theme_bw() +
  # Axes
  ylab("ctDNA fraction") +
  scale_y_continuous(expand = c(0, 0)) +
  expand_limits(y=0.65) +
  xlab("Time after inclusion in the study (days)") +
  scale_x_continuous(expand = c(0, 0)) +
  expand_limits(x = 600) +
  # Disease
  geom_rect(data=rects[which(rects$group == "Red"),], inherit.aes=FALSE, aes(xmin=start, xmax=end, ymin=0,ymax=0.65, group=group), color="transparent", fill="firebrick3", alpha=0.4) +
  geom_rect(data=rects[which(rects$group == "Orange"),], inherit.aes=FALSE, aes(xmin=start, xmax=end, ymin=0,ymax=0.65, group=group), color="transparent", fill="orange2", alpha=0.4) +
  geom_rect(data=rects[which(rects$group == "Green"),], inherit.aes=FALSE, aes(xmin=start, xmax=end, ymin=0,ymax=0.65, group=group), color="transparent", fill="dodgerblue3", alpha=0.4)+
  # Treatment
  geom_vline(xintercept = rects[which(rects$group == "Stripe"),]$start,lty=2, colour = "grey")+
  geom_vline(xintercept = rects[which(rects$group == "Dot"),]$start,lty=3, colour = "grey") +
  # ddPCR and NanoRCS data
  geom_line(aes(color = Technique, linetype = Technique), linewidth = 0.8) +
  scale_color_manual(values=c("grey30","grey30","dodgerblue3","steelblue1"))+
  scale_linetype_manual(values = c(3,1,1,1))+
    geom_errorbar(aes(ymin = min,ymax=max), color = "grey30", linewidth = 0.35)+
  geom_point(aes(shape = Technique_covered, color = Technique),size = 2)+
  scale_shape_manual(values=c(19,17,19,19,19))
#2 Save plot
ggsave(filename = paste(outdir,"GCT02_raw.pdf",sep=""),plot = plot, width = 8.2, height = 3)

