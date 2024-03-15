# Author: M Jager
# Date: December 2023
# Version: 1.3
# R version 4.2.2




#### START UP ####

#1 Directories
#1A Main directory
output_dir = "/path/to/NanoRCS/output/output_figures/"
#1B Input
indir = "/path/to/NanoRCS/output/processed_output/07_PCAWG_TF_limit_of_detection/combined/"
#1C Output
outdir = paste(output_dir,"Suppl Fig 9/raw/",sep="")

#2 Libraries
library(ggplot2) #Version: ‘3.4.0’

#3 Functions
#3A Add the techniques to the df
add.technique.todf <- function(df) {
  # Add a column with the technique to the df
  df$technique = paste(df$sequencing_error_rate," ",df$sequencing_depth,"x",sep="")
  # Add platforms (based on error rate)
  df$technique = gsub("0.00082","NanoRCS",df$technique)
  df$technique = gsub("0.00744","Nanopore",df$technique)
  # Add techniques (based on error rate and coverage)
  df$technique = gsub("NanoRCS 0.8x","NanoRCS Promethion",df$technique)
  df$technique = gsub("NanoRCS 0.04x","NanoRCS Minion",df$technique)
  df$technique = gsub("Nanopore 0.25x","Minion",df$technique)
  df$technique = gsub("Nanopore 3x","Promethion",df$technique)
  # Change order of columns
  df <- cbind(df[,ncol(df)],df[,1:(ncol(df)-1)])
  colnames(df)[1] <- "technique"
  # Return the df
  return(df)
}

#3B Subset df 
subset.df <- function(df,subset) {
  df_selected <- data.frame()
  df_selected = df[which(df$technique == "NanoRCS Minion"),]
  df_selected  = rbind(df_selected ,df[which(df$technique == "NanoRCS Promethion"),])
  df_selected  = rbind(df_selected ,df[which(df$technique == "Minion"),])
  df_selected  = rbind(df_selected ,df[which(df$technique == "Promethion"),])
  return(df_selected)
}

#3C Add detection yes/no
add.detection <- function(df,fdr,tpr) {
  # Get the column names for the t
  detection <- paste("detect_rate.t",fdr,sep="")
  # Add a column for detection yes/no
  df$detect = FALSE
  # Add which tests can detect
  if(!is.na(tpr)) {df[which(df[,detection] >= tpr),]$detect <- TRUE}
  # Return df
  return(df)
}

#3D Get lowest tf
get.lowest.tf <- function(df,cancer) {
  # Add column for unique combination technique and patient
  df$combination = paste(df$technique,df$sample,sep="-")
  # Generate df for lowest_tf output
  lowest_tf <- data.frame()
  # For each unique combination, output the lowest TF
  for(unique in 1:length(unique(df$combination))) {
    tempdf <- df[which(df$combination == unique(df$combination)[unique]),]
    if(nrow(tempdf[which(tempdf$detect == TRUE),]) > 0) {
      lowest_tf <- rbind(lowest_tf,
                         data.frame(technique = unique(tempdf$technique),
                                    tf = min(tempdf[which(tempdf$detect == TRUE),]$tumor_fraction),
                                    cancer = cancer,
                                    sample = unique(tempdf$sample),
                                    SNV = unique(tempdf$SNV_site)))
    }
    if(nrow(tempdf[which(tempdf$detect == TRUE),]) == 0) {
      lowest_tf <- rbind(lowest_tf,
                         data.frame(technique = unique(tempdf$technique),
                                    tf = 1,
                                    cancer = cancer,
                                    sample = unique(tempdf$sample),
                                    SNV = unique(tempdf$SNV_site)))
    }
    remove(tempdf)
  }
  remove(unique)
  # Return lowest_tf dataframe
  return(lowest_tf)
}

#4 Colors for plotting
colors_tumortype <- c("#EC9D40","#DD6FA3")
colors_techs <- c("#84b8fa","#0070FF","#d89ded","#c855f2")

#5 Thresholds
fdr = 0.68
tpr = 0.95




#### LOAD AND PREPARE DATA ####

#1 Load data
files <- list.files(indir)
files <- files[grep("gz.csv",files)]
dfs <- list()
for(n in 1:length(files)) {
  tempdf <- read.csv(paste(indir,files[n],sep=""), header = TRUE)
  if(colnames(tempdf)[1] == "X") {tempdf <- tempdf[,-1]}
  cancertype = head(unlist(strsplit(files[n],split = ".cs")),n=1)
  cancertype = head(unlist(strsplit(cancertype,split = "_")),n=1)
  dfs <- c(dfs,list(tempdf))
  if(length(dfs) == 1) {names(dfs) <- c(names(df),cancertype)}
  if(length(dfs) > 1) {names(dfs) <- c(names(dfs[-n]),cancertype)}
  remove(tempdf,cancertype)
}
remove(n, files)

#2 Add technique and select subset of the df
for(n in 1:length(dfs)) {
  tempdf <- dfs[[n]]
  tempdf <- add.technique.todf(df = tempdf)
  tempdf <- subset.df(df = tempdf)
  dfs[[n]] <- tempdf
  remove(tempdf)
}
remove(n)




##### PLOT LOWEST TF FOR SINGLE COMBINATION OF THRESHOLDS ####

#1 For all dfs at once: get lowest TF
dfs_lowesttf <- list()
for(n in 1:length(dfs)) {
  tempdf <- dfs[[n]]
  #1A Add detection
  tempdf <- add.detection(df = tempdf,
                          fdr = fdr,
                          tpr = tpr)
  #1B Get lowest tf
  tempdf_lowest <- get.lowest.tf(df = tempdf, cancer = names(dfs)[n])
  dfs_lowesttf <- c(dfs_lowesttf,list(tempdf_lowest))
  if(length(names(dfs_lowesttf)) > 0) {names(dfs_lowesttf) <- c(names(dfs_lowesttf)[-n],names(dfs)[n])}
  if(length(names(dfs_lowesttf)) == 0) {names(dfs_lowesttf) <- names(dfs)[n]}
  #1C Close
  remove(tempdf,tempdf_lowest)
}
remove(n)

#2 Combine the dfs with the lowest tumor fractions per patient for plotting
lowest_df <- data.frame()
for(n in 1:length(dfs_lowesttf)) {
  tempdf <- dfs_lowesttf[[n]]
  lowest_df <- rbind(lowest_df,tempdf)
  remove(tempdf)
}
remove(n)

#3 Plot (this plot was not used in the manuscript; but a similar plot was generated with python)
ggplot(lowest_df,
       aes(y=factor(`technique`, levels = c("NanoRCS Promethion","Promethion","NanoRCS Minion","Minion")),x=tf, fill = cancer, color = cancer)) +
  geom_boxplot(linewidth = 0.2,outlier.size=0.2, color = "gray24") +
  theme_bw() +
  xlab("Lowest detectable TF")+
  ggtitle("Lowest TF detectable") +
  ylab("")+
  scale_x_log10() +
  coord_cartesian(xlim=c(0.001,1.0))+
  scale_fill_manual(values = colors_tumortype)+
  theme(axis.text.x = element_text( vjust = 0.5, hjust=1))

#4 Save the data (for plotting the plot in python)
write.table(lowest_df,paste(indir,"Fig5E.txt",sep=""),sep="\t",row.names = FALSE)




##### PLOT FRACTION OF SAMPLES THAT CAN BE DETECTED BELOW TF 0.01 (1%) ####

#1 Tf threshold
tf_threshold <- 0.02

#2 Get TFs until 1%
tfs <- unique(dfs[[1]]$tumor_fraction)
tfs <- sort(tfs)
tfs <- tfs[-which(tfs > tf_threshold)]

#3 Create a dataframe: for each technique at each tf in each cancer type: How many of the cases do you detect
newdf <- data.frame()
for(tf.detect in tfs) {
  tf_df <- data.frame()
  for(patient in unique(lowest_df$sample)) {
    tempdf <- lowest_df[which(lowest_df$sample == patient),]
    tempdf$tf.detect <- tempdf$tf <= tf.detect
    tf_df <- rbind(tf_df,
                   data.frame(
                     tf = tf.detect,
                     technique = tempdf$technique,
                     detect = tempdf$tf.detect,
                     cancer = tempdf$cancer,
                     snv = tempdf$SNV
                   ))
    remove(tempdf)
  }
  remove(patient)
  for(technique in unique(tf_df$technique)) {
    for(ctype in unique(tf_df$cancer)) {
      temp <- tf_df[which(tf_df$technique == technique & tf_df$cancer == ctype),]
      newdf <- rbind(newdf,
                       data.frame(tf = tf.detect,
                                  technique = technique,
                                  detect = nrow(temp[which(temp$detect == TRUE),])/nrow(temp),
                                  cancer = ctype))
      remove(temp)
    }
    remove(ctype)
  }
  remove(tf_df,technique)
}
remove(tf.detect)

newdf$technique <- factor(newdf$technique, levels = c("NanoRCS Minion","NanoRCS Promethion","Minion","Promethion"))

#4 Plot for manuscript and save plot
plot_manuscript_detection <- ggplot(data = newdf,
       aes(x = tf, y = detect, colour = technique)) +
  geom_line(aes(lty = cancer))+
  #geom_vline(xintercept = 0.0035, lty = 2)+
  theme_bw() +
  xlab("Tumor fraction")+
  ylab("Fraction detection")+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_colour_manual(values = colors_techs) +
  theme(legend.position = "right")
ggsave(paste(outdir,"SupplFig9D.pdf",sep=""),plot_manuscript_detection,width = 7,height = 2)





##### Get some numbers #####

eso <- dfs_lowesttf[[1]]
eso <- eso[which(eso$technique == "NanoRCS Promethion"),]
eso <- eso$tf
median(eso)
mean(eso)
min(eso)
min(eso)
max(eso)

ova <- dfs_lowesttf[[2]]
ova <- ova[which(ova$technique == "NanoRCS Promethion"),]
ova <- ova$tf
median(ova)
mean(ova)
min(ova)
max(ova)
