# Author: M Jager
# Date: December 2023
# Version: 1.2
# R version 4.2.2


#### START UP ####

#1 Directories
output_dir = "/path/to/NanoRCS/output/output_figures/"
#1B Input
indir = "/path/to/NanoRCS/output/processed_output/07_PCAWG_TF_limit_of_detection/combined/"
#1C Output
outdir = paste(output_dir,"Suppl Fig 9/raw/",sep="")

#2 Libraries
library(ggplot2) #Version: ‘3.4.0’
library(gridExtra) #Version: '2.3'

#3 Functions
add.nr <- function(df) {
  df$nr <- NA
  templength = NA
  for(sample_nr in length(unique(df$SAMPLE)):1) {
    sample <- unique(df$SAMPLE)[sample_nr]
    nr <- NA
    if(sample_nr == length(unique(df$SAMPLE))) {templength = as.numeric(nchar(as.character(df_mutnr[which(df_mutnr$sample == sample),'nr'])))}
    if(nchar(as.character(df_mutnr[which(df_mutnr$sample == sample),'nr'])) < templength) {
      nr <- as.character(df_mutnr[which(df_mutnr$sample == sample),'nr'])
      nr <- paste(paste(rep("0",templength- nchar(nr)),collapse = ""),nr,sep="")
    }
    if(nchar(as.character(df_mutnr[which(df_mutnr$sample == sample),'nr'])) == templength) {
      nr <- as.character(df_mutnr[which(df_mutnr$sample == sample),'nr'])
    }
    df[which(df$SAMPLE == sample),]$nr <- nr
    remove(nr,sample)
  }
  remove(sample_nr,templength)
  return(df)
}

#4 Colors for plotting
colors_tumortype <- c("#EC9D40","#DD6FA3")
colors_VAF <- c("goldenrod", "blue")




#### Get the data ####

#1 Get the cancer types
types <- list.files(indir)
types <- types[grep("PCAWG.csv",types)]
types <- unlist(strsplit(types,split = "_"))
types <- types[grep("cancer",types)]

#2 Load the data
# (A) Load the tables in a list
# (B) get nr of mutations per sample per tumortype
vafs <- list()
df_mutnr <- data.frame()
# Loop
for(type in types) {
  # Load table
  df <- read.table(paste(indir,"SupplFig9BC_",type,"_PCAWG.csv",sep=""),header = TRUE)
  # And get nr of mutations per sample 
  samples <- unique(df$SAMPLE)
  type_df <- data.frame()
  for(sample in samples) {
    tempdf <- data.frame(type,sample,nrow(df[which(df$SAMPLE == sample),]))
    colnames(tempdf) <- c("type","sample","nr")
    type_df <- rbind(type_df,tempdf)
    remove(tempdf)
  }
  remove(sample)
  type_df <- type_df[order(type_df$nr),]
  df_mutnr <- rbind(df_mutnr, type_df)
  remove(samples, type_df)
  # Sort the dataframe with the VAFS based on MUT nr (patients with lowest nr first)
  new_df <- data.frame()
  samples <- df_mutnr[which(df_mutnr$type == type),]$sample
  for(sample in samples) {
    new_df <- rbind(new_df,df[which(df$SAMPLE == sample),])
  }
  remove(sample,samples)
  # Add tables to list
  vafs <- c(vafs,list(new_df))
  if(length(names(vafs)) == 0) {names(vafs) = c(names(vafs),type)}
  if(length(names(vafs)) > 0) {names(vafs) = c(names(vafs)[-length(names(vafs))],type)}
  # Close
  remove(df,new_df)
}
remove(type)

#3 Add column with total nr of mutations of the patient in which mutation was found
for(dfnr in 1:length(vafs)) {
  tempdf <- vafs[[dfnr]]
  tempdf <- add.nr(tempdf)
  vafs[[dfnr]] <- tempdf
  remove(tempdf)
}
remove(dfnr)




#### Plot ####

#1 Plot: number of mutations
nrplot <- ggplot(rbind(df_mutnr[which(df_mutnr$type == "Esophagus-cancer"),],df_mutnr[which(df_mutnr$type == "Ovarian-cancer"),]), aes(x=nr, lty = type, fill = type)) +
  geom_density(fill = NA)+
  geom_histogram(aes(y = after_stat(density)), position = "identity",alpha = 0.3, binwidth = 2500)+
  scale_fill_manual(values = colors_tumortype)+
  theme_bw()+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  ylab("Patients (density)")+
  xlab("Nr of mutations") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = c(.65, .8))

#2 Plot: VAFs
plot1<- ggplot(vafs[[types = "Esophagus-cancer"]],aes(x=VAF,color = nr))+
  geom_density(fill = NA, linewidth = 0.05)+
  theme_bw()+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  ylab("Mutations\n(density)")+
  xlab("VAF") +
  ggtitle("Esophagus cancer") +
  scale_color_manual(values = colorRampPalette(colors_VAF)(length(unique(vafs[[types = "Esophagus-cancer"]]$SAMPLE))))+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none")
plot2 <- ggplot(vafs[[types = "Ovarian-cancer"]],aes(x=VAF,color = nr))+
  geom_density(fill = NA, linewidth = 0.05)+
  theme_bw()+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  ylab("Mutations\n(density)")+
  xlab("VAF") +
  ggtitle("Ovarian cancer") +
  scale_color_manual(values = colorRampPalette(colors_VAF)(length(unique(vafs[[types = "Ovarian-cancer"]]$SAMPLE))))+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none")
vafplot <- grid.arrange(plot1, plot2 ,nrow = 2)
remove(plot1,plot2)

#3 Save plot
ggsave(file = paste(outdir,"SupplFig9BC.pdf",sep=""), plot = grid.arrange(nrplot,vafplot,ncol=2,widths = c(2,3)), width = 7, height = 3)

