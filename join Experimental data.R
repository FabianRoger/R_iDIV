# install.packages("gdata")
library(gdata)
library(plyr)
library(reshape2)

#import experiment dataframes

# Finland
FIN <- read.xls("Plot_descriptors_-_tree_data_-_Finland.xls", sheet="Raw data")
colnames(FIN)

# Germany
GER <- read.xls("Plot_descriptors_-_tree_data_-_Germany.xls", sheet="Raw data")
colnames(GER)

# Italy
ITA <- read.xls("Plot_descriptors_-_tree_data_-_Italy.xls", sheet="Raw data")
colnames(ITA)

# Poland
POL <- read.xls("Plot_descriptors_-_tree_data_-_Poland.xls", sheet="Raw data")
colnames(POL)

# Rumania
RUM <- read.xls("Plot_descriptors_-_tree_data_-_Romania.xls", sheet="Raw data")
colnames(RUM)

# Spain
SPA <- (read.xls("Plot_descriptors_-_tree_data_-_Spain.xls", sheet="Raw data"))
colnames(SPA)

# extract needed columns
Keep_columns <- c("PlotID","SpeciesCode","Species","BasalArea")

FIN <- FIN[,Keep_columns]
GER <- GER[,Keep_columns]
ITA <- ITA[,Keep_columns]
POL <- POL[,Keep_columns]
RUM <- RUM[,Keep_columns]
SPA <- SPA[,Keep_columns]

###### join datasets

Experiments <- rbind(FIN,GER,ITA,POL,RUM,SPA)

# sum up basal areas for each plot
Experiments <- ddply(Experiments, .(PlotID, SpeciesCode, Species), summarise, 
                     BasalArea = sum(BasalArea))


# add Genus colum and sum up basal are for each Genus
Experiment_G <- Experiments
Experiment_G$Genus <- substr(Experiment_G$SpeciesCode,1,3)

Experiment_G <- ddply(Experiment_G, .(PlotID,Genus), summarise, 
                      BasalArea = sum(BasalArea))

#cast into wide format
Experiment_G_w <- dcast(Experiment_G, PlotID ~ Genus) 
Experiment_G_w[is.na(Experiment_G_w)] <- 0

# rename PlotID to plotcode
colnames(Experiment_G_w)[1] <- "plotcode"

save(Experiment_G_w,file="joined_Experiments_Genus.RData")
