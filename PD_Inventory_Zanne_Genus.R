################################################################################

# Phylogenetic Tree from Zanne et al. 2014 #

# full reference: 

# Zanne, A. E., Tank, D. C., Cornwell, W. K., Eastman, J. M., Smith, S. A., 
# FitzJohn, R. G. (). Three keys to the radiation of angiosperms into freezing 
# environments. Nature, 506(7486), 89–92. doi:10.1038/nature12872

################################################################################

# set your working directory (the trees and Rsession will be saved here)
# setwd()

library("picante")
library("reshape2")
library("data.table")
library("vegan")
library("plyr")
library("gridExtra")


# first time #
######################################################################

# download Tree WARNING: 75.8 MB will be downloaded!

#temp <- tempfile()
#download.file("http://datadryad.org/bitstream/handle/10255/dryad.59003/PhylogeneticResources.zip?sequence=1",temp)
#Tree <- read.tree(unz(temp, "PhylogeneticResources/Vascular_Plants_rooted.dated.tre"))
#unlink(temp)

# save workspace 
#save.image("Tree.RData")

# subsequent times, load Tree from image
######################################################################


# load workspace with Zannes Tree
load("Tree.RData")

# read in the inventory Species list from Sophia
SPInv <- read.table("FunDivTrees.txt",header=T,sep="\t")

# check which species are not in Phylogeny
SPInv[which (! SPInv$PD.code %in% Tree$tip.label),]$PD.code[
  grep("_",SPInv[which (! SPInv$PD.code %in% Tree$tip.label),]$PD.code)]

# 41/240 species or subspecies not present in in Zannes Tree

# check which genus are not in Phylogeny
SPInv[which (! sub("(\\w)_.+","\\1",SPInv$PD.code) %in% 
                 sub("(\\w)_.+","\\1",Tree$tip.label)),]$PD.code

# Heberdenia 
# Visnea

#subset Tree to include only genuses present in Inventory dataset
#and create tree with only one species per genus (first species per genus kept)

# vector of Genuses in inventory
INV.G <- unique(sub("(\\w)_.+","\\1",SPInv$PD.code))

# drop all species from genus not present in inventory dataset
Drop <- Tree$tip.label[-c(which(sub("(\\w)_.+","\\1",Tree$tip.label) %in% INV.G))]
TreeInv<-drop.tip(Tree,Drop)

# keep only the first species per genus
Drop2<- TreeInv$tip.label[-c(match(unique(sub("(\\w)_.+","\\1",TreeInv$tip.label))
              , sub("(\\w)_.+","\\1",TreeInv$tip.label)))]

TreeInvG <- drop.tip(TreeInv,Drop2)

# cut tip labels to genus only
TreeInvG$tip.label <- sub("(\\w)_.+","\\1",TreeInvG$tip.label)

###############################################################################

###################### calculate PD at Genus level ############################

# load inventory species matrix ("species.matrix", structure is dataframe)
# contains all plots in Inventories with unique plotcode (rows) and species /
# genus identifier as columns. Species identifier is of the form GGGSpSpSp
# colums are filled with basal area value of tree / Genus in specific plot

# load species matrix without french data
load("upscale_inventory_species_matrix_ba_ha.RData")
SPM1 <- species.matrix

# load species matrix with french data
load("upscale_inventory_species_matrix_ba_ha-1.RData")
SPM2 <- species.matrix

# load species matrix with exploratory platform but w/o french data
load("upscale_inventory_extrapolation_species_matrix_ba_ha.RData")
SPM3 <- species.matrix

#check if french plots are in SPM3 (are plots that are SPM2 but not in SPM1 
# in SPM3?)

# subset with french data
French <- SPM2[which(!SPM2$plotcode %in% SPM1$plotcode),]

TRUE %in% (French$plotcode %in% SPM3$plotcode)

# add french plotcodes to full matrix wit exploratories
species.matrix <- rbind.fill(SPM3,French)
species.matrix[is.na(species.matrix)] <- 0

           
#######################
#group species to Genus level (add basal are in each plot of all species of the same Genus) #
#######################


# I work with data.table as "melt" from the reshape2 package and "ddply" from the plyr package take very long on 
# dataframe

# transform first to matrix as somehow I can't get rid of the cast_df appendicies of "species.matrix"
specM <- as.matrix(species.matrix[,-1])
dimnames(specM)<-list(species.matrix$plotcode,colnames(species.matrix)[-1])
specDT <- as.data.table(specM)

# adding plotcode to data.table 
specDT[,plotcode := species.matrix$plotcode]

#reshape dataframe in long format
specDTm <- melt(specDT, id="plotcode", variable.name="speccode", value.name="basal_area")

#add Genus column (first three lettres from speccode)
specDTm[, genus := substr(speccode,1,3)]

# sum all basal areas for all species in one genus in one plot
setkeyv(specDTm,c("plotcode","genus"))
spec_Genus <- specDTm[, sum(basal_area), by=list(plotcode,genus)]

#cast datatabel into wide format
spec_Genus_w <- dcast.data.table(spec_Genus, plotcode~genus) ######## data.table with the basal areas for each Genus in each plot

#transform to species matrix
species.matrix.g <- as.matrix(spec_Genus_w[,-1, with=F])
rownames(species.matrix.g) <- spec_Genus_w[[1]]

##################################################
# calculate Phylogentic diversity 
##################################################

# Transform Genus code in TreeInvG to match genus code in species.matrix.g
TreeInvG$tip.label <- toupper(substr(TreeInvG$tip.label,1,3))

# match species in Tree with species present in species.matrix
colnames(species.matrix.g) %in% toupper(substr(TreeInvG$tip.label,1,3))

#### OBS: all Genus in the species.matrix file are present in the Tree and the file contains fewer Genuses than the list above!

TreeInvGp <- prune.sample(species.matrix.g, TreeInvG)

#calculate psv & pse  [Helmus et al 2007; Phylogenetic measures of biodiversity.]
#calculate pd; [Faith D.P. (1992) Conservation evaluation and phylogenetic diversity.]
#calculate mpd (abundance weighted and binary) [Webb et al 2002; Phylogenies and community ecology.]

#### psv: phylogenetic species variability ####

# quantifies how phylogenetic relatedness decreases 
# the variance of a hypothetical unselected/neutral 
# trait shared by all species in a community. 
# bound between 0 and 1 (1 = max variability; all sp. unrelated)

#### pse: phylogenetic species eveness ####

# modified psv that incorporates abundance information
# bound between 0 and 1 (1 = copletely even community with star phylogeny)

psv.result<-psv(species.matrix.g, TreeInvGp)
psv.result$plotcode <- rownames(psv.result)

pse.result<-pse(species.matrix.g, TreeInvGp) # basal are is taken as species abundance
pse.result$plotcode <- rownames(pse.result)


# pd seems to fail if sites with richness 1 are included wherfore I exclude them manually
## OBS pd takes ~2 minutes  to run! ###
pd.result<-pd(species.matrix.g[which(specnumber(species.matrix.g) > 1),], TreeInvGp)
pd.result$plotcode <- rownames(pd.result)

mpd.result <- mpd(species.matrix.g, cophenetic(TreeInvGp), abundance.weighted = F)
mpd.result_aw <- mpd(species.matrix.g, cophenetic(TreeInvGp), abundance.weighted = T)

MPD <- data.frame(plotcode = rownames(species.matrix.g), mpd = mpd.result, 
                  mpd.aw = mpd.result_aw)

# join the different metrics to common dataframe
PhylDiv <- join(psv.result,pse.result,by="plotcode")
PhylDiv <- join(PhylDiv,pd.result,by="plotcode")
PhylDiv <- join(PhylDiv,MPD,by="plotcode")

# reorder columns and exclude redundant ones 
PhylDiv <- PhylDiv[,c(4,1,3,5,7,9,10,2)]
colnames(PhylDiv) <- c("plotcode", "PSV","PSV_var","PSE","PSR","PSR_var","PD","genus_richness" )

#code NAs as 0 (NAs are produced if the genus richness in a plot is 1, in this case phylogenetic diveristy is 0)
PhylDiv[,2:7][is.na(PhylDiv[,2:7])] <- 0

# sort columns and drop redundant SR
PhylDiv <- PhylDiv[,c(4,1,5,7,9,10,2)]

# export results 
write.table(PhylDiv,"PhylDiv_Inv_Genus.txt",sep="\t")


###########################
# compare the correspondance of species level PD metrics with genus level PD metrics
###########################

# which species from the species.matrix are present in the tree from Zanne et al
# compare the species codes (colnames) with the Species in the SPInv file and compare these species names 
# with the tip.lables of Zannes tree.

TreeSpec <- Tree

# Species list of inventories (excluding Genus)
SPInvSP <- SPInv[nchar(as.character(SPInv$Code)) == 6,]

# further subsetted Species list of inventories with only species also present in species.matrix
SPInvSPmat <- SPInvSP[SPInvSP$Code %in% colnames(species.matrix),]

# which species from the species matrix are in TreeSpec # gives row ind from SPInv
SPm <- which(SPInvSPmat$PD.code %in% TreeSpec$tip.label)

# get corresponding columns in species matrix 
CNp <- which(colnames(species.matrix) %in% SPInvSPmat[SPm,]$Code)

# subset species.matrix for species present in Zannes Tree
SPsub <- species.matrix[,c(1,CNp)]

# exclude plots with 0 abundance
SPsub <- SPsub[-c(which(rowSums(SPsub[,-1]) == 0)),]

# make species matrix
SPsubm <- as.matrix(SPsub[,-1])
dimnames(SPsubm)<-list(SPsub$plotcode,colnames(SPsub[,-1]))

#### make Tree with only species present in species matrix

# DropT <- tips to drop
DropT <- TreeSpec$tip.label[which(! TreeSpec$tip.label %in% SPInvSPmat[SPInvSPmat$Code %in% colnames(SPsubm),]$PD.code)]

#drop tips 
TreeSpec <- drop.tip(TreeSpec, DropT)

# transform tiplabel to same code as used in species matrix (GGGSSS)
TreeSpec$tip.label <- toupper(sub("(\\w\\w\\w)\\w+_(\\w\\w\\w)\\w+","\\1\\2",TreeSpec$tip.label))

#check
identical(sort(TreeSpec$tip.label) ,sort(colnames(SPsubm)))

##################
# make corresponding Species.Matrix and Tree collapsed at Genus level
#################

GenDT <- as.data.table(SPsubm)

# adding plotcode to data.table 
GenDT[,plotcode := rownames(SPsubm)]

#reshape dataframe in long format
GenDTl <- melt(GenDT, id="plotcode", variable.name="speccode", value.name="basal_area")

#add Genus column (first three lettres from speccode)
GenDTl[, genus := substr(speccode,1,3)]

# sum all basal areas for all species in one genus in one plot
setkeyv(GenDTl,c("plotcode","genus"))
GenDTlsub <- GenDTl[, sum(basal_area), by=list(plotcode,genus)]

#cast datatabel into wide format
GenDTw <- dcast.data.table(GenDTlsub, plotcode~genus) ######## data.table with the basal areas for each Genus in each plot

#transform to genus matrix
GenDTm <- as.matrix(GenDTw[,-1, with=F])
rownames(GenDTm) <- GenDTw[[1]]

######## subset make corresponding Genus Tree #########

# keep the first species of each Genus and then replace the species name by genus name only 
TreeGen <- drop.tip(TreeSpec, TreeSpec$tip.label[-c(match(unique(substr(TreeSpec$tip.label,1,3)),substr(TreeSpec$tip.label,1,3)))])
TreeGen$tip.label <- substr(TreeGen$tip.label,1,3)

#check
identical(sort(TreeGen$tip.label),sort(colnames(GenDTm)))

################ calcultae Phylogenetic diveristy for the species-based matrix  ############

##### for species level ######

psv.species<-psv(SPsubm, TreeSpec)
psv.species$plotcode <- rownames(psv.species)

pse.species<-pse(SPsubm, TreeSpec) # basal are is taken as species abundance
pse.species$plotcode <- rownames(pse.species)

psr.species<-psr(SPsubm, TreeSpec)
psr.species$plotcode <- rownames(psr.species)

# pd seems to fail if sites with richness 1 are included wherfore I exclude them manually
## OBS pd takes ~2 minutes  to run! ###
pd.species<-pd(SPsubm[which(specnumber(SPsubm) > 1),], TreeSpec)
pd.species$plotcode <- rownames(pd.species)

# join the different metrics to common dataframe
Phyl_species<-join(psv.species,pse.species,by="plotcode")
Phyl_species<-join(Phyl_species,psr.species,by="plotcode")
Phyl_species<-join(Phyl_species,pd.species,by="plotcode")

# reorder columns and exclude redundant ones 
Phyl_species <- Phyl_species[,c(4,1,3,5,7,9,10,2)]
colnames(Phyl_species) <- c("plotcode", "PSV_S","PSV_var_S","PSE_S","PSR_S","PSR_var_S","PD_S","species_richness" )

#code NAs as 0 (NAs are produced if the genus richness in a plot is 1, in this case phylogenetic diveristy is 0)
Phyl_species[,2:7][is.na(Phyl_species[,2:7])] <- 0

################ calcultae Phylogenetic diveristy for the Genus based matrix ############

##### for Genus level ######

psv.genus<-psv(GenDTm, TreeGen)
psv.genus$plotcode <- rownames(psv.genus)

pse.genus<-pse(GenDTm, TreeGen) # basal are is taken as genus abundance
pse.genus$plotcode <- rownames(pse.genus)

psr.genus<-psr(GenDTm, TreeGen)
psr.genus$plotcode <- rownames(psr.genus)

# pd seems to fail if sites with richness 1 are included wherfore I exclude them manually
## OBS pd takes ~2 minutes  to run! ###
pd.genus<-pd(GenDTm[which(specnumber(GenDTm) > 1),], TreeGen)
pd.genus$plotcode <- rownames(pd.genus)

# join the different metrics to common dataframe
Phyl_genus<-join(psv.genus,pse.genus,by="plotcode")
Phyl_genus<-join(Phyl_genus,psr.genus,by="plotcode")
Phyl_genus<-join(Phyl_genus,pd.genus,by="plotcode")

# reorder columns and exclude redundant ones 
Phyl_genus <- Phyl_genus[,c(4,1,3,5,7,9,10,2)]
colnames(Phyl_genus) <- c("plotcode", "PSV_G","PSV_var_G","PSE_G","PSR_G","PSR_var_G","PD_G","genus_richness" )

#code NAs as 0 (NAs are produced if the genus richness in a plot is 1, in this case phylogenetic diveristy is 0)
Phyl_genus[,2:7][is.na(Phyl_genus[,2:7])] <- 0


########### join and compare Phyl_genus & Phyl_species ############

#join
Phyl_comp <- join(Phyl_species, Phyl_genus)

# PSV
PSV <- ggplot(Phyl_comp, aes(x=PSV_S, y=PSV_G))+
  geom_point(alpha=0.4)+
  theme_bw()

# PSE
PSE <- ggplot(Phyl_comp, aes(x=PSE_S, y=PSE_G))+
  geom_point(alpha=0.4)+
  theme_bw()

# PSR
PSR <- ggplot(Phyl_comp, aes(x=PSR_S, y=PSR_G))+
  geom_point(alpha=0.4)+
  theme_bw()

# PD
PD <- ggplot(Phyl_comp, aes(x=PD_S, y=PD_G))+
  geom_point(alpha=0.4)+
  theme_bw()

# S
S <-  ggplot(Phyl_comp, aes(x=species_richness, y=genus_richness))+
  geom_point(alpha=0.4)+
  theme_bw()


grid.arrange(PSV,PSE,PD,PSR,S,nrow=2)

ggplot(Phyl_comp, aes(y=PD_G, x=species_richness))+
  geom_point(alpha=0.4)+
  theme_bw()


