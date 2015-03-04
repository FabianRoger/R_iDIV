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

########## calculate PD at Genus level ########################################

# load inventory species matrix ("species.matrix", structure is dataframe)
# contains all plots in Inventories with unique plotcode (rows) and species /
# genus identifier as columns. Species identifier is of the form GGGSpSpSp
# colums are filled with basal area value of tree / Genus in specific plot

load("upscale_inventory_species_matrix_ba_ha.RData")

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
#### OBS: all Genus in the species.matrix file are present in the Tree and the file contains fewer Genuses than the list above!

TreeInvGp <- prune.sample(species.matrix.g, TreeInvG)

#calculate psv, psr, pse  [Helmus et al 2007; Phylogenetic measures of biodiversity.]
#calculate pd; [Faith D.P. (1992) Conservation evaluation and phylogenetic diversity.]

#### psv: phylogenetic species variability ####

# quantifies how phylogenetic relatedness decreases 
# the variance of a hypothetical unselected/neutral 
# trait shared by all species in a community. 
# bound between 0 and 1 (1 = max variability; all sp. unrelated)

#### psr: phylogenetic species eveness ####

# modified psv that incorporates abundance information
# bound between 0 and 1 (1 = copletely even community with star phylogeny)

#### psr: phylogenetic species richness ####

# psv x richness
# "scales" the richness (S) of the community by it's psv
# bound between 0 and S
# often strongly correlated to S 

psv.result<-psv(species.matrix.g, TreeInvGp)
psv.result$plotcode <- rownames(psv.result)

pse.result<-pse(species.matrix.g, TreeInvGp) # basal are is taken as species abundance
pse.result$plotcode <- rownames(pse.result)

psr.result<-psr(species.matrix.g, TreeInvGp)
psr.result$plotcode <- rownames(psr.result)

# pd seems to fail if sites with richness 1 are included wherfore I exclude them manually
## OBS pd takes ~2 minutes  to run! ###
pd.result<-pd(species.matrix.g[which(specnumber(species.matrix.g) > 1),], TreeInvGp)
pd.result$plotcode <- rownames(pd.result)

# join the different metrics to common dataframe
PhylDiv<-join(psv.result,pse.result,by="plotcode")
PhylDiv<-join(PhylDiv,psr.result,by="plotcode")
PhylDiv<-join(PhylDiv,pd.result,by="plotcode")

# reorder columns and exclude redundant ones 
PhylDiv <- PhylDiv[,c(4,1,3,5,7,9,10,2)]
colnames(PhylDiv) <- c("plotcode", "PSV","PSV_var","PSE","PSR","PSR_var","PD","genus_richness" )

#code NAs as 0 (NAs are produced if the genus richness in a plot is 1, in this case phylogenetic diveristy is 0)
PhylDiv[,2:7][is.na(PhylDiv[,2:7])] <- 0

# export results 
#write.table(PhylDiv,"PhylDiv_Inv_Genus.txt",sep="\t")


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

identical(sort(TreeSpec$tip.label) ,sort(colnames(SPsubm)))


