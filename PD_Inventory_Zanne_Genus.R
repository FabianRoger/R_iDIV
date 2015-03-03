################################################################################

# Phylogenetic Tree from Zanne et al. 2014 #

# full reference: 

# Zanne, A. E., Tank, D. C., Cornwell, W. K., Eastman, J. M., Smith, S. A., 
# FitzJohn, R. G. (). Three keys to the radiation of angiosperms into freezing 
# environments. Nature, 506(7486), 89–92. doi:10.1038/nature12872

################################################################################

# set your working directory (the trees and Rsession will be saved here)
# setwd()

library("ape")
library("ggplot2")
library("picante")
library("reshape2")
library("data.table")

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


# load workspace
load("Tree.RData")

# list of tree species in the experiments (Kreinitz, FORBIO, ORPHEE, Bangor, B-Tree, 
# Satakunta, BioTree) and all the Exploratory Sites (see below for inventory species included)

TreeSp<-c("Abies_alba","Acer_campestre","Acer_platanoides","Acer_pseudoplatanus",
          "Alnus_glutinosa","Betula_pendula","Carpinus_betulus","Castanea_sativa",
          "Fagus_sylvatica","Fraxinus_excelsior","Larix_decidua","Larix_eurolepis",
          "Larix_kaempferi","Larix_sibirica","Ostrya_carpinifolia","Picea_abies",
          "Pinus_nigra","Pinus_pinaster","Pinus_sylvestris","Populus_tremula",
          "Prunus_avium","Pseudotsuga_menziesii","Quercus_cerris","Quercus_faginea",
          "Quercus_ilex","Quercus_petraea","Quercus_pyrenaica","Quercus_robur",
          "Sorbus_aucuparia","Sorbus_torminalis","Tilia_cordata","Ulmus_glabra")


# check which of the Species are in the Tree

TreeSp[which (! TreeSp %in% Tree$tip.label )] 

# The only missing species is  "Larix_eurolepis" which is a hybrid species between
# "Larix_decidua" and "Larix_kaempferi" (both of which are in the tree)
# "Larix_eurolepis" is only present in the FORBIO experiment

#exclude "Larix_eurolepis"

TreeSp<-TreeSp[-c(which (! TreeSp %in% Tree$tip.label ))]

###################################################################################
### extract subtree with only the species in the Experiments and  Exploratories ### 
###################################################################################

# Vector of Tips to drop
Drop <- Tree$tip.label[-c(which(Tree$tip.label %in% TreeSp))]

# substree with all species in Experiments and Exploratories
TreeExp<-drop.tip(Tree,Drop)

# export Subtree
write.tree(TreeExp, "TreeExp.tree")

###############################################################
###  Iventroy Species (taken from FunDivEUROPE_species.csv) ### 
###############################################################

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
# create tree with only one species per genus (first species per genus kept)

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

########## comparision between Genus based PD and Species PD ##################

# load inventory species matrix ("species.matrix", structure is dataframe)
# contains all plots in Inventories with unique plotcode (rows) and species /
# genus identifier as columns. Species identifier is of the form GGGSpSpSp
# colums are filled with basal area value of tree / Genus in specific plot

load("upscale_inventory_species_matrix_ba_ha.RData")


# group species to Genus level (add basal are in each plot of all species of the same Genus) #

# I work with data.table as "melt" from the reshape2 package and "ddply" from the plyr package take very long on 
# dataframe

# transform first to matrix as somehow I can't get rid of the cast_df appendicies of species.matrix 
specM <- as.matrix(species.matrix[,-1])
dimnames(specM)<-list(species.matrix$plotcode,colnames(species.matrix)[-1])
specDT <- as.data.table(specM)

# adding plotcode as data.table 
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

######## calculate Phylogentic diversity #########

# Transform Genus code in TreeInvG to match genus code in species.matrix.g

TreeInvG$tip.label <- toupper(substr(TreeInvG$tip.label,1,3))

# match species in Tree with species present in species.matrix
TreeInvGp <- prune.sample(species.matrix.g, TreeInvG)

psv.result<-psv(species.matrix.g, TreeInvGp)
psv.result$plotcode <- rownames(psv.result)

pse.result<-pse(species.matrix.g, TreeInvGp)
pse.result$plotcode <- rownames(pse.result)

psr.result<-psr(species.matrix.g, TreeInvGp)
psr.result$plotcode <- rownames(psr.result)

# pd seems to fail if sites with richness 1 are included wherfore I exclude them manually

## OBS pd take quite long to run! ###
pd.result<-pd(species.matrix.g[which(specnumber(species.matrix.g) > 1),], TreeInvGp)
pd.result$plotcode <- rownames(pd.result)


PhylDiv<-join(psv.result,pse.result,by="plotcode")
PhylDiv<-join(PhylDiv,psr.result,by="plotcode")
PhylDiv<-join(PhylDiv,pd.result,by="plotcode")

PhylDiv <- PhylDiv[,c(4,1,3,5,7,9,10,2)]
colnames(PhylDiv) <- c("plotcode", "PSV","PSV_var","PSE","PSR","PSR_var","PD","genus_richness" )

#code NAs as 0 (NAs are produced if the genus richness in a plot is 1, in this case phylogenetic diveristy is 0)
PhylDiv[,2:7][is.na(PhylDiv[,2:7])] <- 0

