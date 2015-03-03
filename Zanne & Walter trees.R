################################################################################

# Phylogenetic Tree from Zanne et al. 2014 #

# full reference: 

# Zanne, A. E., Tank, D. C., Cornwell, W. K., Eastman, J. M., Smith, S. A., 
# FitzJohn, R. G. (). Three keys to the radiation of angiosperms into freezing 
# environments. Nature, 506(7486), 89–92. doi:10.1038/nature12872

################################################################################

# set your working directory (the trees and Rsession will be saved here)
setwd("~/Documents/01_PhD/01_Research/04_iDIV workshop")

library("ape")
library("ggplot2")
library("picante")

# first time #
######################################################################

# download Tree WARNING: 75.8 MB will be downloaded!

temp <- tempfile()
download.file("http://datadryad.org/bitstream/handle/10255/dryad.59003/PhylogeneticResources.zip?sequence=1",temp)
Tree <- read.tree(unz(temp, "PhylogeneticResources/Vascular_Plants_rooted.dated.tre"))
unlink(temp)

# save workspace 
save.image("Tree.RData")

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

# check which genus are not in Phylogeny
SPInv[which (! sub("(\\w)_.+","\\1",SPInv$PD.code) %in% 
                 sub("(\\w)_.+","\\1",Tree$tip.label)),]$PD.code

# Heberdenia 
# Visnea

#subset Tree to include only genuses present in Inventory dataset

# vector of Genuses in inventory
INV.G <- unique(sub("(\\w)_.+","\\1",SPInv$PD.code))

#Vectors of Tips to drop
Drop <- Tree$tip.label[-c(which(sub("(\\w)_.+","\\1",Tree$tip.label) %in% INV.G))]

# substree with all species in Experiments and Exploratories
TreeInv<-drop.tip(Tree,Drop)

# create tree with only one species per genus (first species per genus kept)

# vector of tips to drop
Drop2<- TreeInv$tip.label[-c(match(unique(sub("(\\w)_.+","\\1",TreeInv$tip.label))
              , sub("(\\w)_.+","\\1",TreeInv$tip.label)))]

# drop all tips per genus but the first one
TreeInvG <- drop.tip(TreeInv,Drop2)

# cut tip labels to genus only
TreeInvG$tip.label <- sub("(\\w)_.+","\\1",TreeInvG$tip.label)

###############################################################################

########## comparision between Genus based PD and Species PD ##################

load("upscale_inventory_species_matrix_ba_ha.RData")




###############################################################################

#######################     Walter Durka's Tree   ##############################


TreeW<- read.tree("DaPhnE_01/DaPhnE_01.tre")

#check which Species are present in the large Tree
TreeW$tip.label[which(! TreeW$tip.label %in% TreeInv$tip.label)] 

#### only about half the Species are present in both trees

######## subset both trees to common Species ########

# Species present in both Trees
IntSect<- intersect(TreeW$tip.label,TreeInv$tip.label)

TreeSubW <-  drop.tip(TreeInv,TreeInv$tip.label[which (! TreeInv$tip.label %in% IntSect)])

TreeWsub <-   drop.tip(TreeW,TreeW$tip.label[which (! TreeW$tip.label %in% IntSect)])

# create comunity data matrix with Species as colums and samples as rows
# number of Samples
N<-200

CommM<-matrix(data=NA,nrow=N,ncol=length(IntSect))

colnames(CommM)<-TreeWsub$tip.label

# presence absence vector with richness S 
# S = 2 (where biggest differences could be suspected)


S<-2

PrAb <- c(rep(0,(length(IntSect)-S)),rep(1,S))

for (i in 1: N) {
  CommM[i,]<-sample(PrAb,length(PrAb))
}


##### exclude Angiosperm - Gymnosperm pairs ####

# list of Order of Gymnosperms

GymOrd<-c("Abies","Acmopyle","Actinostrobus","Afrocarpus","Agathis","Amentotaxus",
          "Araucaria","Athrotaxis","Austrocedrus","Austrotaxus","Bowenia","Callitris",
          "Calocedrus","Cathaya","Cedrus","Cephalotaxus","Ceratozamia","Chamaecyparis",
          "Columbea","Cryptomeria","Cunninghamia","Cupressus","Cycas","Dacrycarpus",
          "Dacrydium","Dioon","Diselma","Encephalartos","Ephedra","Falcatifolium",
          "Fitzroya","Fokienia","Ginkgo","Glyptostrobus","Gnetum","Halocarpus",
          "Juniperus","Keteleeria","Lagarostrobos","Larix","Lepidothamnus",
          "Lepidozamia","Libocedrus","Macrozamia","Manoao","Margbensonia",
          "Metasequoia","Microbiota","Microcachrys","Microcycas","Microstrobos",
          "Nageia","Neocallitropsis","Nothotsuga","Papuacedrus","Parasitaxus",
          "Phyllocladus","Picea","Pilgerodendron","Pinus","Platycladus",
          "Podocarpus","Prumnopitys","Pseudolarix","Pseudotaxus","Pseudotsuga",
          "Retrophyllum","Sabina","Saxegothaea","Sciadopitys","Sequoia",
          "Sequoiadendron","Stangeria","Sundacarpus","Taiwania","Taxodium","Taxus",
          "Tetraclinis","Thuja","Thujopsis","Thuya","Torreya","Tsuga","Welwitschia",
          "Widdringtonia","Wollemia","Xanthocyparis","Zamia")

isAS<-data.frame(SP1=character(),SP2=character())

# check which pairs contain Angiosperms

for (i in 1:nrow(CommM)){
  isASt <- sub("(\\w+)_\\w+","\\1",colnames(CommM)[which(CommM[i,]==1)]) %in% GymOrd 
  isAS<-rbind(isAS,isASt)
}

# exclude Angiosperm - Gymnosperm pairs form Community matrix

CommM<-CommM[-c(which(rowSums(isAS) == 1)),]

# calculate pd
PDWsub<-pd(CommM,TreeWsub)

# equivalent community matrix for TreeSubW

CommM2<-matrix(data=0,nrow=nrow(CommM),ncol=length(IntSect))
colnames(CommM2)<-TreeSubW$tip.label

for (i in 1:nrow(CommM) ){
  CommM2[i,which(TreeSubW$tip.label %in% colnames(CommM)[which(CommM[i,]==1)])]<-1
}


PDSubW<-pd(CommM2,TreeSubW)


PDcomp<-cbind(PDSubW,PDWsub)

colnames(PDcomp)<-c("PD1","SR1","PD2","SR2")

LM<-summary(lm(PD2~PD1,data=PDcomp))
r<-round(LM$adj.r.squared,3)
p<-signif(LM$coefficients[2,4],3)

pdmaxZ<-round(max(PDcomp[,1]))
pdmaxW<-round(max(PDcomp[,3]))


ggplot(PDcomp, aes(x=PD1,y=PD2)) +
  geom_point(alpha=0.4)+
  annotate("text",y=pdmaxW-20,x=pdmaxZ-100,label=paste("adj. r^2",r,sep=" : "))+
  annotate("text",y=pdmaxW-40,x=pdmaxZ-100,label=paste("p value",p,sep=" : "))+
  labs(title="PD calculated for random Species assamblages (Richness = 2 ) \n
       Zanne et al. ~ Walter et al. ", x="PD Zanne et al",
       y="PD Walter et al")+
  stat_smooth(method="lm")+
  theme_bw(base_size=15)


#################################################################################
          # Smith et al #
#################################################################################

TreeS<-read.tree("Smith_tree.tre")

# check which species are not in Phylogeny
TreeInv[which (! TreeInv %in% TreeS$tip.label)] 

# subset both trees to common species also found in the inventory
ComTrees<-intersect(TreeInv,intersect(TreeS$tip.label,Tree$tip.label))


TreeZsub<-  drop.tip(Tree,Tree$tip.label[which (! Tree$tip.label %in% ComTrees)])

TreeSsub <- drop.tip(TreeS,TreeS$tip.label[which (! TreeS$tip.label %in% ComTrees)])

# create comunity data matrix with Species as colums and samples as rows
# number of Samples
N<-200

CommM<-matrix(data=NA,nrow=N,ncol=length(ComTrees))

colnames(CommM)<-TreeZsub$tip.label

# presence absence vector with richness S 
# S = 2 (where biggest differences could be suspected)


S<-2

PrAb <- c(rep(0,(length(ComTrees)-S)),rep(1,S))

for (i in 1: N) {
  CommM[i,]<-sample(PrAb,length(PrAb))
}

# exclude Angio - Gymno 

isAS<-data.frame(SP1=character(),SP2=character())

# check which pairs contain Angiosperms

for (i in 1:nrow(CommM)){
  isASt <- sub("(\\w+)_\\w+","\\1",colnames(CommM)[which(CommM[i,]==1)]) %in% GymOrd 
  isAS<-rbind(isAS,isASt)
}

# exclude Angiosperm - Gymnosperm pairs form Community matrix

CommM<-CommM[-c(which(rowSums(isAS) == 1)),]

# calculate pd
PDZsub<-psv(CommM,TreeZsub)

# equivalent community matrix for TreeSubW

CommM2<-matrix(data=0,nrow=nrow(CommM),ncol=length(ComTrees))
colnames(CommM2)<-TreeSsub$tip.label

for (i in 1:nrow(CommM) ){
  CommM2[i,which(TreeSsub$tip.label %in% colnames(CommM)[which(CommM[i,]==1)])]<-1
}


PDSsub<-psv(CommM2,TreeSsub)


PDcomp<-cbind(PDSsub,PDZsub)

colnames(PDcomp)<-c("PSV1","SR1","var1","PSV2","SR2","var2")

LM<-summary(lm(PSV2~PSV1,data=PDcomp))
r<-round(LM$adj.r.squared,3)
p<-signif(LM$coefficients[2,4],3)

pdmaxZ<-round(max(PDcomp[,1]))
pdmaxW<-round(max(PDcomp[,4]))


ggplot(PDcomp, aes(x=PSV1,y=PSV2)) +
  geom_point()+
  annotate("text",y=pdmaxW-0.2,x=pdmaxZ-0.3,label=paste("adj. r^2",r,sep=" : "))+
  annotate("text",y=pdmaxW-0.3,x=pdmaxZ-0.3,label=paste("p value",p,sep=" : "))+
  labs(title="PD calculated for random Species assamblages (Richness = 2 ) \n
       Zanne et al. ~ Smith et al. ", x="PD Zanne et al",
       y="PD Smith et al")+
  stat_smooth(method="lm")+
  theme_bw(base_size=15)


