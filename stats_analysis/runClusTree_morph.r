library(multtest)
library(sigclust)
library(ape)
source("/home/xiaoxiaol/work/src/cell-type-analysis/stats_analysis/stepbystep10featsel.r")
source("/home/xiaoxiaol/work/src/cell-type-analysis/stats_analysis/analysis.function.Ephys.r")


#############################################
# Set input & output Dirs
#############################################
PARAMDIR = "/data/informatics/changkyul/Ephys/Data"
DATAINDIR ="/data/mat/xiaoxiaol/data/lims2/0903_filtered_ephys_qc"
OUTDIR = "/data/mat/xiaoxiaol/data/lims2/0903_filtered_ephys_qc"


data_in <- read.csv(paste(DATAINDIR, "preprocessed/features_with_db_tags_fixed.csv", sep="/"), header=TRUE)

print("====================================================")
print(" Read in feature file")
print(" column 'specimen_id', and 'specimen_name' are required")
print("====================================================")
specimen_id <- as.character(data_in[,"specimen_id"])
samplename <- as.character(data_in[,"specimen_name"])
#tmp1 <- get.field(samplename, ";", 1)
#tmp2 <- get.field(samplename, ";", 2)
#samplekey <- paste(tmp1, tmp2, sep=";")  #some specimen name is manually added, does not have complete info, samplekey may contain "NA" for those
samplekey<-samplename


#	Leaving Out non-feature column from the data file
idx.exc <- which(colnames(data_in) %in% c("X","specimen_id","specimen_name","dendrite_type","cre_line","layer","swc_file"))
all_feature_names <- colnames(data_in)[-c(idx.exc)]
#featcolorall <- c(rep("white",length(all_feature_names)))  # all white 
Xall <- matrix(as.numeric(as.matrix(data_in[,all_feature_names])), ncol=length(all_feature_names))
rownames(Xall) <- samplekey
colnames(Xall) <- all_feature_names

#####-------  if NA values present----------------
#print("  NA values are replaced by mean value over all cells")
#Xall.imputed <- t(impute.my(t(Xall)))



#Get the Cre-line after fixing the samplename
cre_lines <- unique(data_in$cre_line)
dendrite_types <- unique(data_in$dendrite_type)


# -----------------------  for creating the legends and color bars
# use the creline tables for consistency with ephy  data
creline_table = read.csv(paste(PARAMDIR, "Cre_line_typeC.06092015.csv", sep="/"))
#creC <- as.character(creline_table[,"cre_line"])
typeC <- substr(as.character(creline_table[,"type"]),1,3)
pchC <- as.numeric(creline_table[,"pch1"])
colorC <- as.character(creline_table[,"color"])
#print(creC)

print("	Get the Spiny Tags")
color.tags<-c("blue","red","green") #aspiny         spiny          sparsely spiny
# color_dendritetype maps specimen id to a color according to its speciemn id
color_dendritetype <- rep("grey", length(specimen_id))
names(color_dendritetype) <- specimen_id
color_dendritetype[data_in$specimen_id] <- color.tags[match(data_in$dendrite_type, dendrite_types)]


print("	Setting up Sample Color & Symbols") 
mypch <- rep(20,length(data_in$cre_line)) 
crecolor <- rainbow(length(cre_lines)+1)[match(data_in$cre_line,cre_lines)] 

print("	Creline based Color Set Up is chosen") 
leg.pch <- rep(20,10) 
leg.crecolor <- rainbow(length(cre_lines)+1)[1:length(cre_lines)]  
leg.tagcolor <- color.tags  # dendrite types

print("	Setting up Legend Color & Symbols") 
leg.pch <- rep(20,10) 
leg.color_t <- c(leg.crecolor, leg.tagcolor)  #creline types, dendrite types
leg.str <- c(cre_lines)
leg.str_t <- c(cre_lines, dendrite_types) 

### to do : have layer info
print("	Now, the layer info is NOT available")

hybrid.GT <- data_in$cre_line
hybrid.leg.str2 <- c(leg.str_t)
hybrid.leg.crecolor2 <- leg.color_t
hybrid.leg.str <- leg.str
hybrid.leg.crecolor <- leg.crecolor
hybrid.crecolor <- crecolor
hybrid.crecolor2 <- t(matrix(c(crecolor,color_dendritetype), ncol=2))
unique.crecolor <- matrix(unique(crecolor),nrow=1) #topo.colors(10)[match(unique(data_in$cre_line),cre3)]
#----- end of setting up legend ------------------



XSel = Xall
ZSel = t(zscore(t(XSel))) 
ZSel35 = ZSel
ZSel35[ZSel< -3.5] <- -3.5
ZSel35[ZSel> 3.5] <- 3.5

print("====================================================")
print("====================================================")



################################### ###################################
# the stopping critria  p-value threshold
pthr=0.01
################################### ###################################

soi.leg <- list()
soi.leg$pch <- leg.pch
soi.leg$crecolor <- leg.crecolor
soi.leg$crecolorL <- leg.crecolor
soi.leg$crecolor2 <- leg.color_t
soi.leg$str <- leg.str
soi.leg$str2 <- leg.str_t


XMorph <- XSel 

str.of.interest <- c( "Morph" )

Node0 <- list()
BCT <- list()
for (soi in str.of.interest[c(1)]) {
  RUNstr <- soi
  
  if (soi == "EphysMorph") {
    Xsoi <- XEphysMorph
  } 
  if (soi == "Ephys") {
    Xsoi <- XEphys
  }
  if (soi == "Morph") {
    Xsoi <- XMorph
  }
  OUTDIRsoi = paste(OUTDIR,  "/", soi, ".", pthr, sep="")
  OUTDIR.soi = paste(OUTDIRsoi, "_e", sep="")
  OUTDIR.soi.str <- paste(soi, ".", pthr, "_e", sep="") 
  
  if(!file.exists(OUTDIR.soi)) { system(paste("mkdir", OUTDIR.soi, "; chmod 777 -R", OUTDIR.soi)) }
  
  XXX <- t(zscore(t(Xsoi)))
  X35 <- XXX
  X35[XXX < -3.5] <- -3.5
  X35[XXX > 3.5] <- 3.5
  
  
  thisStr <- "LDA_GrpSF0.01Nshuffled100.CreMeanLR_0" 
  Node0$NodeStr <- "" 
  Node0$strin=paste(OUTDIR.soi, "/", soi, ".Node", thisStr, sep="")
  Node0$pch <- mypch
  Node0$crecolor <- crecolor
  Node0$hybrid.crecolor2 <- hybrid.crecolor2
  Node0$hybrid.crecolorL <- hybrid.crecolor2
  Node0$GT_cre <- hybrid.GT
  
  Node0$Xsoi <- Xsoi
  Node0$ZSel35 <- X35
  
  print("building Clustering Tree with iterative PCA")
  Node0 <- BuildBinaryClusTree(Node0, soi.leg, Nshuffled=1000, flagDEC="LDA", flagGRP="SEL", flagPthr=pthr, flagPlot=TRUE ) 
  BCT[[soi]] <- Node0
  save(BCT, file=paste(OUTDIR.soi, "/BCT.e", ".Rdata", sep=""))
  
  print("Heatmap with Cluster Specific Genes")
  subfolder <-paste(soi, "Only.Reduced.Regular.LDA.N5.P", pthr, "n1000.Grpby", ".IncPC_diffFOI.Sigclust", sep="")

  
  ## ??? must be due to code copying
 # gather_AssignedID_plotHeatmap_tree (OUTDIR,subfolder , 0.05, Node0, soi.leg, "") 
  gather_AssignedID_plotHeatmap_tree (OUTDIR.soi,"" , 0.05, Node0, soi.leg, "") 
  
  print("Cluster Diagram")
  tree_txt <- get_meanMembership2(Node0, "") #"(L:0.1,(RL:0.2,RR:0.3):0.15);":w
  
  tree <- read.tree(text=paste(tree_txt,";",sep=""))
  pdf(paste(OUTDIR.soi, "Final.Tree.Diagram.pdf", sep="/"))
  plot(tree, dir="d")
  dev.off()
  
}

print("DONE")


