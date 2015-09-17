library(multtest)
library(sigclust)
library(ape)
source("/home/xiaoxiaol/work/src/cell-type-analysis/stats_analysis/stepbystep10featsel.r")
source("/home/xiaoxiaol/work/src/cell-type-analysis/stats_analysis/analysis.function.Ephys.r")


#############################################
# Set input & output Dirs
#############################################
PARAMDIR = "/data/informatics/changkyul/Ephys/Data"
DATAINDIR ="~/work/data/lims2/0903_filtered_ephys_qc"
OUTDIR = "~/work/data/lims2/0903_filtered_ephys_qc"

print("change the following line specifying input file")
i <- 1
din <- read.csv(paste(DATAINDIR, "preprocessed/features_with_db_tags.csv", sep="/"), header=TRUE)

print("====================================================")
print(" Read in feature file")
print(" column 'specimen_id', and 'name' are expected")
print("====================================================")
specimen_id <- as.character(din[,"specimen_id"])
samplename <- as.character(din[,"specimen_name"])
tmp1 <- get.field(samplename, ";", 1)
tmp2 <- get.field(samplename, ";", 2)
samplekey <- paste(tmp1, tmp2, sep=";")  # the same as samplename?


print("	Leaving Out non-feature column from the data file")
idx.exc <- which(colnames(din) %in% c("X","specimen_id","specimen_name","dendrite_type","cre_line","layer","swc_file"))
fn <- colnames(din)[-c(idx.exc)]

Nfeatall <- length(fn)
featnameall <- fn
ephysfeatall <- 1:Nfeatall 
Nephysall <- length(ephysfeatall)
featcolorall <- c(rep("white",length(ephysfeatall)))  # all white 


Xallna <- matrix(as.numeric(as.matrix(din[,featnameall])), ncol=length(featnameall))
rownames(Xallna) <- samplekey
colnames(Xallna) <- featnameall

Xall <-Xallna


print("	Get the Cre-line after  fixing the samplename")
creallin <- samplekey
crealltmp <- get.field(creallin,";", 1)
creall <- gsub("Cre-197628.03.02.01", "Cre", gsub("Cre-197628.06.02.01", "Cre", gsub("Cre-197628.06.01.01", "Cre", crealltmp)))
cre <- unique(creall)

creline_table = read.csv(paste(PARAMDIR, "Cre_line_typeC.06092015.csv", sep="/"))
creC <- as.character(creline_table[,"cre_line"])
typeC <- substr(as.character(creline_table[,"type"]),1,3)
pchC <- as.numeric(creline_table[,"pch1"])
colorC <- as.character(creline_table[,"color"])
print(creC)

print("	Get the Spiny Tags")
spiny_table <- read.csv(paste(PARAMDIR, "Tag.spiny.asipny.08142015.csv", sep="/"))
spiny_type <- gsub("dendrite type - ", "", as.character(spiny_table[,"type"]))
names(spiny_type) <- as.character(spiny_table[,"id"])
idx.t <- match(names(spiny_type), specimen_id)
color.tags <- c("black", "grey", "white", "darkgreen")
colorT <- rep("grey", length(specimen_id))
names(colorT) <- specimen_id
colorT[names(spiny_type)] <- color.tags[match(spiny_type, unique(spiny_type))]



print("	Setting up Sample Color & Symbols") 
mypch <- rep(20,length(creall)) 
crecolor <- rainbow(11+1)[match(creall,creC)] 
print("	Creline based Color Set Up is chosen") 
leg.pch <- rep(20,10) 
leg.crecolor <- rainbow(12)[1:11]  
leg.tagcolor <- color.tags 
print("	Setting up Legend Color & Symbols") 
thisGT <- creall 
leg.pch <- rep(20,10) 
leg.crecolor_t <- c(leg.crecolor, leg.tagcolor) 
leg.crecolor2 <- c(leg.crecolor, leg.tagcolor) 
leg.str <- c(creC) #, 'Excitatory', 'Inhibitory')
leg.str_t <- c(creC, unique(spiny_type)) #, 'Excitatory', 'Inhibitory')
leg.str2 <- c(creC, unique(spiny_type)) #, 'Excitatory', 'Inhibitory')

print("	Now, the layer info is NOT available")

hybrid.GT <- thisGT
hybrid.leg.str2 <- c(leg.str_t)
hybrid.leg.crecolor2 <- leg.crecolor_t
hybrid.leg.str <- leg.str
hybrid.leg.crecolor <- leg.crecolor
hybrid.crecolor <- crecolor
hybrid.crecolor2 <- t(matrix(c(crecolor,colorT), ncol=2))
unique.crecolor <- matrix(unique(crecolor),nrow=1) #topo.colors(10)[match(unique(creall),cre3)]


print("	Colorbar for Ephys, Morph Features")
mask.sameforall <- which(apply(Xall.imputed,2,sd)==0)
# if (length(mask.sameforall)>0) {
#   ephysfeat <- ephysfeatall[-mask.sameforall]
# } else {
#   ephysfeat <- ephysfeatall
# }
# Nephys <- length(ephysfeat)
# featcolor <- c(rep("white",length(ephysfeat)))
# if (length(mask.sameforall) > 0) {
#   Xall <- Xall.imputed[, -mask.sameforall]
# } else {
#   Xall <- Xall.imputed
# }
# 
# print("	Zscore it")
# Zall <- t(zscore(t(Xall))) 
# Zall35 <- Zall
# Zall35[Zall< -3.5] <- -3.5
# Zall35[Zall> 3.5] <- 3.5
# 
# 



XSel = Xall
ZSel = t(zscore(t(XSel))) 
ZSel35 = ZSel
ZSel35[ZSel< -3.5] <- -3.5
ZSel35[ZSel> 3.5] <- 3.5

print(paste("	", i, "-th DATA IS SET"))
print("====================================================")
print("====================================================")
print(paste("	Buildng", i, "-th Tree", sep=""))

set.seed(1)
pthr=0.01

################################### ###################################
################################### ###################################
soi.leg <- list()
soi.leg$pch <- leg.pch
soi.leg$crecolor <- leg.crecolor
soi.leg$crecolorL <- leg.crecolor
soi.leg$crecolor2 <- leg.crecolor2
soi.leg$str <- leg.str
soi.leg$str2 <- leg.str2


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
  OUTDIR.soi = paste(OUTDIRsoi, "_e", i, sep="")
  OUTDIR.soi.str <- paste(soi, ".", pthr, "_e", i, sep="") 
  
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
  Node0 <- BuildBinaryClusTree(Node0, soi.leg, Nshuffled=1000, flagDEC="LDA", flagGRP="SEL", flagPthr=pthr, flagPlot=FALSE ) 
  BCT[[soi]] <- Node0
  save(BCT, file=paste(OUTDIR.soi, "/BCT.e_", i, ".Rdata", sep=""))
  
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


