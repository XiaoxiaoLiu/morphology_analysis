
#source("/data/0351/351_Informatics/Analysis/analysis.function.r")
library(multtest)
library(sigclust)
source("/data/informatics/changkyul/Ephys/Script_Repository/CK/ClusTree.function.r")

#############################################
# Set input & output Dirs
#############################################
PARAMDIR = "/data/informatics/changkyul/Ephys/Data"
DATAINDIR = "/data/informatics/changkyul/Ephys/Data/kylel.Aug31.2015.cluster_tree.nonparam"
OUTDIR = "/data/informatics/changkyul/Ephys/Data/kylel.Aug31.2015.cluster_tree.nonparam"
if(!file.exists(OUTDIR)) { system(paste("mkdir", OUTDIR, "; chmod 777 -R", OUTDIR)) }
       
#for (i in 1:1000) {
    print("change the following line specifying input file")
    i <- 1
    din <- read.csv(paste(DATAINDIR, "/EphysOnly_700n_e", i, ".csv", sep=""), header=TRUE)
    print(paste("now, the data file is ", DATAINDIR, "/EphysOnly_700n_e", i, ".csv", sep=""))

    print("====================================================")
    print(" Read in feature file")
    print(" column 'specimen_id', and 'name' are expected")
    print("====================================================")
    specimen_id <- as.character(din[,"specimen_id"])
    samplename <- as.character(din[,"name"])
    tmp1 <- get.field(samplename, ";", 1)
    tmp2 <- get.field(samplename, ";", 2)
    #tmp22 <- get.field(tmp2, '-', 2)
    samplekey <- paste(tmp1, tmp2, sep=";")

    print("	Leaving Out non-feature column from the data file")
    idx.exc <- which(colnames(din) %in% c("name", "specimen_id", "clusterID", "X"))
    fn <- colnames(din)[-c(idx.exc)]
    
    Nfeatall <- length(fn)
    featnameall <- fn
    ephysfeatall <- 1:Nfeatall 
    Nephysall <- length(ephysfeatall)
    featcolorall <- c(rep("white",length(ephysfeatall)))

    Xallna <- matrix(as.numeric(as.matrix(din[,featnameall])), ncol=length(featnameall))
    rownames(Xallna) <- samplekey
    colnames(Xallna) <- featnameall


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

    print("	Get the Firing Pattern Tags")
    firing_table <- read.csv(paste(PARAMDIR, "Tag.firing.08142015.csv", sep="/"))
    any_burst <-  as.numeric(firing_table[,"any_burst"])
    any_delay <-  as.numeric(firing_table[,"any_delay"])
    any_pause <-  as.numeric(firing_table[,"any_pause"])
    names(any_burst) <- names(any_delay) <- names(any_pause) <- as.character(firing_table[,"specimen_id"])
    idx.a <- match(names(any_burst), specimen_id)

    if (0) {
       print("	burst, delay, pause are appended")
       Xall_any <- cbind(Xallna, any_burst[idx.a], any_delay[idx.a], any_pause[idx.a])
    } else {
       print("	burst, delay, pause are NOT appended")
       Xall_any <- Xallna
    }


    col.na <- apply(Xallna, 2, function(x) { N <- length(x) ; sum(is.na(x)) / N } )
    col.na.0.25.name <- colnames(Xallna)[col.na > 0.25]

    print("	NA values are replaced by mean value over all cells")
    Xall.imputed <- t(impute.my(t(Xall_any)))

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
    if (length(mask.sameforall)>0) {
        ephysfeat <- ephysfeatall[-mask.sameforall]
    } else {
        ephysfeat <- ephysfeatall
    }
    Nephys <- length(ephysfeat)
    featcolor <- c(rep("white",length(ephysfeat)))
    if (length(mask.sameforall) > 0) {
        Xall <- Xall.imputed[, -mask.sameforall]
    } else {
        Xall <- Xall.imputed
    }

    print("	Zscore it")
    Zall <- t(zscore(t(Xall))) 
    Zall35 <- Zall
    Zall35[Zall< -3.5] <- -3.5
    Zall35[Zall> 3.5] <- 3.5
    
    print("	Shorten the feature names")
    featnameall <- colnames(Xall)
    featnameallshort  <-  gsub("initial_access", "ini_acc", gsub("trough", "tr", gsub("threshold_", "thr_", gsub("_square", "_sq", gsub("upstroke_downstroke", "updown_stroke", featnameall)))))
    colnames(Zall35) <- colnames(Zall) <- colnames(Xall) <- featnameallshort

    save(Xall, Zall35, file=paste(OUTDIR, "/XZall.Rdata", sep=""))
    
    idx.ephysfeat <- match(ephysfeat, featnameall)
    ephysfeatshort <- featnameallshort[idx.ephysfeat]
    
    print("	Exclude Features with more than 25% of NA")   
    Feature_NA <-  gsub("initial_access", "ini_acc", gsub("trough", "tr", gsub("threshold_", "thr_", gsub("_square", "_sq", gsub("upstroke_downstroke", "updown_stroke", col.na.0.25.name)))))
    idx_NA_tmp <- match(Feature_NA, colnames(Xall))
    if (any(!is.na(idx_NA_tmp))) {
        idx_NA <- idx_NA_tmp[!is.na(idx_NA_tmp)]
    } else {
        idx_NA <- c()
    }

    print("	Exclude Highly Correlated Features > 0.95") 
    Feature_Leave_Out <- as.character(read.csv(paste(PARAMDIR, "Features.Cor.0.95.leaveout.csv", sep="/"), header=FALSE)[,1])
    Feature_Leave_Out <-  gsub("initial_access", "ini_acc", gsub("trough", "tr", gsub("threshold_", "thr_", gsub("_square", "_sq", gsub("upstroke_downstroke", "updown_stroke", Feature_Leave_Out)))))
    idx_tmp <- match(Feature_Leave_Out, colnames(Xall))
    if (any(!is.na(idx_tmp))) {
        idx <- idx_tmp[!is.na(idx_tmp)]
    } else {
        idx <- c()
    }
    idx_exc <- unique(c(idx, idx_NA))
    
    if (length(idx_exc) > 0) {
        XSel <- Xall[, -c(idx_exc)]
    } else {
        XSel <- Xall
    }


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
    
    
    XEphys <- XSel 
    
    str.of.interest <- c("Ephys", "Morph", "EphysMorph")

    Node0 <- list()
    BCT <- list()
    for (soi in str.of.interest[c(1)]) {
         RUNstr <- soi
    
         if (soi == "EphysMorph") {
             Xsoi <- XEphysMorph
             OUTDIRsoi = paste(OUTDIR,  soi, sep="/")
         } 
         if (soi == "Ephys") {
             Xsoi <- XEphys
             OUTDIRsoi = paste(OUTDIR,  "/", soi, "Only", pthr, sep="")
         }
         if (soi == "Morph") {
             Xsoi <- XMorph
             OUTDIRsoi = paste(OUTDIR,  "/", soi, "Only", pthr, sep="")
         }
         OUTDIR.soi = paste(OUTDIRsoi, "_e", i, sep="")

         if(!file.exists(OUTDIR.soi)) { system(paste("mkdir", OUTDIR.soi, "; chmod 777 -R", OUTDIR.soi)) }
    
         XXX <- t(zscore(t(Xsoi)))
         X35 <- XXX
         X35[XXX < -3.5] <- -3.5
         X35[XXX > 3.5] <- 3.5
    
         soi.leg <- list()
         soi.leg$pch <- leg.pch
         soi.leg$crecolor <- hybrid.leg.crecolor
         soi.leg$crecolorL <- hybrid.leg.crecolor2
         soi.leg$crecolor2 <- hybrid.leg.crecolor2
         soi.leg$str <- hybrid.leg.str
         soi.leg$strL <-  hybrid.leg.str2
               
         thisStr <- "Grp.SFp0.01Nshuffled1000" 
         Node0$NodeStr <- "" 
         Node0$strin=paste(OUTDIR.soi, "/", soi, ".Node", thisStr, sep="")
         Node0$pch <- mypch
         Node0$crecolor <- crecolor
         Node0$hybrid.crecolor2 <- hybrid.crecolor2
         Node0$hybrid.crecolorL <- hybrid.crecolor2
         Node0$GT_cre <- hybrid.GT
         
         Node0$Xsoi <- Xsoi
         Node0$ZSel35 <- X35
    
         Node0 <- BuildBinaryClusTree(Node0, soi.leg, Nshuffled=1000, flagDEC="LDA", flagGRP="SEL", flagPthr=pthr, flagPlot=FALSE ) 
         BCT[[soi]] <- Node0
         save(BCT, file=paste(OUTDIR.soi, "/BCT.e_", i, ".Rdata", sep=""))
    }
    
    print(paste("	", i, "-th ClusterTree has built", sep=""))
    #system(paste("rm -r", OUTDIR.soi))  
#}
