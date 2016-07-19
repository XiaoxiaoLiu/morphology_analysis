options(warn=-1)

library(matrixcalc)

rm_neuron_NAcls <- function (X) {
   idx.NAcls <- which(X[,"cNA"] > sum(X[1,]) / 2)
   idx.NAcls
}

compare_max_k <- function (x,k) {
    maxid <- order(x, decreasing=TRUE)[1]
    if (x[maxid]>=k) {out <- maxid} else {out <- NA}
    out
}

AssignID_CVRF <- function( outdir,  cvrfFN="RFsummary.classification_primary_secondary_2_mc.100.10fold.csv" ) {

    INFILENAME <- paste(outdir, cvrfFN, sep="/")

    cvrfout <- read.csv(INFILENAME, header=TRUE)

    cvrfout.samplename <- as.character(cvrfout[,"cell"])
    cvrfout.label <- as.matrix(cvrfout[,-c(1,2,3,4)])
    rownames(cvrfout.label) <- cvrfout.samplename
    
    cvrfout.label.woNA <- cvrfout.label[-rm_neuron_NAcls(cvrfout.label), ]
    
    cvrfout.maxlabel <- colnames(cvrfout.label.woNA)[apply(cvrfout.label.woNA, 1, function(x) { order(x, decreasing=TRUE)[1] })]
    names(cvrfout.maxlabel) <- rownames(cvrfout.label.woNA)
    return(cvrfout.maxlabel)
}

test.cv.rf <- function(datver="09102015", norm.dat, cl, n.bin=5,padj.th = 0.01, n.markers=10,othervals=c()) {

    require(randomForest)
    require(e1071)
    bins = sample(1:n.bin, length(cl),replace=T)
    bins=unlist(tapply(names(cl), cl, function(x){
      if(length(x) > n.bin){
        tmp=rep_len(1:n.bin, length(x))
      }else{
        tmp = sample(1:n.bin, length(x))
      }
    #  print(c(length(x), length(tmp)))
      setNames(tmp[sample(length(tmp))], x)
    }))

    if (datver=="09102015") {
       tmp1 <- get.field(names(bins), ";", 1)
       tmp12 <- get.field(tmp1, ".", 2)
       tmp2 <- get.field(names(bins), ";", 2)
       names(bins) <- paste(tmp12, tmp2, sep=";")
    } else {
       tmp1 <- get.field(names(bins), ";", 1)
       tmp2 <- get.field(names(bins), ";", 2)
       tmp22 <- get.field(tmp2, "-", 2)
       names(bins) <- paste(tmp1, tmp22, sep=";")
    }
  #names(cl) <- paste(cl,names(cl), sep=".")
  #bins <- bins[names(cl)]
 # print(bins) 
    rf.pred.cl = setNames(rep(NA, length(cl)), names(cl))
    rf.pred.prob = matrix(0, nrow=length(cl), ncol=length(unique(cl)))
    rownames(rf.pred.prob) <- names(cl)
    colnames(rf.pred.prob) <- unique(cl)
    tmp2<-c()
    for(i in 1:n.bin){
      if (n.bin==1) {
        select.cells = names(cl)
      } else {
        select.cells = names(cl)[bins!=i]
      }
      de.df = DE.genes.pw(norm.dat[,select.cells], paste("cl",cl[select.cells],sep=""))
      de.genes = sapply(de.df, function(x){ x = x[order(x$padj),]
                                            up = x$padj < padj.th & x$lfc > 1
                                            down = x$padj < padj.th & x$lfc < -1
                                            c(head(row.names(x)[up],n.markers), head(row.names(x)[down],n.markers))
                                          }, simplify=F)
      #markers = as.numeric(unique(unlist(de.genes)))
      markers = unique(unlist(de.genes))
      if (length(markers)>1) {
          rf.result<-randomForest(as.matrix(t(norm.dat[markers, select.cells])),as.factor(cl[select.cells]),ntree=1000)
      
          tmp = predict(rf.result, as.matrix(t(norm.dat[markers, names(cl)[bins==i]])),type="prob")
          rf.pred.prob[bins==i,colnames(tmp)] <- tmp
          rf.pred.cl[bins ==i] <-  apply(tmp, 1, which.max)
          if (length(othervals)>0) {
            tmp2 = predict(rf.result, as.matrix(t(norm.dat[markers, names(othervals)])),type="prob")
            tmp2 = apply(tmp2,1,which.max)
          }  
      } else {
          rf.pred.cl[bins == i] = NA;
      }
    }
    if ((length(othervals)>0) & (length(tmp2)>0)) {
       names(tmp2)<-names(othervals)
    } else {
       tmp2<-rep(NA,length(othervals))
       names(tmp2)<-names(othervals)
    }
    return(list(binvals=bins,rf=list(class=rf.pred.cl, posterior = rf.pred.prob,classother=tmp2))) 
}

DE.genes.pw <- function(dat,cl){
  require(limma)  
  design=model.matrix(~0+ as.factor(cl))
  colnames(design)=levels(as.factor(cl))
  fit = lmFit(dat , design=design)
  tmp.cl = colnames(design)
  de.df=list()
  for(i in 1:(length(tmp.cl)-1)){
    for(j in (i+1):length(tmp.cl)){
      x=tmp.cl[i]
      y=tmp.cl[j]
      ctr <<- paste(x, "- ", y)
      contrasts.matrix <- makeContrasts(ctr,  levels=design)
      fit2 = contrasts.fit(fit, contrasts.matrix)
      fit2 = eBayes(fit2)
      padj = apply(fit2$p.value, 2, p.adjust)
      lfc = coef(fit2)
      pair=paste(x,y,sep="_")
      de.df[[pair]]=data.frame(padj=padj[,1],pval=fit2$p.value[,1],lfc=lfc[,1])
      row.names(de.df[[pair]])= row.names(dat)
    }
  }
  return(de.df)
}

###Function to run cross-validation on clusters using Random Forest. This code divides
###the overall sample set into groups of 20%, and then fits a random forest to the 80%
###and predicts the membership of the remaining 20%. All comparisons are pairwise among
###clusters (to minimize limma's inclusion of non-discriminatory genes in DE analysis)
###and then non-dominated clusters are assigned to each sample.
RUN_CVRF <- function ( AssignedID_FN, Feature, outdir, flag.plot=TRUE, datver="09102015" ) {
 
    numruns<-100;  ###number of cross-validation runs
    pthresh<-0.05  ###threshold p-value for differentially expressed genes used for classification
    keepmarkers<-10  ###number of differentially expressed genes used for classification among each pair
    outfile <- paste("RFsummary.full_assignment_rf_mc.", numruns, ".10fold.RData", sep="")   ###output Rdata file with all cross-validation info
    outfile2<- paste("RFsummary.classification_primary_secondary_2_mc.", numruns, ".10fold.csv", sep="")  ###output csv with primary and secondary membership
    outfile3<- paste("RFsummary.10run_mc.", numruns, ".10fold.pdf", sep="")  ###output csv with primary and secondary membership
    outfile4<- paste("RFsummary.classification_Any.Same.Assignment.", numruns, ".10fold.csv", sep="")  ###output csv with primary and secondary membership
            
 
    filename.Id <- AssignedID_FN #"/data/informatics/changkyul/Ephys/RUN09102015/Data.Simulation/seed9woRampEphysOnly.PC0.01.11.csv"
    out <- read.csv(filename.Id, header=TRUE)
    Id.name <- gsub("~sim", "",as.character(out[,1]))
    Id.cre <- get.field(as.character(out[,1]), ";",1)
    Id.id <- as.character(out[,2])
    names(Id.id) <- Id.name
    
    XXX <- Feature
    rownames(XXX) <- gsub("-197628.06.01.01;NA", ";197628.06.01.01", gsub("-197628.06.02.01;NA", ";197628.06.02.01", gsub("-197628.03.02.01;NA", ";197628.03.02.01", rownames(XXX))))
    nfeat <- ncol(XXX)
    idx <- match(rownames(XXX), names(Id.id))
             
    XXX.ClusterID <- Id.id[idx]
    unique.ClusterID <- sort(unique(XXX.ClusterID))
    XXX.label <- match(XXX.ClusterID, unique.ClusterID)
    names(XXX.label) <- names(Id.id)[idx]
    unique.label <- unique(XXX.label)
         
    set.seed(0) ###random number seed
    #excludeclusters<-c(0,11,13,2,4,7)  ###clusters to exclude in the cross-validation, first round
    excludeclusters<-c(0) ###clusters to exclude in the cross-validation, second round
    cross.validation=T;  ###T=cross-validation (i.e. removing 20% of cells), F=classification of new data
         
    rpkmcount<-t(XXX)
    allcl<-paste0("c",XXX.label)
    names(allcl) <- names(XXX.label) 
    
    shuffmat<-c();
    for (ii in (unique(allcl))) {
           shuffmat<-cbind(shuffmat,rpkmcount[,names(allcl)[which(allcl==ii)[1]]])
    } 
    rpkmcount2<-cbind(rpkmcount,shuffmat)
    tempvec2<-rep("Out",ncol(shuffmat))
    names(tempvec2)<-paste("Out",1:ncol(shuffmat),sep=";")
    allcl<-c(allcl,tempvec2)
    clusttab<-table(allcl)
    allclust<-setdiff(unique(allcl),paste0("c",excludeclusters))    
    allclust<-allclust[order(allclust)]
    allcells<-names(allcl[allcl!="Out"])
    colnames(rpkmcount2)<-names(allcl)
         
    ####remove 20% of cells, run cross-validation against every pair of clusters###
    full_assignment_rf_all<-list();
    if (cross.validation) {
           nbins<-10;
    } else {
           nbins<-1;
    }
     
    for (ii in 1:numruns) {
              full_assignment_rf_all[[ii]] <-test.cv.rf(datver, rpkmcount2[,],allcl,n.bin=nbins,padj.th=pthresh,n.markers=keepmarkers,othervals=c())
    }
     

    for (ii in 1:numruns) {
        if (ii==1) {
            Ncluster <- length(unique(allcl))
            NN <- length(full_assignment_rf_all[[ii]]$rf$class)
            fulltab <- matrix(0, nrow=NN, ncol=Ncluster)
            rownames(fulltab) <- names(full_assignment_rf_all[[ii]]$rf$class)
            colnames(fulltab) <- sort(unique(allcl))
        }
     
        ii.names <- names(full_assignment_rf_all[[ii]]$rf$class)
        ii.class <- full_assignment_rf_all[[ii]]$rf$class
        for (n in 1:NN) {
            fulltab[ii.names[n], ii.class[n]] <- fulltab[ii.names[n], ii.class[n]]  + 1
        }
    }      
         
    colnames(fulltab) <- paste0("c",unique.ClusterID[as.numeric(gsub("c","",colnames(fulltab)))])
    nc <- length(colnames(fulltab))
    nx <- length(XXX.ClusterID)
    reordered <- c(order(colnames(fulltab)[1:(nc-1)]) , nc)
    fulltab.reordered <- fulltab[1:nx, reordered]
       
         
    save(full_assignment_rf_all,fulltab, fulltab.reordered, file=paste(outdir, outfile, sep="/"))
 
 
    fulltab.reordered[fulltab.reordered > 100] <- 100
    primvec <- apply(fulltab.reordered, 1, function(x) { order(x,decreasing=TRUE)[c(1)] })
    secondvec <- apply(fulltab.reordered, 1, function(x) { order(x,decreasing=TRUE)[c(2)] })
    outdat<-data.frame(cell=allcells,primvec=primvec,secondvec=secondvec,fulltab.reordered)
    outdat<-outdat[outdat$cell %in% colnames(rpkmcount2),]
    write.csv(outdat,file=paste(outdir, outfile2, sep="/"))
    
    creline_table = read.csv("/data/informatics/changkyul/Ephys/Data/Cre_line_typeC.09102015.csv", header=TRUE)
    creC <- as.character(creline_table[,"cre_line"])
    mycre <- get.field(rownames(fulltab.reordered), ";",1)
    crecolor <- rainbow(11+1)[match(mycre,creC)]
    leg.pch <- rep(15,11)
    leg.str <- creC
    leg.crecolor <- rainbow(12)[1:11] 
    
    myheat.color <- rev(heat.colors(101)) ; myheat.color[101] <- "black"
    if (flag.plot) {
        pdf(paste(outdir, outfile3, sep="/"), height=8, width=16)
        heatmap.3(t(fulltab.reordered), ColSideColors=crecolor, Colv=TRUE, Rowv=FALSE, trace="none", keysize=0.8, margins=c(10,10), col=myheat.color, margin=c(10,10), dendrogram="none", main="Cluster Validation ; 100 run of 10-fold cross validationn")
        legend("bottomleft", col=leg.crecolor, pch=rep(15,11), legend=leg.str,cex=0.75)
    }    
    print(" summary table")
    print(" fulltab.reordered  (702 x Ncluster) ")
    idx.cNA <- match("cNA", colnames(fulltab.reordered))
    print("Samples Always Assigned to NA")
    print(rownames(fulltab.reordered)[fulltab.reordered[, idx.cNA]==10])

    fulltab.reordered.woNA <- fulltab.reordered[, -idx.cNA]

    tmpR <- matrix(rep(apply(fulltab.reordered.woNA, 1, sum), ncol(fulltab.reordered.woNA)), nrow=nrow(fulltab.reordered.woNA))
    fulltab.reordered.woNA.R <- round(100 * fulltab.reordered.woNA / tmpR)
    fulltab.reordered.woNA.R[ is.na(fulltab.reordered.woNA.R) ] <- 0
    class.N2R <- t(apply(fulltab.reordered.woNA.R, 2, function(x) { Ngt0 <- sum(x>0); N100 <- sum(x==100) ; c(Ngt0,N100,round(100*N100/Ngt0)) } ))
    colnames(class.N2R) <- c("N.Any", "N.Same", "Same(%)")
    print(class.N2R)
    write.csv(class.N2R, file=paste(outdir, outfile4, sep="/"))

    colnames(fulltab.reordered.woNA.R) <- paste(class.N2R[,3],colnames(fulltab.reordered.woNA.R),sep="%")
    if (flag.plot) {
        heatmap.3(t(fulltab.reordered.woNA.R), ColSideColors=crecolor, Colv=TRUE, Rowv=FALSE, trace="none", keysize=0.8, margins=c(10,10), col=myheat.color, margin=c(10,10), dendrogram="none", main="Cluster Validation ; 100 run of 10-fold cross validation")
        legend("bottomleft", col=leg.crecolor, pch=rep(15,11), legend=leg.str,cex=0.75)
        dev.off()
    }
}




#########################################

SetUp_Feature_ClusTree <- function (dataFN, dataver, OUTDIR0) {

    #############################################
    # Fill in missing values 
    #############################################
    #DATADir = "/data/informatics/changkyul/Ephys/Data"
    #FNin = "ephys_features_20150311.filtered.csv"
    #FNin = "ephys_features_manual_passed_20150416_v2.csv"
    #dataFN <- paste(DATADir, FNin, sep="/")
    
    if(!file.exists(OUTDIR0)) { system(paste("mkdir", OUTDIR0, "; chmod 777 -R", OUTDIR0)) }
    print(dataver) 
    if (dataver=="09102015") {
        FEAT <- fillNA.local(dataFN, dataver)
    
        din.na <- read.csv(FEAT$filename_NA, header=TRUE)
        col.na <- apply(din.na, 2, function(x) { N <- length(x) ; sum(is.na(x)) / N } )
        col.na.4.name <- colnames(din.na)[col.na > 0.25]
        
        din <- read.csv(FEAT$filename_NA_filled, header=TRUE)
        specimen_id <- as.character(din[,"specimen_id"])
        samplename <- as.character(din[,"name"])
        newcre <- as.character(din[,"newcre"])
        creall <- gsub("_POS", "", gsub("_NEG", "", newcre))
        names(creall) <- newcre
        color.posneg <- rep("white", length(newcre))
        color.posneg[grep("NEG", newcre)] <- "black"
        tmp1 <- get.field(samplename, ";", 1)
        tmp2 <- get.field(samplename, ";", 2)
        tmp22 <- get.field(tmp2, '-', 2)
        samplekey <- paste(creall, specimen_id, sep=";")

        idx_passed <- grep("passed_long", colnames(din))
        fn <- colnames(din)[-c(1:3, idx_passed)]

    } 
    if (dataver=="Morph") {
        din <- read.csv(dataFN, header=TRUE)
        samplename <- as.character(din[,"specimen_name"])
        tmp1 <- get.field(samplename, ";", 1)
        creall <- tmp1
        tmp2 <- get.field(samplename, ";", 2)
        tmp22 <- get.field(tmp2, '-', 2)
        samplekey <- paste(creall, tmp22, sep=";")

        newcre <- paste(creall, "_POS", sep="")
        names(creall) <- newcre
        color.posneg <- rep("white", length(newcre))
        color.posneg[grep("NEG", newcre)] <- "black"

        fn <- colnames(din)[-1]
    }
    
    print("Set up The Feature : Leaving Out First 10 column of data file")
    Nfeatall <- length(fn)
    featnameall <- fn[1:Nfeatall]
    ephysfeatall <- Nfeatall 
    Nephysall <- length(ephysfeatall)
    featcolorall <- c(rep("white",length(ephysfeatall)))
    
    Xallna <- matrix(as.numeric(as.matrix(din[,featnameall])), ncol=length(featnameall))
    rownames(Xallna) <- samplekey
    colnames(Xallna) <- featnameall
    
    print("NA values, if any,  are replaced by mean value over all cells")
    Xall.imputed <- t(impute.my(t(Xallna)))
    
    mask.sameforall <- which(apply(Xall.imputed,2,sd)==0)
    ephysfeat <- ephysfeatall[-mask.sameforall]
    Nephys <- length(ephysfeat)
    featcolor <- c(rep("white",length(ephysfeat)))
    if (length(mask.sameforall) > 0) {
        Xall <- Xall.imputed[, -mask.sameforall]
    } else {
        Xall <- Xall.imputed
    }
    rownames(Xall) <- rownames(Xallna)
    
    #Zall <- t(zscore.mask(Xall, din.na)) 
    Zall <- t(zscore(t(Xall))) 
    rownames(Zall) <- rownames(Xall)
    Zall35 <- Zall
    Zall35[Zall< -3.5] <- -3.5
    Zall35[Zall> 3.5] <- 3.5
    
    
    if (dataver=="09102015") {

        featnameall <- colnames(Xall)
        featnameallshort  <-  gsub("initial_access", "ini_acc", gsub("trough", "tr", gsub("threshold_", "thr_", gsub("_square", "_sq", gsub("upstroke_downstroke", "updown_stroke", featnameall)))))
        colnames(Zall35) <- colnames(Zall) <- colnames(Xall) <- featnameallshort
        
        idx.ephysfeat <- match(ephysfeat, featnameall)
        ephysfeatshort <- featnameallshort[idx.ephysfeat]
        
        featnamealle <- colnames(din.na)
        featnameallshorte  <-  gsub("initial_access", "ini_acc", gsub("trough", "tr", gsub("threshold_", "thr_", gsub("_square", "_sq", gsub("upstroke_downstroke", "updown_stroke", featnamealle)))))
        colnames(din.na) <- featnameallshorte
        
        
        Feature_NA <-  gsub("initial_access", "ini_acc", gsub("trough", "tr", gsub("threshold_", "thr_", gsub("_square", "_sq", gsub("upstroke_downstroke", "updown_stroke", col.na.4.name)))))
        
        idx_NA <- as.vector(na.omit(match(Feature_NA, colnames(Xall))))
        
        print("data is set & saved")
        save(Xall, Zall35, specimen_id, file=paste(OUTDIR0, "/XZall.Rdata", sep=""))
        
        if (flagPlot) {
            print("first heatmap")
            fn=paste(OUTDIR0, 'First.Heatmap.All.ephys.only.1.pdf', sep="/")
            pdf(fn, height=10, width=20)
            hs42 <- heatmap.3(t(Zall35), hclustfun=hhclust, trace="none", Colv=TRUE, ColSideColors=hybrid.crecolor2, 
                      Rowv=FALSE, RowSideColors=featcolor,
                      main=paste("Data-driven Clustering, Ephys"), keysize=0.8, 
                      margins=c(10,13), cexCol=0.5,cexRow=1, col=rbg)
            legend("bottomleft", bg="white", legend=leg.str2, pch=15, col=leg.crecolor2, cex=0.75) #, bty="n")
            dev.off()
        
            print("PCA ")
            titlestr='PCA of Cells : zscore, |z|<3.5'
            outfile=paste(OUTDIR0, 'First.PCA.pdf', sep="/")
            pdf(outfile)
            myPCA <- plotPCA(Zall35, titlestr, mypch, leg.pch, crecolor, leg.crecolor2, leg.str2, outfile) 
            pairs(myPCA$x[,1:8], col=crecolor, pch=mypch)
            dev.off()
        
            print("MDS ")
            outfile=paste(OUTDIR0, 'Top.MDS.pdf', sep="/")
            pdf(outfile)
            titlestr='MDS of Cells : zscore, |z|<3.5'
            plotMDS(Zall35, titlestr, mypch, leg.pch, crecolor, leg.crecolor2, leg.str2) 
            dev.off()
        }
        
        print("Feature Reduction")
        Cor.allna <- cor(Xall, use="pairwise.complete.obs")
        Cor.all <- Cor.allna
        Cor.all[is.na(Cor.all)] <- 0
        
        
        Feature_Leave_Out <- as.character(read.csv("/data/informatics/changkyul/Ephys/RUN07012015.Extended/Features.Cor.0.95.leaveout.csv", header=FALSE)[,1])
        
        idx_cor <- match(Feature_Leave_Out, colnames(Xall))
        idx_exc <- unique(c(idx_cor, idx_NA))
        
        XSel <- Xall[, -c(idx_exc)]
    } else {

        XSel <- Xall
        colnames(XSel) <- colnames(Xall) 

    }

    rownames(XSel) <- rownames(Xall)
    ZSel = t(zscore(t(XSel))) 
    rownames(ZSel) <- rownames(XSel)
    #ZSel = t(zscore.mask(XSel, din.na)) 
    ZSel35 = ZSel
    ZSel35[ZSel< -3.5] <- -3.5
    ZSel35[ZSel> 3.5] <- 3.5
    
    #selfeatcolor <- featcolor[-idx]
    save(Xall, XSel, ZSel35, file=paste(OUTDIR0, "/XZall.Rdata", sep=""))

    OUT <- list()
    OUT$Feature <- XSel
    OUT$creall  <- creall    

    return(OUT)
}

SetUp_Plot_ClusTree <- function( creall ) { 

       print("plotting parameters set up")
       DATADir = "/data/informatics/changkyul/Ephys/Data"
       creline_table = read.csv(paste(DATADir, "Cre_line_typeC.09102015.csv", sep="/"))
       creC <- as.character(creline_table[,"cre_line"])
       typeC <- substr(as.character(creline_table[,"type"]),1,3)
       pchC <- as.numeric(creline_table[,"pch1"])
       colorC <- as.character(creline_table[,"color"])
       cre <- creC
        
       print("Setting up Color & Symbols")
       mycolor <- rainbow(12)
       mypch <- rep(20,length(creall))
       crecolor <- mycolor[match(creall,creC)]
       print("Creline based Color Set Up is chosen")
       leg.pch <- rep(20,12)
       leg.crecolor <- mycolor
       leg.str <- c(creC) #, 'Excitatory', 'Inhibitory')
       leg.crecolor2 <- c(leg.crecolor, "white", "black")
       leg.str2 <- c(creC, 'Cre_Positive', 'Cre_Negative')
    
       print("Uncomment when the layer info is available")

       newcre <- names(creall)
       color.posneg <- rep("white", length(newcre))
       color.posneg[grep("NEG", newcre)] <- "black"

       hybrid.GT <- creall
       hybrid.leg.str2 <- c(leg.str)
       hybrid.leg.str <- leg.str
       hybrid.leg.crecolor <- leg.crecolor
       hybrid.crecolor <- crecolor
       hybrid.crecolor2 <- t(matrix(c(crecolor,color.posneg), ncol=2))
       unique.crecolor <- matrix(unique(crecolor),nrow=1) #topo.colors(10)[match(unique(creall),cre3)]
    
       plotparams <- list()
    
       leg <- list()
       leg$pch <- leg.pch
       leg$crecolor <- leg.crecolor
       leg$crecolorL <- leg.crecolor
       leg$crecolor2 <- leg.crecolor2
       leg$str <- leg.str
       leg$str2 <- leg.str2
        
       plotparams$leg <- leg
       plotparams$pch <- mypch
       plotparams$crecolor <- crecolor
       plotparams$crecolorL <- crecolor
       plotparams$hybrid.crecolor <- hybrid.crecolor
       plotparams$hybrid.crecolor2 <- hybrid.crecolor2
       plotparams$hybrid.crecolorL <- hybrid.crecolor2
       plotparams$GT_cre <- hybrid.GT
    
   return(plotparams)
}


fillNA <- function( X ) {
    XXX <- as.numeric(X)
    Xmean <- mean(XXX,na.rm=TRUE)
    OUT <- XXX
    OUT[is.na(XXX) | is.nan(XXX)] <- Xmean
    OUT
}

fill.local.missing <- function( XXX ) {
    X <- as.numeric(XXX)
    t <- 1:6 
    OUT <- X
    Y <- X
    idx.na <- which(is.na(X))
    N.na <- length(idx.na)
    idx.nonna <- which(!is.na(X))
    N.nonna <- length(idx.nonna)
#    print(N.nonna)
    if (N.nonna==6 || N.nonna==0) {OUT <- X}
    if (N.nonna==1) {OUT <- rep(X[idx.nonna],6)}
    if (N.nonna==2) {
     #   print(idx.nonna)
        t1 <- t[idx.nonna[1]] ; x1 <- X[t1]
        t2 <- t[idx.nonna[2]] ; x2 <- X[t2]
        dt <- t2-t1
        dx <- x2-x1
        slope <- dx/dt
        midt <- (t1+t2)/2
        midx <- (x1+x2)/2

        for (i in 1:6) {
            di <- (i-midt) 
            xi <- midx + (slope*di)
            Y[i] <- xi
        }
        OUT[idx.na] <- Y[idx.na]
    } 

    if (N.nonna==3) {
     #   print(idx.nonna)
        t1 <- t[idx.nonna[1]] ; x1 <- X[t1]
        t2 <- t[idx.nonna[2]] ; x2 <- X[t2]
        t3 <- t[idx.nonna[3]] ; x3 <- X[t3]
     #   print(c(t1,t2,t3))
     #   print(c(x1,x2,x3))
        dt <- t3-t1 
        dx <- x3-x1
        slope <- dx/dt
        midt <- (t1+t2+t3)/3
        midx <- (x1+x2+x3)/3

        for (i in 1:6) {
            di <- (i-midt) 
            xi <- midx + (slope*di)
            Y[i] <- xi
        }
        OUT[idx.na] <- Y[idx.na]
    } 
    if (N.nonna==4) {
     #   print(idx.nonna)
        t1 <- t[idx.nonna[1]] ; x1 <- X[t1]
        t2 <- t[idx.nonna[2]] ; x2 <- X[t2]
        t3 <- t[idx.nonna[3]] ; x3 <- X[t3]
        t4 <- t[idx.nonna[4]] ; x4 <- X[t4]
        dt <- t4-t1 
        dx <- x4-x1
        slope <- dx/dt
        midt <- (t1+t2+t3+t4)/4
        midx <- (x1+x2+x3+x4)/4

        for (i in 1:6) {
            di <- (i-midt) 
            xi <- midx + (slope*di)
            Y[i] <- xi
        }
        OUT[idx.na] <- Y[idx.na]
    } 
    if (N.nonna==5) {
 #       print(idx.nonna)
        t1 <- t[idx.nonna[1]] ; x1 <- X[t1]
        t2 <- t[idx.nonna[2]] ; x2 <- X[t2]
        t3 <- t[idx.nonna[3]] ; x3 <- X[t3]
        t4 <- t[idx.nonna[4]] ; x4 <- X[t4]
        t5 <- t[idx.nonna[5]] ; x5 <- X[t5]
        dt <- t5-t1 
        dx <- x5-x1
#print(c(t1,t2,t3,t4,t5))
#print(c(x1,x2,x3,x4,x5))
        slope <- dx/dt
        midt <- (t1+t2+t3+t4+t5)/5
        midx <- (x1+x2+x3+x4+x5)/5
#print(c(dt,dx,slope,midt,midx))

        for (i in 1:6) {
            di <- (i-midt) 
            xi <- midx + (slope*di)
            Y[i] <- xi
#            print(c(i, di, Y[i]))
        }
        OUT[idx.na] <- Y[idx.na]
#       print(c(idx.na, Y[idx.na]))
    } 

#print(OUT)
    OUT
}


fillNA.local <- function (dataFN, dataver="09102015") {

    #DATADir = "/data/informatics/changkyul/Ephys/Data/Ephys.Feature.09102015.csv"
    FeatureFile.Nathan <- dataFN
    
    mytable <- read.csv(FeatureFile.Nathan, header=TRUE)
    idx.rowno <- 1
    idx.name <- match("name", colnames(mytable))
    idx.line <- match("line", colnames(mytable))
    idx.specimen_id <- match("specimen_id", colnames(mytable))
    idx.cre_line <- match("cre_line", colnames(mytable))

    name <- as.character(mytable[,idx.name])
    specimen_id <- as.character(mytable[, idx.specimen_id])

    if (dataver=="09102015") {
        idx.cre_reporter_positive <- match("cre_reporter_positive", colnames(mytable))
        cre_reporter <- rep("POS", nrow(mytable)) 
        cre_reporter[which(as.character(mytable[, idx.cre_reporter_positive])=="FALSE")] <- "NEG"
    
        newcre <- paste(as.character(mytable[, idx.cre_line]), cre_reporter, sep="_")
        idx.nonfeat <- c(idx.rowno, idx.name, idx.line, idx.specimen_id, idx.cre_line, idx.cre_reporter_positive)
    } 
    if (dataver=="07012015") {
        newcre <- as.character(mytable[, idx.cre_line])
        idx.nonfeat <- c(idx.rowno, idx.name, idx.line, idx.specimen_id, idx.cre_line)
    }

    filled_local_na_table <- as.matrix(mytable[,-idx.nonfeat])
    print("sub NA with local mean or mean-4sd")

    idx0 <- grep("long_square_0", colnames(filled_local_na_table))
    featname0 <- gsub("_0", "",colnames(filled_local_na_table)[idx0])
    filled_na_table <- filled_local_na_table
    for (f in 1:length(idx0)) {
         f.idx <- grep(featname0[f], colnames(filled_local_na_table))
         print(paste(f, featname0[f]))
    #     print(colnames(filled_local_na_table)[f.idx])
         f.val <- filled_local_na_table[, f.idx]
         f.local <- t(apply(f.val, 1, fill.local.missing))
         filled_na_table[, f.idx] <- f.local
    }
    
    filled_local_table <- apply(filled_na_table, 2, fillNA)

    filename_NA <- gsub( ".csv", "_NA.csv", FeatureFile.Nathan)
    write.csv(data.frame(name, specimen_id, newcre, filled_na_table), file=filename_NA, row.names=FALSE)
    filename_NA_filled <- gsub(".csv", "_NA_filled.csv", FeatureFile.Nathan)
    write.csv(data.frame(name, specimen_id, newcre, filled_local_table), file=filename_NA_filled,row.names=FALSE)
    OUT <- list()
    OUT$filename_NA <- filename_NA
    OUT$filename_NA_filled <- filename_NA_filled
    
    return(OUT)
}



BuildClusTree <- function ( Xin, params, OUTDIR ) {

    print("Run with R version 3.1.0 (2014-04-10)") 

    library(sigclust)
    library(multtest)
    library( ape )

    seedk     <- params$rndseed
    grp       <- params$grouping
    pthr      <- params$pthr
    partition <- params$partition
    flagZ     <- params$flagZ
    DEC       <- params
    pre       <- params$PRE
    soi       <- "Ephys"
    nsim      <- params$Nshuffled
    sparse    <- params$sparse
    featsel   <- params$featsel

    plotparams <- params$plot
    
    set.seed(seedk)
    if (params$sparse) {hcluststr <- 'Sparse'} else { hcluststr <- 'Regular' }
    thisStr <-  paste(soi, ".seed", seedk, ".", hcluststr, ".P", pthr, ".Grpby", grp, ".FSby", featsel, ".IncPC.", partition, sep="")
                   
    OUTDIR.soi = paste(OUTDIR,  thisStr, sep="/")
    if (!file.exists(OUTDIR.soi)) { system(paste("mkdir", OUTDIR.soi, "; chmod 777 -R", OUTDIR.soi)) }
                       
    if (flagZ) {
        ZZZ <- t(zscore(t(Xin)))
        Z35 <- ZZZ
        Z35[ZZZ < -3.5] <- -3.5
        Z35[ZZZ > 3.5] <- 3.5
    } else {
        Z35 <- Xin
    }
                       
    Node0         <- list()
    Node0$NodeStr <- "Node0" 
    Node0$strin   <- paste(OUTDIR.soi, "/Node", sep="")
    Node0$Xsoi    <- Xin
    Node0$ZSel35  <- Z35

    Node0$pch              <- plotparams$pch 
    Node0$GT_cre           <- plotparams$GT_cre
    Node0$crecolor         <- plotparams$crecolor
    Node0$crecolorL        <- plotparams$crecolor
    Node0$hybrid.crecolor  <- plotparams$hybrid.crecolor
    Node0$hybrid.crecolor2 <- plotparams$hybrid.crecolor2
    Node0$hybrid.crecolorL <- plotparams$hybrid.crecolor2

    leg <- plotparams$leg 

    print(paste("Building Tree using Iterative", params$partition))
                       
    Node0 <- BuildBinaryClusTree(Node0, plotparams, Nshuffled=nsim, flagDEC="LDA", flagGRP=grp, flagPthr=pthr, flagSparse=sparse, flagPlot=FALSE, flagDIP=FALSE, flagMembership=TRUE, flagPartition=partition, flagSEL=featsel, flagPRE=pre)
    
    save(Node0, file=paste(OUTDIR.soi, "Node0.Rdata", sep="/"))
                   
    print("Heatmap with Cluster Specific Genes")
    clusterFN <- gather_AssignedID_plotHeatmap_tree (OUTDIR, thisStr, 0.05, Node0, plotparams$leg, "", flagPlot=TRUE) 
    
    print("DONE")
    return(clusterFN)
    
}


SetUpChild <- function (Nodei, idx, thisStr) {

    NodeOUT <- list()
    NodeOUT$NodeStr <- paste(Nodei$NodeStr, thisStr, sep="")
    NodeOUT$strin <- paste(Nodei$strin, thisStr, sep="")
    NodeOUT$pch <- Nodei$pch[idx]
    NodeOUT$crecolor <- Nodei$crecolor[idx]
    NodeOUT$hybrid.crecolor2 <- Nodei$hybrid.crecolor2[, idx]
    NodeOUT$hybrid.crecolorL <- Nodei$hybrid.crecolorL[, idx]
    NodeOUT$GT_cre <- Nodei$GT_cre[idx]

    if (length(idx)>1) {
        NodeOUT$Xsoi <- Nodei$Xsoi[idx,]
    } else {
        NodeOUT$Xsoi <- t(as.matrix(Nodei$Xsoi[idx,],nrow=1,ncol=ncol(Nodei$Xsoi)))
        rownames(NodeOUT$Xsoi) <- rownames(Nodei$Xsoi)[idx]
        colnames(NodeOUT$Xsoi) <- colnames(Nodei$Xsoi)
    }
    NodeOUT.ZSel <- t(zscore(t(NodeOUT$Xsoi))) 
    NodeOUT$ZSel35 <- NodeOUT.ZSel
    NodeOUT$ZSel35[NodeOUT.ZSel< -3.5] <- -3.5
    NodeOUT$ZSel35[NodeOUT.ZSel> 3.5] <- 3.5

   NodeOUT

}

Assign_NodeID <- function(Nodei) {
    if (nrow(Nodei$Xsoi)==1) { print(paste(rownames(Nodei$Xsoi), Nodei$Assigned[rownames(Nodei$Xsoi)])) }
    Nodei$Assigned[rownames(Nodei$Xsoi)] <- Nodei$NodeStr

    write.csv(data.frame(samplename=rownames(Nodei$Xsoi), clusterID=Nodei$Assigned[rownames(Nodei$Xsoi)]), file=paste(Nodei$strin, ".AssignedID.csv", sep=""))
}


BuildBinaryClusTree <- function(Nodei, plotparams, Nshuffled=1000, flagDEC="LDA", flagGRP="SEL", flagPthr=0.05, flagSparse=FALSE, flagPlot=TRUE, flagDIP=FALSE, flagMembership=TRUE, flagPartition="PCA", flagSEL="NA", flagPRE="sigclust", flagMinN=5, flagDebug=FALSE ) {

    if(flagDebug) {
        Nodei <- Nodei$Right
        Nshuffled=1000
        flagDEC="LDA"
        flagGRP="SEL"
        flagPthr=0.05
        flagSparse=FALSE
        flagPlot=FALSE 
        flagDebug=TRUE 
        flagMembership=TRUE
    }

    soi.leg <- plotparams$leg

    TN.Sample.thr <- flagMinN
    print("##########################################################")
    print(paste("Numbef of samples :", nrow(Nodei$Xsoi),  "&   Number of Features : ", ncol(Nodei$Xsoi)))
    if (!is.matrix(Nodei$Xsoi) || nrow(Nodei$Xsoi) <= TN.Sample.thr) {
        Nodei$idxR <- c()
        Nodei$idxL <- c()
        Nodei$terminal <- TRUE
        Assign_NodeID (Nodei)
        print(paste("  ", Nodei$NodeStr, ": terminal node...STOP # Not enough Sample", nrow(Nodei$Xoi))) 

    } else {
        print(paste("   BuildBinaryClusTree", Nodei$NodeStr))
        Nodei <- clustering2decR(Nodei, Nodei$Xsoi, Nodei$GT_cre, Nodei$ZSel35, Nodei$strin, Nodei$pch, soi.leg$pch, Nodei$crecolor, Nodei$hybrid.crecolorL, 
                            soi.leg$crecolor, soi.leg$crecolorL, soi.leg$str, soi.leg$strL, Nodei$NodeStr, Nshuffled, flagDEC, flagGRP, flagPthr, flagSparse, flagPlot, flagDIP, flagMembership, flagPartition, flagSEL, flagPRE, flagDebug) 
    
        print(paste("   BinaryClusTree is Built with ", Nodei$NodeStr))

        if (Nodei$nPC > 0) {
            Nodei.table <- table(Nodei$grouping)
            if (0) {
                if (Nodei.table[1] > Nodei.table[2]) {
                     Nodei$idxR <- which(Nodei$grouping==1)
                     Nodei$idxL <- which(Nodei$grouping==2)
                } else {
                     Nodei$idxR <- which(Nodei$grouping==2)
                     Nodei$idxL <- which(Nodei$grouping==1)
                }
            } else {
                Nodei$idxR <- which(Nodei$grouping==2)
                Nodei$idxL <- which(Nodei$grouping==1)
                Nodei$idxRPC <- which(Nodei$groupingPC==2)
                Nodei$idxLPC <- which(Nodei$groupingPC==1)
            }
            Nodei$terminal <- FALSE
            save(Nodei, file=paste(Nodei$strin, ".Rdata", sep=""))
            print("##########################################################")
            print("Left") 
            Nodei$Left <- SetUpChild( Nodei, Nodei$idxL, "L" )
            Nodei$Left <- BuildBinaryClusTree ( Nodei$Left, plotparams, Nshuffled, flagDEC, flagGRP, flagPthr, flagSparse, flagPlot, flagDIP, flagMembership, flagPartition, flagSEL, flagPRE, flagMinN, flagDebug )
    
            print("Right") 
            NodeRStr <- paste(Nodei$Str, "R", sep="")
            Nodei$Right <- SetUpChild( Nodei, Nodei$idxR, "R" )
            Nodei$Right <- BuildBinaryClusTree ( Nodei$Right, plotparams, Nshuffled, flagDEC, flagGRP, flagPthr, flagSparse, flagPlot, flagDIP, flagMembership, flagPartition, flagSEL, flagPRE, flagMinN, flagDebug )
    
        } else {
            print(paste(Nodei$NodeStr, ": terminal node...STOP")) 
            Nodei$idxR <- c()
            Nodei$idxL <- c()
            Nodei$terminal <- TRUE
            Assign_NodeID (Nodei)
        }
    
    }

    Nodei
}

apply_Ttest <- function(X, G) {
    tmpXX <- X  +  X*runif(length(X), 0.01, 0.02)
    if (var(tmpXX)==0) {
        pout <- 1
    } else { 
        myt <- t.test(tmpXX ~ G, var.equal=FALSE)
        pout <-  myt$p.value
    }
    pout
}

apply_TtestM12 <- function(X, G) {
    tmpXX <- X  +  X*runif(length(X), 0.01, 0.02)
    if (var(tmpXX)==0) {
        pout <- c(X[1], X[1])
    } else { 
        myt <- t.test(tmpXX ~ G, var.equal=FALSE)
        pout <-  myt$estimate
    }
    pout
}
    
order_feat_byJackstraw <- function(Zin, nPC, Xin) {
    out <- jackstraw(t(Zin), PC=1:nPC, r=nPC, B=100)

    TadjP <- out$p.value
    porder <- order(TadjP, decreasing=FALSE)

    OUT <- list()
    OUT$order <- porder
    OUT$pval.ordered <- TadjP[porder] 
    OUT$pval.lt0.01 <- sum(TadjP[porder]<0.01, na.rm=TRUE) 
    OUT$pval.lt0.05 <- sum(TadjP[porder]<0.05, na.rm=TRUE) 

    OUT
}

order_feat_byJackstraw <- function(Zin, nPC, Xin) {
    out <- jackstraw(t(Zin), PC=1:nPC, r=nPC, B=100)

    TadjP <- out$p.value
    porder <- order(TadjP, decreasing=FALSE)

    OUT <- list()
    OUT$order <- porder
    OUT$pval.ordered <- TadjP[porder] 
    OUT$pval.lt0.01 <- sum(TadjP[porder]<0.01, na.rm=TRUE) 
    OUT$pval.lt0.05 <- sum(TadjP[porder]<0.05, na.rm=TRUE) 

    OUT
}


order_feat_byNA <- function(Zin, Gin, Xin) {

    TadjP <- rep(0.001, ncol(Zin))
    porder <- order(TadjP, decreasing=FALSE)

    OUT <- list()
    OUT$order <- porder
    OUT$pval.ordered <- TadjP[porder] 
    OUT$pval.lt0.01 <- sum(TadjP[porder]<0.01, na.rm=TRUE) 
    OUT$pval.lt0.05 <- sum(TadjP[porder]<0.05, na.rm=TRUE) 

    OUT
}

get_signif_PC <- function(ZSel35, ppp, Nsim=1000) {

if (Nsim==0) {
    print("Using Broken-Stick method")
    screeplot.prcomp(ppp, bstick=TRUE)
    mask.signif <- bstick.prcomp(ppp) < ppp$sdev^2

} else {
    maxvarR <- 0
    print("Using max(PC1(shuffled))")
    ivarR <- rep(0,Nsim)
    for (i in 1:Nsim) {
         ZSel35.shuffledi <- t(apply(ZSel35, 2, function(x) { x[sample(length(x))]}))
         ivarR[i] <- summary(prcomp(ZSel35.shuffledi))$importance[2,1]
    }
    meanvarR <- mean(ivarR)
    sdvarR <- sd(ivarR)
    maxvarR <- max(ivarR)

    mean3.5sdvarR <- meanvarR + (3.5*sdvarR)
    if (maxvarR > mean3.5sdvarR) thrvarR <- maxvarR
    else thrvarR <- mean3.5sdvarR
    print(c(maxvarR, mean3.5sdvarR, thrvarR))
    mask.signif <- summary(prcomp(ZSel35))$importance[2,] > thrvarR 
}

mask.signif
}


binary_partition_PC <- function (Pin, nPC, meth, flagSparse=FALSE) {

debug <- 0
if (debug==1) {
     Pin <- orig.Pin
     PDFFN <- "Node0.PCA.pval.cindex.pdf" 
     nPC <- 5
     meth="ward.D"
}
     newdistmat<-matrix(0,nrow=nrow(Pin$x),ncol=nrow(Pin$x));
     coordmat <- matrix(0,nrow=nrow(Pin$x),ncol=nPC)
    
     print("====== cell-cell distance matrix ======")
#    pdf("testnPC1.pdf")
     for (ii in 1:nPC) {
          tempdist<-outer(Pin$x[,ii],Pin$x[,ii],"-");
          if (meth=='ward.D') {
              newdistmat<-newdistmat+ (abs(tempdist)*(summary(Pin)$importance[2,ii]))^2  
          } else { 
              newdistmat<-newdistmat+ (abs(tempdist)*(summary(Pin)$importance[2,ii]))
          }
          coordmat[,ii]<-Pin$x[,ii]*summary(Pin)$importance[2,ii]
#         hrrreatmap.2(zscore(newdistmat), trace="none", Rowv=FALSE, Colv=FALSE, col=rbg)
          print(paste("PC", ii, " eigenvalue=", summary(Pin)$importance[2,ii]))

     }

     if (flagSparse) {
         perm.out1 <- HierarchicalSparseCluster.permute(Pin$x[,1:nPC], wbounds=c(1.5,2:6), nperms=5)
         hhh_sparse <- HierarchicalSparseCluster(dists=perm.out1$dists, wbound=perm.out1$bestw, method="complete", dissimilarity="squared.distance")
         hhh_PCA <- hhh_sparse$hc
     } else {
         hhh_PCA <- hclust(as.dist(newdistmat),method="ward.D")
     }
     print("====== binary partition L R and its significance (binary_partition_PC)======")
     hhh_grpLR <- cutree(hhh_PCA,k=2) ;
     hhh_ttt <- table(hhh_grpLR)

     OUT <- list()
     if (min(hhh_ttt) <= 1) {
         OUT$cindex <- 1 
         OUT$pval   <- 1
         print(paste("Cindex=", OUT$cindex, "    pval=", OUT$pval))

     } else {
         if (nPC==1) {
             sigval <- sigclust(cbind(coordmat,1+coordmat), 100, labflag=0, icovest=2) 
         } else {
             sigval <- sigclust(coordmat, 100, labflag=0, icovest=2) 
         }
         OUT$cindex <- sigval@xcindex
         OUT$pval  <-  sigval@pvalnorm
         print(paste("Cindex=", sigval@xcindex, "    pval=", sigval@pvalnorm))
     }

     grp12 <- hhh_grpLR
     if (0) {
         if (sum(hhh_grpLR==1) < sum(hhh_grpLR==2)) {
             grp12 <- 3 - hhh_grpLR
         } else {
             grp12 <- hhh_grpLR
         }
     }

     print("====== DONE! binary partition L R and its significance======")
     OUT$grpLR <- grp12
     OUT$hhh_PCA <- hhh_PCA

     OUT
}


binary_partition_incPC <- function (Zin, Pin, nPC, meth, flagSparse=FALSE, flagDebug=FALSE) {

if (flagDebug) {
     Pin <- orig.Pin
     PDFFN <- "Node0.PCA.pval.cindex.pdf" 
     nPC <- 5
     meth="ward.D"
}

     newdistmat<-matrix(0,nrow=nrow(Pin$x),ncol=nrow(Pin$x));
     coordmat <- matrix(0,nrow=nrow(Pin$x),ncol=nPC)
    
#    pdf("testnPC1.pdf")
     print("====== with increasing number of PC, pick most partitioining nPC  ======")
     print(paste('nPC=', nPC, 'ncol(Pin$x)=', ncol(Pin$x)))
     #if (nPC==1) nPC=nPC + 1
     min.pval <- 0.01
     min.nPC <- 1
     min.grp12 <- c()
     min.hhh_PCA <- c()
     min.cindex <- 0.0
     for (ii in 1:nPC) {
          tempdist<-outer(Pin$x[,ii],Pin$x[,ii],"-");
          if (meth=='ward.D') {
              newdistmat<-newdistmat+ (abs(tempdist)*(summary(Pin)$importance[2,ii]))^2  
          } else { 
              newdistmat<-newdistmat+ (abs(tempdist)*(summary(Pin)$importance[2,ii]))
          }
          coordmat[,ii]<-Pin$x[,ii]*summary(Pin)$importance[2,ii]

          
          if (ii > 1) {
              if (flagSparse) {
                  perm.out1 <- HierarchicalSparseCluster.permute(Pin$x[,1:nPC], wbounds=c(1.5,2:6), nperms=5)
                  hhh_sparse <- HierarchicalSparseCluster(dists=perm.out1$dists, wbound=perm.out1$bestw, method="complete", dissimilarity="squared.distance")
                  hhh_PCA <- hhh_sparse$hc
              } else {
                  hhh_PCA <- hclust(as.dist(newdistmat),method="ward.D")
              }
              hhh_grpLR <- cutree(hhh_PCA,k=2) ;
              hhh_ttt <- table(hhh_grpLR)
 
            if (min(hhh_ttt) > 1) {
                  sigval <- sigclust(coordmat, 100, labflag=0, icovest=2) 
                  if (0) {
                     print("kyle_fem")
                     femval <- kyle_fem(as.matrix(Zin), hhh_grpLR) 
                     print(paste("sigclust", sigval@pvalnorm, sigval@xcindex, ", femval", femval@icl, femval@k1icl, femval@pvalnorm, femval@k))
                     print(paste("sigclust", sigval@pvalnorm, sigval@xcindex))
                  }
                  if (sigval@pvalnorm < min.pval) {
                      min.pval <- sigval@pvalnorm
                      min.cindex <- sigval@xcindex
                      min.grp12 <- hhh_grpLR 
                      min.hhh_PCA <- hhh_PCA 
                      min.nPC <- ii
                  }
             } else {
                  min.pval <- NA 
                  min.cindex <- NA
              }
          #print(paste("PC", ii, " eigenvalue=", summary(Pin)$importance[2,ii], "min.pval=",min.pval))
         }
     }

     OUT <- list()
     OUT$cindex  <- min.cindex
     OUT$pval    <- min.pval
     OUT$grpLR   <- min.grp12
     OUT$hhh_PCA <- min.hhh_PCA
     OUT$nPC <- min.nPC
     OUT$nPCsig <- min.nPC

     print(paste("Cindex=", OUT$cindex, "    pval=", OUT$pval))

     OUT
}


binary_partition_incPC_debug <- function (Zin, Pin, nPC, meth, flagSparse=FALSE, flagDebug=FALSE) {

if (flagDebug) {
     Pin <- orig.Pin
     PDFFN <- "Node0.PCA.pval.cindex.pdf" 
     nPC <- 5
     meth="ward.D"
}

     newdistmat<-matrix(0,nrow=nrow(Pin$x),ncol=nrow(Pin$x));
     coordmat <- matrix(0,nrow=nrow(Pin$x),ncol=nPC)
    
#    pdf("testnPC1.pdf")
     print("====== with increasing number of PC, pick most partitioining nPC  ======")
     #print(paste('nPC=', nPC, 'ncol(Pin$x)=', ncol(Pin$x)))
     #if (nPC==1) nPC=nPC + 1
     min.pval <- 1
     min.nPC <- 1
     min.grp12 <- c()
     min.hhh_PCA <- c()
     min.cindex <- 0.0
     for (ii in 1:nPC) {
          tempdist<-outer(Pin$x[,ii],Pin$x[,ii],"-");
          if (meth=='ward.D') {
              newdistmat<-newdistmat+ (abs(tempdist)*(summary(Pin)$importance[2,ii]))^2  
          } else { 
              newdistmat<-newdistmat+ (abs(tempdist)*(summary(Pin)$importance[2,ii]))
          }
          coordmat[,ii]<-Pin$x[,ii]*summary(Pin)$importance[2,ii]
      #    print(ii)

          
          if (ii > 1) {
              if (flagSparse) {
                  perm.out1 <- HierarchicalSparseCluster.permute(Pin$x[,1:nPC], wbounds=c(1.5,2:6), nperms=5)
                  hhh_sparse <- HierarchicalSparseCluster(dists=perm.out1$dists, wbound=perm.out1$bestw, method="complete", dissimilarity="squared.distance")
                  hhh_PCA <- hhh_sparse$hc
              } else {
                  hhh_PCA <- hclust(as.dist(newdistmat),method="ward.D")
              }
              hhh_grpLR <- cutree(hhh_PCA,k=2) ;
              hhh_ttt <- table(hhh_grpLR)
      #      print(hhh_ttt) 
            if (min(hhh_ttt) > 1) {
                  sigval <- sigclust(coordmat, 100, labflag=0, icovest=2) 
                  if (1) {
      #               print("kyle_fem")
                     femval <- kyle_fem(as.matrix(Zin), hhh_grpLR) 
                     print(paste("sigclust", sigval@pvalnorm, sigval@xcindex, ", femval", femval@icl, femval@k1icl, femval@pvalnorm, femval@k))
                  } else {
                     print(paste("sigclust", sigval@pvalnorm, sigval@xcindex))
                  }
                  if (sigval@pvalnorm < min.pval) {
                      min.pval <- sigval@pvalnorm
                      min.cindex <- sigval@xcindex
                      min.grp12 <- hhh_grpLR 
                      min.hhh_PCA <- hhh_PCA 
                      min.nPC <- ii
                  }
             } else {
                  min.pval <- NA 
                  min.cindex <- NA
              }
          #print(paste("PC", ii, " eigenvalue=", summary(Pin)$importance[2,ii], "min.pval=",min.pval))
         }
     }

     OUT <- list()
     OUT$cindex  <- min.cindex
     OUT$pval    <- min.pval
     OUT$grpLR   <- min.grp12
     OUT$hhh_PCA <- min.hhh_PCA
     OUT$nPC <- min.nPC
     OUT$nPCsig <- min.nPC

#     print(paste("Cindex=", OUT$cindex, "    pval=", OUT$pval))

     OUT
}

part2grp_PCA <- function(Pin, Zin, PDFFN, Nsim=1000, flagSparse, flagPlot, flagSEL, flagDebug=FALSE) {

if (flagDebug) {
    Pin <- orig.Pin
    Zin <- orig.ZSel35
    PDFFN <- "Node0.PCA.pval.cindex.pdf" 
}
    meth <- "ward.D"
    sigval_min <- 1.0
    xcindex_min <- 1.0
   
    print(paste("      get the number of significant PCs by shuffling", Nsim, " ======") )
    if (flagSEL=="jackstraw") {
         print(dim(Zin))
         permPC <- permutationPA(Zin, B = 100, threshold = 0.01)
         if (is.na(permPC$r)) { NN1 <- 0
         } else {
             NN1 <- permPC$r
         }
    } else {
        NN1 <- sum(get_signif_PC(Zin, Pin, 100))
        print(paste("NN1=",NN1))

        flagINC = FALSE
       # if ( NN1 ==  0 ) {
            if (ncol(Pin$x) < 10) { 
                NN1 <- ncol(Pin$x) 
            } else {
                NN1 <- 10
            }
       # }
    }
    #if (NN1 > 10) NN1 <- 10
    print(paste("      Number of Signigicant PC's", NN1))

    if ( NN1 > 0 ) {
        partLR <- binary_partition_incPC (Zin, Pin, NN1, "ward.D", flagSparse, flagDebug) 
        print("Grouping is done....")
        OUT <- list()
        OUT$grp12 <- partLR$grpLR 
        OUT$PCA_grp12 <- partLR$grpLR
        OUT$hhh_PCA <- partLR$hhh_PCA
        OUT$cindex <- partLR$cindex
        OUT$pval <- partLR$pval
        OUT$gap <- 1
        OUT$nPC <- partLR$nPC
        OUT$nPCsig <- partLR$nPC
    } else {
        OUT <- list()
        OUT$nPC <- 0
        OUT$nPCsig <- 0
        OUT$cindex <- NA
        OUT$gap <- 0 
        OUT$pval <- NA
    }

    OUT
}

     
binary_partition_DLM <- function (Zin, Pin, nPC, meth, flagSparse=FALSE, flagDebug=FALSE) {

if (flagDebug) {
     Pin <- orig.Pin
     PDFFN <- "Node0.PCA.pval.cindex.pdf" 
     nPC <- 5
     meth="ward.D"
}
     if (ncol(Pin$x) < 10) { 
         nPC <- ncol(Pin$x) 
     } else {
         nPC <- 10
     }

     newdistmat<-matrix(0,nrow=nrow(Pin$x),ncol=nrow(Pin$x));
     coordmat <- matrix(0,nrow=nrow(Pin$x),ncol=nPC)
    
#    pdf("testnPC1.pdf")
     print("====== binary partition L R and its significance======")
     min.pval <- 1
     min.nPC <- 1
     min.grp12 <- c()
     min.hhh_PCA <- c()
     min.cindex <- 0.0
     for (ii in 1:nPC) {
          tempdist<-outer(Pin$x[,ii],Pin$x[,ii],"-");
          if (meth=='ward.D') {
              newdistmat<-newdistmat+ (abs(tempdist)*(summary(Pin)$importance[2,ii]))^2  
          } else { 
              newdistmat<-newdistmat+ (abs(tempdist)*(summary(Pin)$importance[2,ii]))
          }
          coordmat[,ii]<-Pin$x[,ii]*summary(Pin)$importance[2,ii]

          
          if (ii > 1) {
              if (flagSparse) {
                  perm.out1 <- HierarchicalSparseCluster.permute(Pin$x[,1:nPC], wbounds=c(1.5,2:6), nperms=5)
                  hhh_sparse <- HierarchicalSparseCluster(dists=perm.out1$dists, wbound=perm.out1$bestw, method="complete", dissimilarity="squared.distance")
                  hhh_PCA <- hhh_sparse$hc
              } else {
                  hhh_PCA <- hclust(as.dist(newdistmat),method="ward.D")
              }
              hhh_grpLR <- cutree(hhh_PCA,k=2) ;
              hhh_ttt <- table(hhh_grpLR)
 
            if (min(hhh_ttt) > 1) {
                  femval <- kyle_fem(as.matrix(Zin), hhh_grpLR) 
                  sigval <- sigclust(coordmat, 100, labflag=0, icovest=2) 
                  print(paste("sigclust", sigval@pvalnorm, sigval@xcindex, ", femval", femval@icl, femval@k1icl, femval@pvalnorm, femval@k))
                  if (femval@pvalnorm < min.pval) {
                      min.pval <- femval@pvalnorm
                      min.cindex <- femval@icl
                      min.grp12 <- hhh_grpLR 
                      min.hhh_PCA <- hhh_PCA 
                      min.nPC <- ii
                  }
             } else {
                  min.pval <- NA 
                  min.cindex <- NA
              }
          #print(paste("PC", ii, " eigenvalue=", summary(Pin)$importance[2,ii], "min.pval=",min.pval))
         }
     }

     OUT <- list()
     OUT$cindex  <- min.cindex
     OUT$pval    <- min.pval
     OUT$grpLR   <- min.grp12
     OUT$hhh_PCA <- min.hhh_PCA
     OUT$nPC <- min.nPC
     OUT$nPCsig <- min.nPC

     print(paste("Cindex=", OUT$cindex, "    pval=", OUT$pval))

     OUT
}

part2grp_DLM <- function(Pin, Zin, PDFFN, Nsim=1000, flagSparse, flagPlot, flagDebug=FALSE) {

if (flagDebug) {
    Pin <- orig.Pin
    Zin <- orig.ZSel35
    PDFFN <- "Node0.PCA.pval.cindex.pdf" 
}
    meth <- "ward.D"
    sigval_min <- 1.0
    xcindex_min <- 1.0
   
    print(paste("      get the number of significant PCs by shuffling", Nsim, " ======") )
    NN1 <- ncol(Pin$x)
    #NN1 <- sum(get_signif_PC(Zin, Pin, Nsim))
    print(paste("      Nfeature", NN1))

    if (NN1 > 0) {
        #partLR <- binary_partition_PC (Pin, NN1, "ward", flagSparse) 
        partLR <- binary_partition_DLM (Zin, Pin, NN1, "ward.D", flagSparse, flagDebug) 
        print("Grouping is done....")
        OUT <- list()
        OUT$grp12 <- partLR$grpLR 
        OUT$PCA_grp12 <- partLR$grpLR
        OUT$hhh_PCA <- partLR$hhh_PCA
        OUT$cindex <- partLR$cindex
        OUT$pval <- partLR$pval
        OUT$gap <- 1
        OUT$nPC <- partLR$nPC
        OUT$nPCsig <- partLR$nPC

    } else {

        OUT <- list()
        OUT$nPC <- 0
        OUT$nPCsig <- 0
        OUT$cindex <- NA
        OUT$gap <- 0 
        OUT$pval <- NA
    }

    OUT
}

get_grp_mean_cov <- function(Xin, grp12, Npc, featnames) {

    if (is.list(Xin) & "x" %in% names(Xin)) {
        PCn <- Xin$x
        loadings <- Xin$rotation
    } else {
        PCn <- Xin
        loadings <- c()
    }

    idx1 <- which(grp12==1)
    idx2 <- which(grp12==2)

    X1 <- as.matrix(PCn[idx1,1:Npc])
    X2 <- as.matrix(PCn[idx2,1:Npc])

    mu1 <- apply(X1, 2, mean)
    cov1 <- cov(X1)

    mu2 <- apply(X2, 2, mean)
    cov2 <- cov(X2)

    OUT <- list()
    OUT$mu1 <- mu1
    OUT$cov1 <- cov1
    OUT$mu2 <- mu2
    OUT$cov2 <- cov2
    OUT$loadings <- loadings
    OUT$featnames <- featnames

    OUT
}


cal_12membership_mean_sd <- function (X, param, Nin) {

    N = length(param$mu1)
    if (Nin>N) print(paste("Error", Nin, "is bigger than", N))

    p1 <- p2 <- normp1 <- normp2 <- rep(0,Nin)
    for (i in 1:Nin) {
         if (param$mu1[i] < param$mu2[i]) {
             p1[i] <- 1 - pnorm(X[i], param$mu1[i], param$sd1[i])
             p2[i] <- pnorm(X[i], param$mu2[i], param$sd2[i])
         } else {
             p1[i] <- pnorm(X[i], param$mu1[i], param$sd1[i])
             p2[i] <- 1 - pnorm(X[i], param$mu2[i], param$sd2[i])
         }
         normp1[i] <- p1[i]/(p1[i] + p2[i])
         normp2[i] <- p2[i]/(p1[i] + p2[i])
    }

    norm1 <- norm(as.matrix(normp1), type="F")
    norm2 <- norm(as.matrix(normp2), type="F") 
    mem1 <- norm1 / (norm1 + norm2)
    mem2 <- 1-mem1
    
    OUT <- list()
    OUT$p1 <- p1
    OUT$p2 <- p2
    OUT$normp1 <- normp1
    OUT$normp2 <- normp2
    OUT$mem1 <- mem1
    OUT$mem2 <- mem2

    OUT 
}


cal_12membership_PC <- function (XP, param, Nin, flagPC=TRUE) {
    if (flagPC) { 
        if (is.list(XP) & "x" %in% names(XP)) {
            P <- XP$x
        } else {
            P <- XP[, param$featnames] %*% param$loadings
        }
    } else {
        P <- XP[, param$featnames] 
    }
    X <- P[, 1:Nin]
  
    if (is.matrix(X)) {
        Nsample <- nrow(X)
    } else {
        Nsample <- length(X)
    }
    Nfeat = length(param$mu1)
    if (Nin>Nfeat) print(paste("Error", Nin, "is bigger than", Nfeat))

    if (0) {
        MU1 <- t(matrix(rep(param$mu1, Nsample), ncol=Nsample))
        MU2 <- t(matrix(rep(param$mu2, Nsample), ncol=Nsample))
        if (is.singular.matrix(param$cov1) | is.singular.matrix(param$cov2)) {
            print("param$cov1 or param$cov2 is singular")
            param$cov1 <- param$cov2 <- (param$cov1 + param$cov2)/2
        }
        DET1 <- det(param$cov1)
        INVCOV1 <- solve(param$cov1)
        DET2 <- det(param$cov2)
        INVCOV2 <- INVCOV1 #solve(param$cov2)
    
        dy1 = t(X - MU1)
        dy2 = t(X - MU2)
    
        p1 <- p2 <- mem1 <- mem2 <- rep(0,Nsample)
        for (i in 1:Nsample) {
             p1[i] <- exp(((t(dy1[,i]) %*% INVCOV1) %*% dy1[,i])*(-1/2)) / ((2*pi)^(Nfeat/2) * (1/(DET1)^(1/2)))   
             p2[i] <- exp(((t(dy2[,i]) %*% INVCOV2) %*% dy2[,i])*(-1/2)) / ((2*pi)^(Nfeat/2) * (1/(DET2)^(1/2)))   
             mem1[i] <- p1[i]#/(p1[i] + p2[i])
             mem2[i] <- p2[i]#/(p1[i] + p2[i])
        }
     } else {
        q1 <- pmnorm(X, param$mu1, param$cov1)
        q2 <- pmnorm(X, param$mu2, param$cov2)
     }

    OUT <- list()
    OUT$p1 <- q1
    OUT$p2 <- q2
    OUT$mem1 <- q1 / (q1 + q2)
    OUT$mem2 <- q2 / (q1 + q2)

    return(cbind(OUT$mem1, OUT$mem2))
}

cal_probability_nc <- function (Nodei, X) {

    Nin <- length(Nodei$membership.params$mu1) 
    out <- cal_12membership_PC (X, Nodei$membership.params, Nin, flagPC=TRUE) 
    return(out) 
}

cal_membership_node <- function (Nodei, X, nodestr, partitions=c("L","R"), branches=c("Left", "Right") ) {

    print(nodestr)
    NC <- length(branches) 
    if (nchar(nodestr) > 1) {
        
       nodestr.i        <- substr(nodestr, 1, 1)
       nodestr.i_remain <- substr(nodestr, 2, nchar(nodestr))

       branch.i <- which(partitions==nodestr.i)
       membership.i <- cal_probability_nc (Nodei, X)[,branch.i] 
       membership <- membership.i * cal_membership_node (Nodei[[branches[branch.i]]], X, nodestr.i_remain, partitions, branches)


    } else {
       nodestr.i <- nodestr
       branch.i <- which(partitions==nodestr.i)
       membership <- cal_probability_nc (Nodei, X)[,branch.i] 
      
    }

    return(membership)
}

clustering2decR <- function (ORIG, orig.XSel,orig.GT_cre, orig.ZSel35, orig.strin, mypch, leg.pch, crecolor, crecolor2, leg.crecolor, leg.crecolor2, leg.str, leg.str2, nodestr, nsim, flagDEC="LDA", flagGRP="SEL", flagPthr=0.05, flagSparse=FALSE, flagPlot=TRUE, flagDIP=FALSE, flagMembership=FALSE, flagPartition="PCA", flagSEL="jackstraw", flagPRE=TRUE, flagDebug=FALSE) { 

     if (flagDebug) {
        ORIG <- Nodei
        orig.GT_cre <- Nodei$GT_cre
        orig.ZSel35 <- Nodei$ZSel35
        orig.XSel <- Nodei$Xsoi
        orig.strin <- "hello"
        mypch <- Nodei$pch
        leg.pch <- soi.leg$pch
        crecolor <- Nodei$crecolor
        crecolor2 <- Nodei$hybrid.crecolorL 
        leg.crecolor <- soi.leg$crecolor
        leg.crecolor2 <- soi.leg$crecolorL
        leg.str <-  soi.leg$str
        leg.str2 <- soi.leg$strL
        nodestr <- "debug"
        nsim <- 100
        flagDEC="LDA"
        flagGRP="SEL"
        flagPthr=0.05 
        flagSparse=FALSE
        flagPlot=FALSE    
     }
     set.seed(1)
     OUT <- ORIG

     OUT$nPC  <- 0 #print("this one is added")
     #print(dim(orig.ZSel35))

     pdffn <- paste(orig.strin, ".pca.pdf", sep="") 
     if (!flagPlot) save(ORIG, file="ORIG.Rdata") 
     if (flagPlot) pdf(pdffn)

     if (is.na(flagPRE)) {
         print("      using all features without Pre-selection")
         foi = 1:ncol(orig.ZSel35)

     } else {

         print("      selecting features with bimodal distribution")
         nfeat <- ncol(orig.ZSel35) 
         print(paste("           number of cells :", nrow(orig.ZSel35), ",     number of features :", ncol(orig.ZSel35)))

         if (flagPRE=="diptest") {
             foi <- c() 
             for (f in c(1:nfeat)) {
                  frun <- dip.test(orig.ZSel35[,f], simulate.p.value=TRUE, B=nsim)
                  if (frun$p.value < flagPthr) foi <- c(foi, f)
             } 
         } 
         if (flagPRE=="sigclust") {
             varfeat <- apply(orig.ZSel35, 2, sd)
             zero.sd <- sum(varfeat==0)

             if (zero.sd==0) { f0 <- 2} else { f0 <- 1 }
             foi <- c() 
             if (zero.sd!=nfeat) {
                 varorder <- order(varfeat, decreasing=FALSE)
                 for (f in c((zero.sd+f0):nfeat)) {
                      Zf <- orig.ZSel35[,c(varorder[1],varorder[f])]
                      fsig <- sigclust(Zf, nsim, labflag=0, label=0, icovest=2)
                      if (fsig@pvalnorm < flagPthr) {
                          foi <- c(foi, varorder[f]) 
                      }
                 }
             }
         }
         print(paste("           number of interesting features = ", length(foi), "out of", nfeat))
     }

     print("     Do pca")
     if (length(foi)>0) {

         ZSel35 <- orig.ZSel35[, foi]
         XSel <- orig.XSel[, foi]
         orig.pin <- prcomp(ZSel35, cor=TRUE)
    
         if (flagPartition=="PCA") {
             print("     Partition in pca")
             pdffn <- paste(orig.strin, ".pca.pdf", sep="") 
             orig.grp <- part2grp_PCA(orig.pin, ZSel35, pdffn, nsim, flagSparse, flagPlot, flagSEL, flagDebug)
         } else {
             pdffn <- paste(orig.strin, ".DLM.pdf", sep="") 
             orig.grp <- part2grp_DLM(orig.pin, ZSel35, pdffn, nsim, flagSparse, flagPlot, flagDebug)
         }
         if (flagPlot) dev.off()
    
         pdffn <- paste(orig.strin, ".pca.clustering.pdf", sep="")

         print(paste("           (orig.grp$nPCsig, orig.grp$pval, orig.grp$pval)=",orig.grp$nPCsig, round(orig.grp$pval,3), round(orig.grp$pval,3))) 
         if ((orig.grp$nPCsig > 0) && !is.na(orig.grp$cindex) && !is.na(orig.grp$pval) && (orig.grp$pval< flagPthr)) {

             clusterPCA (orig.pin, orig.grp$nPCsig, crecolor2, leg.crecolor2, leg.str2, pdffn, flagPlot) 
             print("            any features with discrimination power given binary partition?")
             if (flagSEL=="NA") idx.ordered <- order_feat_byNA (ZSel35, orig.grp$nPCsig, XSel)
             if (flagSEL=="jackstraw") idx.ordered <- order_feat_byJackstraw (ZSel35, orig.grp$nPCsig, XSel)
             if (flagSEL=="ttest") idx.ordered <- order_feat_byTtest (ZSel35, orig.grp$PCA_grp12, XSel)

             ZSel35.pordered <- ZSel35[, idx.ordered$order]
             tmp.pordered <- ZSel35.pordered 
             colnames(tmp.pordered) <- paste(signif(-log10(idx.ordered$pval.ordered),2), colnames(ZSel35)[idx.ordered$order])
             if (flagPthr==0.01) nadjp <- idx.ordered$pval.lt0.01
             if (flagPthr==0.05) nadjp <- idx.ordered$pval.lt0.05 
             print(paste("           ", nadjp, "features with adjp<", flagPthr)) 

             OUT$ttest.feat <- idx.ordered 

             if (nadjp > 0) {
        
                 write.csv(data.frame(feature=colnames(ZSel35)[idx.ordered$order[1:nadjp]], log10P.neg=signif(-log10(idx.ordered$pval.ordered[1:nadjp]),4)), file= paste(orig.strin, '.PCApartition.AllFeatures.adjP', flagPthr, '.sorted.csv', sep=""), row.names=FALSE)
            
                 if (flagPlot) {
                     pdf(paste(orig.strin, '.heatmap.pcapartition.pdf', sep=""), height=5+round(ncol(ZSel35)/5), width=8+round(nrow(ZSel35)/10))
                     hh <- heatmap.3(t(tmp.pordered), hclustfun=hhclust, Colv=as.dendrogram(orig.grp$hhh_PCA), Rowv=TRUE, trace="none",  dendrogram="column",
                          ColSideColors=crecolor2, main=paste("clustering of cells @", nodestr, " with ", orig.grp$nPCsig," pc's \n", 
                          " cindex =", signif(orig.grp$cindex,3), " & -log10(p) =", signif(-log10(orig.grp$pval),3), "\n",
                          ncol(ZSel35), " features", sep=""), 
                          keysize=0.8, margins=c(10,10), cexcol=0.5,cexrow=1, col=rbg)
                     legend("bottomleft", bg="white", legend=leg.str2, pch=15, col=leg.crecolor2, cex=0.75, bty="n")
                 }
    
                 hhh_allFeat <- hclust(dist(tmp.pordered),method="ward.D")
                 SelFeat_grpLR <- cutree(hhh_allFeat,k=2)
            
                 if (nadjp > 1) {
                     adjp0.01.ZSel35 <- ZSel35.pordered[, 1:nadjp]
                     tmp.adjp.pordered <- tmp.pordered[, 1:nadjp]
        
                     if (flagPlot) {    
                         heatmap.3(t(tmp.adjp.pordered), hclustfun=hhclust, Colv=as.dendrogram(orig.grp$hhh_PCA), Rowv=TRUE, trace="none",  dendrogram="column",
                                  ColSideColors=crecolor2, main=paste("clustering of cells @", nodestr, " with ", orig.grp$nPCsig," pc's \n", 
                                  " cindex =", signif(orig.grp$cindex,3), " & -log10(p) =", signif(-log10(orig.grp$pval),3), "\n",
                                  " showing ", nadjp, " features w/ t.test adjp <", flagPthr, sep=""), 
                                  keysize=0.8, margins=c(10,10), cexcol=0.5,cexrow=1, col=rbg)
                         legend("bottomleft", bg="white", legend=leg.str2, pch=15, col=leg.crecolor2, cex=0.75, bty="n")
                         hh <- heatmap.3(t(tmp.adjp.pordered), hclustfun=hhclust, Colv=TRUE, Rowv=TRUE, trace="none",  dendrogram="column",
                                  ColSideColors=crecolor2, main=paste("clustering of cells @", nodestr, " with ", orig.grp$nPCsig," pc's \n", 
                                  " cindex =", signif(orig.grp$cindex,3), " & -log10(p) =", signif(-log10(orig.grp$pval),3), "\n",
                                  " clustered with ", nadjp, " features w/ t.test adjp <", flagPthr, sep=""), 
                                  keysize=0.8, margins=c(10,10), cexcol=0.5,cexrow=1, col=rbg)
                         legend("bottomleft", bg="white", legend=leg.str2, pch=15, col=leg.crecolor2, cex=0.75, bty="n")
                     }
    
                     hhh_SelFeat <- hclust(dist(tmp.adjp.pordered),method="ward.D")
                     hhh_1 <- hclust(dist(t(tmp.adjp.pordered)),method="ward.D")
                     SelFeat_grpLR <- cutree(hhh_SelFeat,k=2)
                     OUT$hhh_SEL <- hhh_SelFeat 
                     OUT$hhh_SEL_rowInd <- hhh_SelFeat$order 
                     OUT$hhh_SEL_colInd <- hhh_1$order 
                 } 
                 if (flagPlot) dev.off()
    
           if (flagMembership) { 
                 print("          estimate Membership parameter for two clusters")
                 membership.params <- get_grp_mean_cov (orig.pin, orig.grp$PCA_grp12, orig.grp$nPCsig, colnames(ZSel35))
                 print("          calculate memebership score")
                 #membershiplrs <- apply(as.matrix(orig.pin$x[,1:orig.grp$nPCsig]), 1, cal_12membership_mean_sd, membership.params, orig.grp$nPCsig)
                 membershiplrs <- cal_12membership_PC(orig.pin, membership.params, orig.grp$nPCsig, TRUE)
           }
                
            
                 OUT$groupingSelFeat  <- SelFeat_grpLR 
                 OUT$groupingPC       <- orig.grp$PCA_grp12 
                # print(table( OUT$groupingSelFeat ))
                # print(table( orig.grp$PCA_grp12 ))
                 if (flagGRP == "SEL") OUT$grouping  <- OUT$groupingSelFeat
                 if (flagGRP == "PC") OUT$grouping  <- orig.grp$PCA_grp12 
                 OUT$hhh_PCA <- orig.grp$hhh_PCA 
                 OUT$nPC     <- orig.grp$nPCsig 
    
    
                 print("          get up to 10 best features dividing these two groups")
                
                 orig.class <- orig.grp$PCA_grp12
                 order_pc1_loading <- order(abs(orig.pin$rotation[,"PC1"]), decreasing=TRUE)
                 orig_feat_pcloading <- rownames(orig.pin$rotation)[order_pc1_loading]
                 if (flagDEC=="SVM") {
                 print("          by SVM decision rule")
                     orig.pcloading.cv4roc <- mystepbystepSVMcv(orig.class, orig.class, ZSel35, orig_feat_pcloading, 
             	        length(orig_feat_pcloading), 4, 1.5, paste(orig.strin, ".roc.gty2.trainy2.feature.thr.1.5.pdf", sep=""), 0 )
                 } 
                 if (flagDEC=="LDA") {
                 print("          by LDA decision rule")
                     orig.pcloading.cv4roc <- mystepbystepcv(orig.class, orig.class, ZSel35, orig_feat_pcloading, 
             	        length(orig_feat_pcloading), 4, 1.5, paste(orig.strin, ".roc.gty2.trainy2.feature.thr.1.5.pdf", sep=""), 0 )
                 } 
    
                 selected_featnames10 <- orig.pcloading.cv4roc$stepbystep.featname[,4+1]
                 #print(selected_featnames10) 
                 sf10.ZSel35 <- ZSel35[, selected_featnames10]
                 print("          parameter for two clusters are estimated")
                 sf10.params <- get_grp_mean_cov (sf10.ZSel35, orig.grp$PCA_grp12, orig.pcloading.cv4roc$nf, colnames(sf10.ZSel35))
            
                 OUT$nPC         <- orig.grp$nPCsig 
                 OUT$loadingsPC  <- orig.pin$rotation[,1:orig.grp$nPCsig] 
                 OUT$groupingPC  <- orig.grp$PCA_grp12 
                 OUT$grouping_pval_cindex_pc <- c(orig.grp$pvalue, orig.grp$cindex)
            
           #if (flagMembership) {
           if (0) {
                 print("calculate memebership score")
                 sf10.membershiplrs <- list()
                 sf10.significance <- matrix(0,nrow=2,ncol=orig.pcloading.cv4roc$nf)
                 rownames(sf10.significance) <- c("pvalue", "cindex")
                 for (nsf in 1:orig.pcloading.cv4roc$nf) {
                      #print(paste('number of feature=', nsf))
                      #tmp.membershiplrs <- apply(as.matrix(sf10.ZSel35[,1:nsf]), 1, cal_12membership_mean_sd, sf10.params, nsf)
                      tmp.membershiplrs <- cal_12membership_PC(sf10.ZSel35[,1:nsf], sf10.params, nsf, FALSE)
                      tmp.call <- matrix(1,1, ncol(tmp.membershiplrs))
                      tmp.call[1,tmp.membershiplrs[1,] < tmp.membershiplrs[2,]] <- 2
                      sf10.membershiplrs[[nsf]] <- t(rbind(tmp.membershiplrs, tmp.call))
            
                      if (nsf==1) {
                          tmp.pcindex <- sigclust(sf10.ZSel35[,c(1,1)],2*ncol(tmp.call), labflag=1, label=tmp.call,icovest=2) 
                      } else {
                          tmp.pcindex <- sigclust(sf10.ZSel35[,1:nsf],2*ncol(tmp.call), labflag=1, label=tmp.call,icovest=2) 
                      }
                      sf10.significance[, nsf] <- c(tmp.pcindex@pvalnorm, tmp.pcindex@xcindex)
                   }
                   OUT$params_sf10 <- sf10.params
                   OUT$featname_sf10 <- selected_featnames10
                   OUT$membership_sf10 <- sf10.membershiplrs 
                   OUT$membership.params_sf10 <- sf10.params 
                   OUT$grouping_pval_cindex_sf10 <- sf10.significance 
          } else {

                 OUT$params_sf10 <- c()
                 OUT$featname_sf10 <- c()
                 OUT$membership_sf10 <- c() 
                 OUT$grouping_pval_cindex_sf10 <- c() 
          } 
                 OUT$params_pc   <- params
                 OUT$membership <- membershiplrs
                 OUT$membership.params <- membership.params
    
                 
                 if (flagPlot) {
                     print("calculate creline composition & draw pie chart")
    
                     print(length(orig.GT_cre)) 
                     print(length(OUT$grouping)) 
                     OUT$ct <- table(orig.GT_cre, OUT$grouping)
                     print(OUT$ct)
                     write.csv(OUT$ct, file=paste(orig.strin, ".cre_line.composition.csv", sep=""))
    
                     pdf(paste(orig.strin, ".cre_line.composition.pie.pdf", sep=""))
                     par(mfrow=c(2,1))
                     if (0) {
                         if (sum(OUT$ct[,1]) > sum(OUT$ct[,2])) {
                             piel <- 2 ; pier <- 1
                         } else {
                             piel <- 1 ; pier <- 2
                         }
                     } else {
                         piel <- 1 ; pier <- 2
                     }
                     pie( OUT$ct[,piel], labels = names(OUT$ct[,piel]), 
                          col=leg.crecolor[match(names(OUT$ct[,piel]), leg.str)],
                          main=paste("Cell Composition @ ", nodestr, "L", sep=""))
                     pie( OUT$ct[,pier], labels = names(OUT$ct[,pier]), 
                          col=leg.crecolor[match(names(OUT$ct[,pier]), leg.str)],
                          main=paste("Cell Composition @ ", nodestr, "R", sep=""))
                     dev.off()
                 }
    
             } else {
                 print("                !!!! terminal node # no feature with between two groups")
                 OUT$nPC  <- 0
             }
                 
         } else {
             if ((orig.grp$gap)==0) {
                  print("                !!!! terminal node # gap stat says STOP ")
             } 
             if (is.na(orig.grp$pval) && is.na(orig.grp$cindex)) {
                  print("                !!!! terminal node # no pc with eigenvalue higher than shuffled")
             } else {
                  print("                !!!! terminal node # lr branching is 1:(n-1)")
             }
             OUT$nPC  <- 0
         }
     } else {
         OUT$nPC  <- 0
     }
    
     OUT
}


mystepbystepSVMcv <- function(gt, y, zin, featnameallshort, nfeat, nfold, thr, rocpdffn, debug=0, flagPlot=TRUE) {
    library(e1071)

#    if (debug) {
#        gt <- orig.class
#        y <- orig.class
#        zin <- orig.ZSel35
#        featnameallshort <- orig_feat_pcloading
#        nfeat <-38
#        nfold <- 4
#        thr <- 1.5 
#        rocpdffn <- 'ck2test.0.pdf'
#    }

    if (nfeat <= 10) {
        nf=nfeat
    } else {
        nf=10
    }

    print(paste("               number of max feature to choose ", nf))
    print(paste("               number of columns in data ", ncol(zin)))

    y2.idx <- which(y==2)
    y2.n <- length(y2.idx)
    
    y1.idx <- which(y==1)
    y1.n <- length(y1.idx)
    
    err1.cv <- err2.cv <- err3.cv <- err4.cv <- err5.cv <- err6.cv <- err7.cv <- err8.cv <- err9.cv <- err10.cv <- rep(0,nfold+1) 
    err1.resub <- err2.resub <- err3.resub <- err4.resub <- err5.resub <- err6.resub <- err7.resub <- err8.resub <- err9.resub <- err10.resub <- rep(0,nfold+1) 
    ct1.cv <- ct2.cv <- ct3.cv <- ct4.cv <- ct5.cv <- ct6.cv <- ct7.cv <- ct8.cv <- ct9.cv <- ct10.cv <- matrix(0,ncol=nfold+1,nrow=2)
    ct1.resub <- ct2.resub <- ct3.resub <- ct4.resub <- ct5.resub <- ct6.resub <- ct7.resub <- ct8.resub <- ct9.resub <- ct10.resub <- matrix(0,ncol=nfold+1,nrow=2) 
    sig1.cv <- sig2.cv <- sig3.cv <- sig4.cv <- sig5.cv <- sig6.cv <- sig7.cv <- sig8.cv <- sig9.cv <- sig10.cv <- matrix(0,ncol=nfold+1,nrow=2)
    sig1.resub <- sig2.resub <- sig3.resub <- sig4.resub <- sig5.resub <- sig6.resub <- sig7.resub <- sig8.resub <- sig9.resub <- sig10.resub <- matrix(0,ncol=nfold+1,nrow=2) 
        
    stepbystep.featname = matrix("hello", 10, nfold+1)
    nsample = nrow(zin)
    
    mymod1 = 1:y1.n %% nfold
    mymod2 = 1:y2.n %% nfold
    out <- list()

    print(paste("               ", nfold, "fold CV + Resub."))
    for (f in 1:(nfold+1)) {
        #idxtr = c(sample(y2.idx,f.samplesize), y1.idx) 
        #dxte = c(setdiff(1:nsample,idxtr), y1.idx)
        if (f > nfold) {
           idxtr <- c(y2.idx, y1.idx)
           idxte <- idxtr
        } else {
           mask1 <- mymod1 == (f-1)
           mask2 <- mymod2 == (f-1)
           idx1 <- idx2 <- c()
           if (sum(mask1) > 0) {
               idx1 <- y1.idx[which(mask1)]
           } else {
               idx1 <- sample(y1.idx,1)
           }
           if (sum(mask2) > 0) {
               idx2 <- y2.idx[which(mask2)]
           } else {
               idx2 <- sample(y2.idx,1)
           }
           
           idxte <- c(idx1, idx2)
           idxtr <- setdiff(1:nsample,idxte)
        }

        trin <- zin[idxtr,]
        tein <- zin[idxte,]
        iintr  <- matrix(1,length(idxtr),1)
        iinte  <- matrix(1,length(idxte),1)
        ytr <- y[idxtr]
        yte <- y[idxte]
        gttr <- gt[idxtr]
        gtte <- gt[idxte]
        
       # if (debug) print("1")

 ##############################################################    
       # print(1)
       # print(paste("size of gttr, ytr",length(ytr)))
       # print(paste(gttr, ytr))

        if (ncol(zin)>=1 && nf>=1) {
           minerr1=1000
           #print(colnames(trin))
           #print(match(featnameallshort,colnames(trin)))
           for (i in 1:nfeat) {
           #     print(paste(i,featnameallshort[i]))
                xi <- trin[, featnameallshort[i]]
                DATAi <- data.frame(X=xi, Y=ytr)

                svmi <- svm(Y ~ ., data=DATAi)
                #svmiResult <- tune(svm, Y ~ ., data=DATAi, ranges=list(epsilon=c(seq(0,1,0.1)),cost=2^(2:9)))
                #svmi <- svmiResult$best.model
                ipred <- predict(svmi, DATAi)
                ierr <- rmse(DATAi$Y-ipred)
                if (ierr < minerr1) {
                    minpred <- ipred
                    minerr1 <- ierr 
                    min1 <- i 
                    minfiti <- svmi
                }
           }
      print("min1 and featnameallshort(min1)")
      print(paste(min1, featnameallshort[c(min1)]))
           xj <- tein[, featnameallshort[c(min1)]]
           DATAj <- data.frame(X=xj, Y=yte)

           est.resub <- minpred
           est.cv    <- predict(minfiti, DATAj)
           #err.cv    <- rmse(DATAj$Y - predict(minfiti, DATAj))

           ct1.cv[1:2,f]    = get_sen_spc(est.cv, gtte, thr)
           ct1.resub[1:2,f] = get_sen_spc(est.resub, gttr, thr)
        
           stepbystep.featname[1,f] = featnameallshort[c( min1)]
     
        }
        

 ##############################################################    
        if (ncol(zin)>=2 && nf>=2) {   
           minerr2=1000;
           for (i in 1:nfeat) {
                if (!(i %in% c(min1))) {
           #     print(paste(i,featnameallshort[i]))
           #         print(paste(featnameallshort[c(min1, i)]))
                    xi <- trin[, featnameallshort[c(min1, i)]]
 
                    DATAi <- data.frame(X=xi, Y=ytr)
                    svmi <- svm(Y ~ ., data=DATAi)
                   # svmiResult <- tune(svm, Y ~ ., data=DATAi, ranges=list(epsilon=seq(0,1,0.1),cost=2^(2:9)))
                   # svmi <- svmiResult$best.model
                    ipred <- predict(svmi, DATAi)
                    ierr <- rmse(DATAi$Y-ipred)
                    if (ierr < minerr2) {
                        minpred <- ipred
                        minerr2 <- ierr 
                        min2 <- i 
                        minfiti <- svmi
                    }
               }
           }
           xj <- tein[, featnameallshort[c(min1, min2)]]
           DATAj <- data.frame(X=xj, Y=yte)
    
           est.resub <- minpred
           est.cv    <- predict(minfiti, DATAj)
           #err.cv    <- rmse(DATAj$Y - predict(minfiti, DATAj))
    
           ct2.cv[1:2,f]    = get_sen_spc(est.cv, gtte, thr)
           ct2.resub[1:2,f] = get_sen_spc(est.resub, gttr, thr)
        
           stepbystep.featname[1:2,f] = featnameallshort[c( min1,min2)]
        }

 ##############################################################    
        if (debug) print("3")
        if (ncol(zin)>=3 && nf>=3) {   
           minerr3=1000;
           for (i in 1:nfeat) {
                if (!(i %in% c(min1,min2))) {

                    xi <-trin[, featnameallshort[c(min1,min2,i)]]
 
                    DATAi <- data.frame(X=xi, Y=ytr)
                    svmi <- svm(Y ~ ., data=DATAi)
                    #svmiResult <- tune(svm, Y ~ ., data=DATAi, ranges=list(epsilon=seq(0,1,0.1),cost=2^(2:9)))
                    #svmi <- svmiResult$best.model
                    ipred <- predict(svmi, DATAi)
                    ierr <- rmse(DATAi$Y-ipred)
                    if (ierr < minerr3) {
                        minpred <- ipred
                        minerr3 <- ierr 
                        min3 <- i 
                        minfiti <- svmi
                    }
               }
           }

           xj <- tein[, featnameallshort[c(min1, min2, min3)]]
           DATAj <- data.frame(X=xj, Y=yte)
    
           est.resub <- minpred
           est.cv    <- predict(minfiti, DATAj)
           #err.cv    <- rmse(DATAj$Y - predict(minfiti, DATAj))
    
           ct3.cv[1:2,f]    = get_sen_spc(est.cv, gtte, thr)
           ct3.resub[1:2,f] = get_sen_spc(est.resub, gttr, thr)
        
           stepbystep.featname[1:3,f] = featnameallshort[c( min1,min2,min3)]
           
        }
 ##############################################################    
        if (ncol(zin)>=4 && nf>=4) {   
           minerr4=1000;
           for (i in 1:nfeat) {
                if (!(i %in% c(min1,min2,min3))) {
                    xi <-trin[, featnameallshort[c(min1,min2,min3,i)]]
 
                    DATAi <- data.frame(X=xi, Y=ytr)
                    svmi <- svm(Y ~ ., data=DATAi)
                    #svmiResult <- tune(svm, Y ~ ., data=DATAi, ranges=list(epsilon=seq(0,1,0.1),cost=2^(2:9)))
                    #svmi <- svmiResult$best.model
                    ipred <- predict(svmi, DATAi)
                    ierr <- rmse(DATAi$Y-ipred)
                    if (ierr < minerr4) {
                        minpred <- ipred
                        minerr4 <- ierr 
                        min4 <- i 
                        minfiti <- svmi
                    }
               }
           }

           xj <- tein[, featnameallshort[c(min1, min2, min3, min4)]]
           DATAj <- data.frame(X=xj, Y=yte)
    
           est.resub <- minpred
           est.cv    <- predict(minfiti, DATAj)
           #err.cv    <- rmse(DATAj$Y - predict(minfiti, DATAj))
    
           ct4.cv[1:2,f]    = get_sen_spc(est.cv, gtte, thr)
           ct4.resub[1:2,f] = get_sen_spc(est.resub, gttr, thr)
        
           stepbystep.featname[1:4,f] = featnameallshort[c( min1,min2,min3,min4)]
           
        }
 ##############################################################    
        if (ncol(zin)>=5 && nf>=5) {   
           minerr5=1000;
           for (i in 1:nfeat) {
                if (!(i %in% c(min1,min2,min3,min4))) {
                    xi <-trin[, featnameallshort[c(min1,min2,min3,min4,i)]]
 
                    DATAi <- data.frame(X=xi, Y=ytr)
                    svmi <- svm(Y ~ ., data=DATAi)
                    #svmiResult <- tune(svm, Y ~ ., data=DATAi, ranges=list(epsilon=seq(0,1,0.1),cost=2^(2:9)))
                    #svmi <- svmiResult$best.model
                    ipred <- predict(svmi, DATAi)
                    ierr <- rmse(DATAi$Y-ipred)
                    if (ierr < minerr5) {
                        minpred <- ipred
                        minerr5 <- ierr 
                        min5 <- i 
                        minfiti <- svmi
                    }
               }
           }

           xj <- tein[, featnameallshort[c(min1, min2, min3, min4, min5)]]
           DATAj <- data.frame(X=xj, Y=yte)
    
           est.resub <- minpred
           est.cv    <- predict(minfiti, DATAj)
           #err.cv    <- rmse(DATAj$Y - predict(minfiti, DATAj))
    
           ct5.cv[1:2,f]    = get_sen_spc(est.cv, gtte, thr)
           ct5.resub[1:2,f] = get_sen_spc(est.resub, gttr, thr)
        
           stepbystep.featname[1:5,f] = featnameallshort[c( min1,min2,min3,min4,min5)]

        }
 ##############################################################    
        if (ncol(zin)>=6 && nf>=6) {   
           minerr6=1000;
           for (i in 1:nfeat) {
                if (!(i %in% c(min1,min2,min3,min4,min5))) {
                    xi <-trin[, featnameallshort[c(min1,min2,min3,min4,min5,i)]]
 
                    DATAi <- data.frame(X=xi, Y=ytr)
                    svmi <- svm(Y ~ ., data=DATAi)
                    #svmiResult <- tune(svm, Y ~ ., data=DATAi, ranges=list(epsilon=seq(0,1,0.1),cost=2^(2:9)))
                    #svmi <- svmiResult$best.model
                    ipred <- predict(svmi, DATAi)
                    ierr <- rmse(DATAi$Y-ipred)
                    if (ierr < minerr6) {
                        minpred <- ipred
                        minerr6 <- ierr 
                        min6 <- i 
                        minfiti <- svmi
                    }
               }
           }

           xj <- tein[, featnameallshort[c(min1, min2, min3, min4, min5, min6)]]
           DATAj <- data.frame(X=xj, Y=yte)
    
           est.resub <- minpred
           est.cv    <- predict(minfiti, DATAj)
           #err.cv    <- rmse(DATAj$Y - predict(minfiti, DATAj))
    
           ct6.cv[1:2,f]    = get_sen_spc(est.cv, gtte, thr)
           ct6.resub[1:2,f] = get_sen_spc(est.resub, gttr, thr)
        
           stepbystep.featname[1:6,f] = featnameallshort[c( min1,min2,min3,min4,min5,min6)]
           
        }
 ##############################################################    
        if (ncol(zin)>=7 && nf>=7) {   
           minerr7=1000;
        
           for (i in 1:nfeat) {
                if (!(i %in% c(min1,min2,min3,min4,min5,min6))) {
                    xi <-trin[, featnameallshort[c(min1,min2,min3,min4,min5,min6,i)]]
 
                    DATAi <- data.frame(X=xi, Y=ytr)
                    svmi <- svm(Y ~ ., data=DATAi)
                    #svmiResult <- tune(svm, Y ~ ., data=DATAi, ranges=list(epsilon=seq(0,1,0.1),cost=2^(2:9)))
                    #svmi <- svmiResult$best.model
                    ipred <- predict(svmi, DATAi)
                    ierr <- rmse(DATAi$Y-ipred)
                    if (ierr < minerr7) {
                        minpred <- ipred
                        minerr7 <- ierr 
                        min7 <- i 
                        minfiti <- svmi
                    }
               }
           }

           xj <- tein[, featnameallshort[c(min1, min2, min3, min4, min5, min6, min7)]]
           DATAj <- data.frame(X=xj, Y=yte)
    
           est.resub <- minpred
           est.cv    <- predict(minfiti, DATAj)
           #err.cv    <- rmse(DATAj$Y - predict(minfiti, DATAj))
    
           ct7.cv[1:2,f]    = get_sen_spc(est.cv, gtte, thr)
           ct7.resub[1:2,f] = get_sen_spc(est.resub, gttr, thr)
        
           stepbystep.featname[1:7,f] = featnameallshort[c( min1,min2,min3,min4,min5,min6,min7)]
           
        }
 ##############################################################    
        if (ncol(zin)>=8 && nf>=8) {   
           minerr8=1000;
           for (i in 1:nfeat) {
                if (!(i %in% c(min1,min2,min3,min4,min5,min6,min7))) {
                    xi <-trin[, featnameallshort[c(min1,min2,min3,min4,min5,min6,min7,i)]]
 
                    DATAi <- data.frame(X=xi, Y=ytr)
                    svmi <- svm(Y ~ ., data=DATAi)
                    #svmiResult <- tune(svm, Y ~ ., data=DATAi, ranges=list(epsilon=seq(0,1,0.1),cost=2^(2:9)))
                    #svmi <- svmiResult$best.model
                    ipred <- predict(svmi, DATAi)
                    ierr <- rmse(DATAi$Y-ipred)
                    if (ierr < minerr8) {
                        minpred <- ipred
                        minerr8 <- ierr 
                        min8 <- i 
                        minfiti <- svmi
                    }
               }
           }

           xj <- tein[, featnameallshort[c(min1, min2, min3, min4, min5, min6, min7, min8)]]
           DATAj <- data.frame(X=xj, Y=yte)
    
           est.resub <- minpred
           est.cv    <- predict(minfiti, DATAj)
           #err.cv    <- rmse(DATAj$Y - predict(minfiti, DATAj))
    
           ct8.cv[1:2,f]    = get_sen_spc(est.cv, gtte, thr)
           ct8.resub[1:2,f] = get_sen_spc(est.resub, gttr, thr)
        
           stepbystep.featname[1:8,f] = featnameallshort[c( min1,min2,min3,min4,min5,min6,min7,min8)]
           
        }
 ##############################################################    
        if (ncol(zin)>=9 && nf>=9) {   
           minerr9=1000;
           for (i in 1:nfeat) {
                if (!(i %in% c(min1,min2,min3,min4,min5,min6,min7,min8))) {
                    xi <-trin[, featnameallshort[c(min1,min2,min3,min4,min5,min6,min7,min8,i)]]
 
                    DATAi <- data.frame(X=xi, Y=ytr)
                    svmi <- svm(Y ~ ., data=DATAi)
                    #svmiResult <- tune(svm, Y ~ ., data=DATAi, ranges=list(epsilon=seq(0,1,0.1),cost=2^(2:9)))
                    #svmi <- svmiResult$best.model
                    ipred <- predict(svmi, DATAi)
                    ierr <- rmse(DATAi$Y-ipred)
                    if (ierr < minerr9) {
                        minpred <- ipred
                        minerr9 <- ierr 
                        min9 <- i 
                        minfiti <- svmi
                    }
               }
           }

           xj <- tein[, featnameallshort[c(min1, min2, min3, min4, min5, min6, min7, min8, min9)]]
           DATAj <- data.frame(X=xj, Y=yte)
    
           est.resub <- minpred
           est.cv    <- predict(minfiti, DATAj)
           #err.cv    <- rmse(DATAj$Y - predict(minfiti, DATAj))
    
           ct9.cv[1:2,f]    = get_sen_spc(est.cv, gtte, thr)
           ct9.resub[1:2,f] = get_sen_spc(est.resub, gttr, thr)
        
           stepbystep.featname[1:9,f] = featnameallshort[c( min1,min2,min3,min4,min5,min6,min7,min8,min9)]
           
        }
 ##############################################################    
        if (ncol(zin)>=10 && nf>=10) {   
           minerr10=1000;
           for (i in 1:nfeat) {
                if (!(i %in% c(min1,min2,min3,min4,min5,min6,min7,min8,min9))) {
                    xi <-trin[, featnameallshort[c(min1,min2,min3,min4,min5,min6,min7,min8,min9,i)]]
 
                    DATAi <- data.frame(X=xi, Y=ytr)
                    svmi <- svm(Y ~ ., data=DATAi)
                    #svmiResult <- tune(svm, Y ~ ., data=DATAi, ranges=list(epsilon=seq(0,1,0.1),cost=2^(2:9)))
                    #svmi <- svmiResult$best.model
                    ipred <- predict(svmi, DATAi)
                    ierr <- rmse(DATAi$Y-ipred)
                    if (ierr < minerr10) {
                        minpred <- ipred
                        minerr10 <- ierr 
                        min10 <- i 
                        minfiti <- svmi
                    }
               }
           }

           xj <- tein[, featnameallshort[c(min1, min2, min3, min4, min5, min6, min7, min8, min9, min10)]]
           DATAj <- data.frame(X=xj, Y=yte)
    
           est.resub <- minpred
           est.cv    <- predict(minfiti, DATAj)
           #err.cv    <- rmse(DATAj$Y - predict(minfiti, DATAj))
    
           ct10.cv[1:2,f]    = get_sen_spc(est.cv, gtte, thr)
           ct10.resub[1:2,f] = get_sen_spc(est.resub, gttr, thr)
        
           stepbystep.featname[1:10,f] = featnameallshort[c( min1,min2,min3,min4,min5,min6,min7,min8,min9,min10)]

        } 
        
    }  # nfold
        
        
    avg.ct1.cv <- apply(ct1.cv, 1, mean)
    avg.ct2.cv <- apply(ct2.cv, 1, mean) 
    avg.ct3.cv <- apply(ct3.cv, 1, mean) 
    avg.ct4.cv <- apply(ct4.cv, 1, mean) 
    avg.ct5.cv <- apply(ct5.cv, 1, mean) 
    avg.ct6.cv <- apply(ct6.cv, 1, mean) 
    avg.ct7.cv <- apply(ct7.cv, 1, mean) 
    avg.ct8.cv <- apply(ct8.cv, 1, mean) 
    avg.ct9.cv <- apply(ct9.cv, 1, mean) 
    avg.ct10.cv <- apply(ct10.cv, 1, mean) 
    
    sd.ct1.cv <- apply(ct1.cv, 1, sd) 
    sd.ct2.cv <- apply(ct2.cv, 1, sd) 
    sd.ct3.cv <- apply(ct3.cv, 1, sd) 
    sd.ct4.cv <- apply(ct4.cv, 1, sd) 
    sd.ct5.cv <- apply(ct5.cv, 1, sd) 
    sd.ct6.cv <- apply(ct6.cv, 1, sd) 
    sd.ct7.cv <- apply(ct7.cv, 1, sd) 
    sd.ct8.cv <- apply(ct8.cv, 1, sd) 
    sd.ct9.cv <- apply(ct9.cv, 1, sd) 
    sd.ct10.cv <- apply(ct10.cv, 1, sd) 
   
    avg.ct1.resub <- apply(ct1.resub, 1, mean) 
    avg.ct2.resub <- apply(ct2.resub, 1, mean) 
    avg.ct3.resub <- apply(ct3.resub, 1, mean) 
    avg.ct4.resub <- apply(ct4.resub, 1, mean) 
    avg.ct5.resub <- apply(ct5.resub, 1, mean) 
    avg.ct6.resub <- apply(ct6.resub, 1, mean) 
    avg.ct7.resub <- apply(ct7.resub, 1, mean) 
    avg.ct8.resub <- apply(ct8.resub, 1, mean) 
    avg.ct9.resub <- apply(ct9.resub, 1, mean) 
    avg.ct10.resub <- apply(ct10.resub, 1, mean) 
        
    sd.ct1.resub <- apply(ct1.resub, 1, sd) 
    sd.ct2.resub <- apply(ct2.resub, 1, sd) 
    sd.ct3.resub <- apply(ct3.resub, 1, sd) 
    sd.ct4.resub <- apply(ct4.resub, 1, sd) 
    sd.ct5.resub <- apply(ct5.resub, 1, sd) 
    sd.ct6.resub <- apply(ct6.resub, 1, sd) 
    sd.ct7.resub <- apply(ct7.resub, 1, sd) 
    sd.ct8.resub <- apply(ct8.resub, 1, sd) 
    sd.ct9.resub <- apply(ct9.resub, 1, sd) 
    sd.ct10.resub <- apply(ct10.resub, 1, sd) 
        
        
    stepbystep.feat.mean.resub <- cbind( avg.ct1.resub , avg.ct2.resub , avg.ct3.resub , avg.ct4.resub , avg.ct5.resub , avg.ct6.resub , avg.ct7.resub , avg.ct8.resub , avg.ct9.resub , avg.ct10.resub)
    stepbystep.feat.mean.cv <- cbind(avg.ct1.cv , avg.ct2.cv, avg.ct3.cv, avg.ct4.cv, avg.ct5.cv, avg.ct6.cv, avg.ct7.cv, avg.ct8.cv, avg.ct9.cv, avg.ct10.cv) 
        
    stepbystep.feat.sd.resub <- cbind(sd.ct1.resub , sd.ct2.resub , sd.ct3.resub , sd.ct4.resub , sd.ct5.resub , sd.ct6.resub , sd.ct7.resub , sd.ct8.resub , sd.ct9.resub , sd.ct10.resub)
    stepbystep.feat.sd.cv <- cbind(sd.ct1.cv , sd.ct2.cv, sd.ct3.cv, sd.ct4.cv, sd.ct5.cv, sd.ct6.cv, sd.ct7.cv, sd.ct8.cv, sd.ct9.cv, sd.ct10.cv)
        
    spcsen.cv <- matrix(0.0,nrow=2,ncol=nf)
    spcsen.resub <- matrix(0.0,nrow=2,ncol=nf)
    sd.spcsen.cv <- matrix(0.0,nrow=2,ncol=nf)
    sd.spcsen.resub <- matrix(0.0,nrow=2,ncol=nf)
 if (nf>=1) {    
    spcsen.cv[,1]    <- c(mean(1-ct1.cv[2,1:nfold]), mean(ct1.cv[1,1:nfold]))
    spcsen.resub[,1]    <- c(mean(1-ct1.resub[2,1:nfold]), mean(ct1.resub[1,1:nfold]))
    sd.spcsen.cv[,1]    <- c(sd(1-ct1.cv[2,1:nfold]), sd(ct1.cv[1,1:nfold]))
    sd.spcsen.resub[,1]    <- c(sd(1-ct1.resub[2,1:nfold]), sd(ct1.resub[1,1:nfold]))
    } 
 if (nf>=2) {    
    spcsen.cv[,2]    <- c(mean(1-ct2.cv[2,1:nfold]), mean(ct2.cv[1,1:nfold]))
    spcsen.resub[,2]    <- c(mean(1-ct2.resub[2,1:nfold]), mean(ct2.resub[1,1:nfold]))
    sd.spcsen.cv[,2]    <- c(sd(1-ct2.cv[2,1:nfold]), sd(ct2.cv[1,1:nfold]))
    sd.spcsen.resub[,2]    <- c(sd(1-ct2.resub[2,1:nfold]), sd(ct2.resub[1,1:nfold]))
    } 
 if (nf>=3) {    
    spcsen.cv[,3]    <- c(mean(1-ct3.cv[2,1:nfold]), mean(ct3.cv[1,1:nfold]))
    spcsen.resub[,3]    <- c(mean(1-ct3.resub[2,1:nfold]), mean(ct3.resub[1,1:nfold]))
    sd.spcsen.cv[,3]    <- c(sd(1-ct3.cv[2,1:nfold]), sd(ct3.cv[1,1:nfold]))
    sd.spcsen.resub[,3]    <- c(sd(1-ct3.resub[2,1:nfold]), sd(ct3.resub[1,1:nfold]))
    } 
 if (nf>=4) {    
    spcsen.cv[,4]    <- c(mean(1-ct4.cv[2,1:nfold]), mean(ct4.cv[1,1:nfold]))
    spcsen.resub[,4]    <- c(mean(1-ct4.resub[2,1:nfold]), mean(ct4.resub[1,1:nfold]))
    sd.spcsen.cv[,4]    <- c(sd(1-ct4.cv[2,1:nfold]), sd(ct4.cv[1,1:nfold]))
    sd.spcsen.resub[,4]    <- c(sd(1-ct4.resub[2,1:nfold]), sd(ct4.resub[1,1:nfold]))
    } 
 if (nf>=5) {    
    spcsen.cv[,5]    <- c(mean(1-ct5.cv[2,1:nfold]), mean(ct5.cv[1,1:nfold]))
    spcsen.resub[,5]    <- c(mean(1-ct5.resub[2,1:nfold]), mean(ct5.resub[1,1:nfold]))
    sd.spcsen.cv[,5]    <- c(sd(1-ct5.cv[2,1:nfold]), sd(ct5.cv[1,1:nfold]))
    sd.spcsen.resub[,5]    <- c(sd(1-ct5.resub[2,1:nfold]), sd(ct5.resub[1,1:nfold]))
    } 
 if (nf>=6) {    
    spcsen.cv[,6]    <- c(mean(1-ct6.cv[2,1:nfold]), mean(ct6.cv[1,1:nfold]))
    spcsen.resub[,6]    <- c(mean(1-ct6.resub[2,1:nfold]), mean(ct6.resub[1,1:nfold]))
    sd.spcsen.cv[,6]    <- c(sd(1-ct6.cv[2,1:nfold]), sd(ct6.cv[1,1:nfold]))
    sd.spcsen.resub[,6]    <- c(sd(1-ct6.resub[2,1:nfold]), sd(ct6.resub[1,1:nfold]))
    } 
 if (nf>=7) {    
    spcsen.cv[,7]    <- c(mean(1-ct7.cv[2,1:nfold]), mean(ct7.cv[1,1:nfold]))
    spcsen.resub[,7]    <- c(mean(1-ct7.resub[2,1:nfold]), mean(ct7.resub[1,1:nfold]))
    sd.spcsen.cv[,7]    <- c(sd(1-ct7.cv[2,1:nfold]), sd(ct7.cv[1,1:nfold]))
    sd.spcsen.resub[,7]    <- c(sd(1-ct7.resub[2,1:nfold]), sd(ct7.resub[1,1:nfold]))
    } 
 if (nf>=8) {    
    spcsen.cv[,8]    <- c(mean(1-ct8.cv[2,1:nfold]), mean(ct8.cv[1,1:nfold]))
    spcsen.resub[,8]    <- c(mean(1-ct8.resub[2,1:nfold]), mean(ct8.resub[1,1:nfold]))
    sd.spcsen.cv[,8]    <- c(sd(1-ct8.cv[2,1:nfold]), sd(ct8.cv[1,1:nfold]))
    sd.spcsen.resub[,8]    <- c(sd(1-ct8.resub[2,1:nfold]), sd(ct8.resub[1,1:nfold]))
    } 
 if (nf>=9) {    
    spcsen.cv[,9]    <- c(mean(1-ct9.cv[2,1:nfold]), mean(ct9.cv[1,1:nfold]))
    spcsen.resub[,9]    <- c(mean(1-ct9.resub[2,1:nfold]), mean(ct9.resub[1,1:nfold]))
    sd.spcsen.cv[,9]    <- c(sd(1-ct9.cv[2,1:nfold]), sd(ct9.cv[1,1:nfold]))
    sd.spcsen.resub[,9]    <- c(sd(1-ct9.resub[2,1:nfold]), sd(ct9.resub[1,1:nfold]))
    } 
 if (nf>=10) {    
    spcsen.cv[,10]    <- c(mean(1-ct10.cv[2,1:nfold]), mean(ct10.cv[1,1:nfold]))
    spcsen.resub[,10]    <- c(mean(1-ct10.resub[2,1:nfold]), mean(ct10.resub[1,1:nfold]))
    sd.spcsen.cv[,10]    <- c(sd(1-ct10.cv[2,1:nfold]), sd(ct10.cv[1,1:nfold]))
    sd.spcsen.resub[,10]    <- c(sd(1-ct10.resub[2,1:nfold]), sd(ct10.resub[1,1:nfold]))
    } 

if (debug) {
    print("spcsen.cv")
    print(spcsen.cv)
    print("sd.spcsen.cv")
    print(sd.spcsen.cv)
    print("spcsen.resub")
    print(spcsen.resub)
    print("sd.spcsen.resub")
    print(sd.spcsen.resub)
}
    pvalcindex.cv <- matrix(0.0,nrow=2,ncol=nf)
    pvalcindex.resub <- matrix(0.0,nrow=2,ncol=nf)
    sd.pvalcindex.cv <- matrix(0.0,nrow=2,ncol=nf)
    sd.pvalcindex.resub <- matrix(0.0,nrow=2,ncol=nf)
    
if (debug) {
    print("stepbystep.featname")
    print(stepbystep.featname)
}
    myout <- list()
    myout$stepbystep.featname <- stepbystep.featname[1:nf,]
    myout$stepbystep.feat.mean.resub <- stepbystep.feat.mean.resub
    myout$stepbystep.feat.sd.resub <- stepbystep.feat.sd.resub
    myout$stepbystep.feat.mean.cv <- stepbystep.feat.mean.cv
    myout$stepbystep.feat.sd.cv <- stepbystep.feat.sd.cv
    myout$spcsen.cv <- spcsen.cv
    myout$spcsen.resub <- spcsen.resub
    myout$sd.spcsen.cv <- sd.spcsen.cv
    myout$sd.spcsen.resub <- sd.spcsen.resub
    myout$pvalcindex.cv <- pvalcindex.cv
    myout$pvalcindex.resub <- pvalcindex.resub
    myout$sd.pvalcindex.cv <- sd.pvalcindex.cv
    myout$sd.pvalcindex.resub <- sd.pvalcindex.resub
    myout$nf <- nf
    
if (flagPlot) {print("	making plots ") }
    pdf(rocpdffn, width=8, height=8)
    par(mfrow=c(1,2))
    plot(1:nf, spcsen.cv[1,1:nf], pch=16, col="red", main="false positive rate \n (1 - specificity)", ylab="false positive rate", xlab="number of features", ylim=c(0,0.5), type="o", cex=0.8)
    arrows(1:nf, y0=spcsen.cv[1,1:nf]-sd.spcsen.cv[1,1:nf], y1=spcsen.cv[1,1:nf]+sd.spcsen.cv[1,1:nf], code=3, angle=90, length=.05, col="red")
    lines(1:nf, spcsen.resub[1,1:nf], pch=16, col="blue", type="o")
    arrows(1:nf, y0=spcsen.resub[1,1:nf]-sd.spcsen.resub[1,1:nf], y1=spcsen.resub[1,1:nf]+sd.spcsen.resub[1,1:nf], code=3, angle=90, length=.05, col="blue")
    grid()
    legend("topleft", legend=c(paste(nfold, " fold cv"), "resub"), pch=c(16,16), col=c("red", "blue"))
    plot(1:nf, spcsen.cv[2,1:nf], pch=16, col="red", main="sensitivity", ylab="sensitivity", xlab="number of features", ylim=c(0.5, 1), type="o", cex=0.8)
    arrows(1:nf, y0=spcsen.cv[2,1:nf]-sd.spcsen.cv[2,1:nf], y1=spcsen.cv[2,1:nf]+sd.spcsen.cv[2,1:nf], code=3, angle=90, length=.05, col="red")
    lines(1:nf, spcsen.resub[2,1:nf], pch=16, col="blue", type="o")
    arrows(1:nf, y0=spcsen.resub[2,1:nf]-sd.spcsen.resub[2,1:nf], y1=spcsen.resub[2,1:nf]+sd.spcsen.resub[2,1:nf], code=3, angle=90, length=.05, col="blue")
    grid()
    legend("bottomright", legend=c(paste(nfold, "fold cv"), "resub"), pch=c(16,16), col=c("red", "blue"))

    plot(1:nf, pvalcindex.cv[1,1:nf], pch=16, col="red", main="single cluster : pvalue", ylab="-log10(pval)", xlab="number of features", ylim=c(0,max(c(pvalcindex.cv[1,1:nf]))), type="o", cex=0.8)
    arrows(1:nf, y0=pvalcindex.cv[1,1:nf]-sd.pvalcindex.cv[1,1:nf], y1=pvalcindex.cv[1,1:nf]+sd.pvalcindex.cv[1,1:nf], code=3, angle=90, length=.05, col="red")
    lines(1:nf, pvalcindex.resub[1,1:nf], pch=16, col="blue", type="o")
    arrows(1:nf, y0=pvalcindex.resub[1,1:nf]-sd.pvalcindex.resub[1,1:nf], y1=pvalcindex.resub[1,1:nf]+sd.pvalcindex.resub[1,1:nf], code=3, angle=90, length=.05, col="blue")
    grid()
    legend("topleft", legend=c(paste(nfold, " fold cv"), "resub"), pch=c(16,16), col=c("red", "blue"))
    plot(1:nf, pvalcindex.cv[2,1:nf], pch=16, col="red", main="cluster_index", ylab="cindex", xlab="number of features", ylim=c(0, 1), type="o", cex=0.8)
    arrows(1:nf, y0=pvalcindex.cv[2,1:nf]-sd.pvalcindex.cv[2,1:nf], y1=pvalcindex.cv[2,1:nf]+sd.pvalcindex.cv[2,1:nf], code=3, angle=90, length=.05, col="red")
    lines(1:nf, pvalcindex.resub[2,1:nf], pch=16, col="blue", type="o")
    arrows(1:nf, y0=pvalcindex.resub[2,1:nf]-sd.pvalcindex.resub[2,1:nf], y1=pvalcindex.resub[2,1:nf]+sd.pvalcindex.resub[2,1:nf], code=3, angle=90, length=.05, col="blue")
    grid()
    legend("bottomright", legend=c(paste(nfold, "fold cv"), "resub"), pch=c(16,16), col=c("red", "blue"))
    dev.off()
#
if (0) {
    print("in mystepbystep") 
    print(sss.str)
    print("dim(out[[sss.str]]$spcsen.cv")
    print(dim(out[[sss.str]]$spcsen.cv))
    print("dim(out[[sss.str]]$spcsen.resub")
    print(dim(out[[sss.str]]$spcsen.resub))
}


if (debug) {print("	returning from mystepbystepcv") }
myout


}


mystepbystepcv <- function(gt, y, zin, featnameallshort, nfeat, nfold, thr, rocpdffn, debug=0, flagPlot=TRUE) {

    if (debug) {
        gt <- orig.class
        y <- orig.class
        zin <- orig.ZSel35
        featnameallshort <- orig_feat_pcloading
        nfeat <-38
        nfold <- 4
        thr <- 1.5 
        rocpdffn <- 'ck2test.11.pdf'
    }

    if (nfeat <= 10) {
        nf=nfeat
    } else {
        nf=10
    }
#print(paste("number of max feature to choose ", nf))
#print(paste("number of columns in data ", ncol(zin)))

    y2.idx <- which(y==2)
    y2.n <- length(y2.idx)
    
    y1.idx <- which(y==1)
    y1.n <- length(y1.idx)
    
    err1.cv <- err2.cv <- err3.cv <- err4.cv <- err5.cv <- err6.cv <- err7.cv <- err8.cv <- err9.cv <- err10.cv <- rep(0,nfold+1) 
    err1.resub <- err2.resub <- err3.resub <- err4.resub <- err5.resub <- err6.resub <- err7.resub <- err8.resub <- err9.resub <- err10.resub <- rep(0,nfold+1) 
    ct1.cv <- ct2.cv <- ct3.cv <- ct4.cv <- ct5.cv <- ct6.cv <- ct7.cv <- ct8.cv <- ct9.cv <- ct10.cv <- matrix(0,ncol=nfold+1,nrow=2)
    ct1.resub <- ct2.resub <- ct3.resub <- ct4.resub <- ct5.resub <- ct6.resub <- ct7.resub <- ct8.resub <- ct9.resub <- ct10.resub <- matrix(0,ncol=nfold+1,nrow=2) 
    sig1.cv <- sig2.cv <- sig3.cv <- sig4.cv <- sig5.cv <- sig6.cv <- sig7.cv <- sig8.cv <- sig9.cv <- sig10.cv <- matrix(0,ncol=nfold+1,nrow=2)
    sig1.resub <- sig2.resub <- sig3.resub <- sig4.resub <- sig5.resub <- sig6.resub <- sig7.resub <- sig8.resub <- sig9.resub <- sig10.resub <- matrix(0,ncol=nfold+1,nrow=2) 
        
    stepbystep.featname = matrix("hello", 10, nfold+1)
    nsample = nrow(zin)
    
    mymod1 = 1:y1.n %% nfold
    mymod2 = 1:y2.n %% nfold
    out <- list()

    for (f in 1:(nfold+1)) {
        #print(paste("fold ",f))
        #idxtr = c(sample(y2.idx,f.samplesize), y1.idx) 
        #dxte = c(setdiff(1:nsample,idxtr), y1.idx)
        if (f > nfold) {
           idxtr = c(y2.idx, y1.idx)
           idxte = idxtr
        } else {
           mask1 <- mymod1 == (f-1)
           mask2 <- mymod2 == (f-1)
           idx1 <- idx2 <- c()
           if (sum(mask1) > 0) idx1 <- y1.idx[which(mask1)]
           if (sum(mask2) > 0) idx2 <- y2.idx[which(mask2)]
           idxte = c(idx1, idx2)
           idxtr = setdiff(1:nsample,idxte)
        }

        trin = zin[idxtr,]
        tein = zin[idxte,]
        iintr  = matrix(1,length(idxtr),1)
        iinte  = matrix(1,length(idxte),1)
        ytr = y[idxtr]
        yte = y[idxte]
        gttr = gt[idxtr]
        gtte = gt[idxte]
#   print("c(length(idxte), length(idxtr))")        
#   print(c(length(idxte), length(idxtr)))        
   if (length(idxte)>1 && length(idxtr)>1) {

        if (debug) { print(0) }
        if (ncol(zin)>=1 && nf>=1) {
           minerr1=1
           for (i in 1:nfeat) {
                xi = trin[, featnameallshort[i]]
                fiti = lm(ytr ~ xi)
                ierr = mean(abs(fiti$residuals))
                if (ierr < minerr1) {
                    minerr1 = ierr 
                    min1 = i 
                    minfiti = fiti
                }
           }
           if (1) {
              err1.resub[f] = minerr1 ;
              est.resub = minfiti$coefficients %*%  t(cbind(iintr, trin[, featnameallshort[c(min1)]]))
              est.cv    = minfiti$coefficients %*%  t(cbind(iinte, tein[, featnameallshort[c(min1)]]))
              err1.cv[f]    = mean(abs(est.cv - gtte))
              ct1.cv[1:2,f]    = get_sen_spc(est.cv, gtte, thr)
              ct1.resub[1:2,f] = get_sen_spc(est.resub, gttr, thr)
           } else {
              err1.resub[f] = get_resub_err(fiti, ytr) ;
              err1.cv[f]    = get_cv_err(fiti, tein[, featnameallshort[c(min1)]], yte)
           }
        
          
           if (debug) print("1")
           if (0) {
               tmp.resub <- sigclust(trin[, featnameallshort[c(min1,min1)]], 100, labflag=1, label=as.vector(est.resub>thr), icovest=2)
               tmp.cv <- sigclust(tein[, featnameallshort[c(min1,min1)]], 100, labflag=1, label=as.vector(est.cv>thr), icovest=2)
               sig1.resub[1:2,f] <- c(tmp.resub@xcindex, -log10(tmp.resub@pvalnorm))
               sig1.cv[1:2,f] <- c(tmp.cv@xcindex, -log10(tmp.cv@pvalnorm))

               sig1.resub[1:2,f] <- c(0,0)
               sig1.cv[1:2,f] <- c(0,0)
           }
           stepbystep.featname[1,f] = featnameallshort[c( min1)]
     
        }
        if (debug) print(1)
        
        if (ncol(zin)>=2 && nf>=2) {   
           minerr2=1;
           x1 = trin[, featnameallshort[min1]]
           for (i in 1:nfeat) {
                if (!(i %in% min1)) {
                    xi = trin[, featnameallshort[i]]
                    fiti = lm(ytr ~ x1 + xi)
                    ierr = mean(abs(fiti$residuals))
                    if (ierr < minerr2) {
                        minerr2 = ierr 
                        min2 = i 
                    minfiti = fiti
                    }
                }
           }
           if (1) {
              err2.resub[f] = minerr2 ;
              est.resub = minfiti$coefficients %*%  t(cbind(iintr, trin[, featnameallshort[c(min1,min2)]]))
              est.cv    = minfiti$coefficients %*%  t(cbind(iinte, tein[, featnameallshort[c(min1,min2)]]))
              err2.cv[f]    = mean(abs(est.cv - gtte))
              ct2.cv[1:2,f]    = get_sen_spc(est.cv, gtte, thr)
              ct2.resub[1:2,f] = get_sen_spc(est.resub, gttr, thr)
           } else {
              err2.resub[f] = get_resub_err(fiti, ytr) ;
              err2.cv[f]    = get_cv_err(fiti, tein[, featnameallshort[c(min1,min2)]], yte)
           } 
           
           if (debug) print("2")

           if (0) {
              tmp.resub <- sigclust(trin[, featnameallshort[c(min1,min2)]], 100, labflag=1, label=as.vector(est.resub>thr), icovest=2)
              sig2.resub[1:2,f] <- c(tmp.resub@xcindex, -log10(tmp.resub@pvalnorm))
              tmp.cv <- sigclust(tein[, featnameallshort[c(min1,min2)]], 100, labflag=1, label=as.vector(est.cv>thr), icovest=2)
              sig2.cv[1:2,f] <- c(tmp.cv@xcindex, -log10(tmp.cv@pvalnorm))
           }
           stepbystep.featname[1:2,f] = featnameallshort[c( min1,min2)]
           
        }
        if (ncol(zin)>=3 && nf>=3) {   
           minerr3=1;
           x1 = trin[, featnameallshort[min1]]
           x2 = trin[, featnameallshort[min2]]
           for (i in 1:nfeat) {
                if (!(i %in% c(min1,min2))) {
                    xi = trin[, featnameallshort[i]]
                    fiti = lm(ytr ~ x1 + x2 + xi)
                    ierr = mean(abs(fiti$residuals))
                    if (ierr < minerr3) {
                        minerr3 = ierr 
                        min3 = i 
                    minfiti = fiti
                    }
                }
           }
           if (1) {
           err3.resub[f] = minerr3 ;
           est.resub = minfiti$coefficients %*%  t(cbind(iintr, trin[, featnameallshort[c(min1,min2,min3)]]))
           est.cv    = minfiti$coefficients %*%  t(cbind(iinte, tein[, featnameallshort[c(min1,min2,min3)]]))
           err3.cv[f]    = mean(abs(est.cv - gtte))
           ct3.cv[1:2,f]    = get_sen_spc(est.cv, gtte, thr)
           ct3.resub[1:2,f] = get_sen_spc(est.resub, gttr, thr)
           } else {
           err3.resub[f] = get_resub_err(fiti, ytr) ;
           err3.cv[f]    = get_cv_err(fiti, tein[, featnameallshort[c(min1,min2,min3)]], yte)
           } 

           if (debug) print("3")
           if (0) {
              tmp.resub <- sigclust(trin[, featnameallshort[c(min1,min2,min3)]], 100, labflag=1, label=as.vector(est.resub>thr), icovest=2)
              sig3.resub[1:2,f] <- c(tmp.resub@xcindex, -log10(tmp.resub@pvalnorm))
              tmp.cv <- sigclust(tein[, featnameallshort[c(min1,min2,min3)]], 100, labflag=1, label=as.vector(est.cv>thr), icovest=2)
              sig3.cv[1:2,f] <- c(tmp.cv@xcindex, -log10(tmp.cv@pvalnorm))
           }
           stepbystep.featname[1:3,f] = featnameallshort[c( min1,min2,min3)]
           
        }
        if (ncol(zin)>=4 && nf>=4) {   
           minerr4=1;
           x1 = trin[, featnameallshort[min1]]
           x2 = trin[, featnameallshort[min2]]
           x3 = trin[, featnameallshort[min3]]
           for (i in 1:nfeat) {
                if (!(i %in% c(min1,min2,min3))) {
                    xi = trin[, featnameallshort[i]]
                    fiti = lm(ytr ~ x1 + x2 + x3 + xi)
                    ierr = mean(abs(fiti$residuals))
                    if (ierr < minerr4) {
                        minerr4 = ierr 
                        min4 = i 
                    minfiti = fiti
                    }
                }
           }
           if (1) {
           err4.resub[f] = minerr4 ;
           est.resub = minfiti$coefficients %*%  t(cbind(iintr, trin[, featnameallshort[c(min1,min2,min3,min4)]]))
           est.cv    = minfiti$coefficients %*%  t(cbind(iinte, tein[, featnameallshort[c(min1,min2,min3,min4)]]))
           err4.cv[f]    = mean(abs(est.cv - gtte))
           ct4.cv[1:2,f]    = get_sen_spc(est.cv, gtte, thr)
           ct4.resub[1:2,f] = get_sen_spc(est.resub, gttr, thr)
           } else {
           err4.resub[f] = get_resub_err(fiti, ytr) ;
           err4.cv[f]    = get_cv_err(fiti, tein[, featnameallshort[c(min1,min2,min3,min4)]], yte)
           }

if (debug) print("4")
if (0) {
           tmp.resub <- sigclust(trin[, featnameallshort[c(min1,min2,min3,min4)]], 100, labflag=1, label=as.vector(est.resub>thr), icovest=2)
           sig4.resub[1:2,f] <- c(tmp.resub@xcindex, -log10(tmp.resub@pvalnorm))
           tmp.cv <- sigclust(tein[, featnameallshort[c(min1,min2,min3,min4)]], 100, labflag=1, label=as.vector(est.cv>thr), icovest=2)
           sig4.cv[1:2,f] <- c(tmp.cv@xcindex, -log10(tmp.cv@pvalnorm))
}
           stepbystep.featname[1:4,f] = featnameallshort[c( min1,min2,min3,min4)]
           
        }
        if (ncol(zin)>=5 && nf>=5) {   
           minerr5=1;
           x1 = trin[, featnameallshort[min1]]
           x2 = trin[, featnameallshort[min2]]
           x3 = trin[, featnameallshort[min3]]
           x4 = trin[, featnameallshort[min4]]
           for (i in 1:nfeat) {
                if (!(i %in% c(min1,min2,min3,min4))) {
                    xi = trin[, featnameallshort[i]]
                    fiti = lm(ytr ~ x1 + x2 + x3 + x4 + xi)
                    ierr = mean(abs(fiti$residuals))
                    if (ierr < minerr5) {
                        minerr5 = ierr 
                        min5 = i 
                    minfiti = fiti
                    }
                }
           }
           if (1) {
           err5.resub[f] = minerr5 ;
           est.resub = minfiti$coefficients %*%  t(cbind(iintr, trin[, featnameallshort[c(min1,min2,min3,min4,min5)]]))
           est.cv    = minfiti$coefficients %*%  t(cbind(iinte, tein[, featnameallshort[c(min1,min2,min3,min4,min5)]]))
           err5.cv[f]    = mean(abs(est.cv - gtte))
           ct5.cv[1:2,f]    = get_sen_spc(est.cv, gtte, thr)
           ct5.resub[1:2,f] = get_sen_spc(est.resub, gttr, thr)
           } else {
           err5.resub[f] = get_resub_err(fiti, ytr) ;
           err5.cv[f]    = get_cv_err(fiti, tein[, featnameallshort[c(min1,min2,min3,min4,min5)]], yte)
           } 
           
if (debug) print("5")
if (0) {
           tmp.resub <- sigclust(trin[, featnameallshort[c(min1,min2,min3,min4,min5)]], 100, labflag=1, label=as.vector(est.resub>thr), icovest=2)
           sig5.resub[1:2,f] <- c(tmp.resub@xcindex, -log10(tmp.resub@pvalnorm))
           tmp.cv <- sigclust(tein[, featnameallshort[c(min1,min2,min3,min4,min5)]], 100, labflag=1, label=as.vector(est.cv>thr), icovest=2)
           sig5.cv[1:2,f] <- c(tmp.cv@xcindex, -log10(tmp.cv@pvalnorm))
}
           stepbystep.featname[1:5,f] = featnameallshort[c( min1,min2,min3,min4,min5)]

        }
        if (ncol(zin)>=6 && nf>=6) {   
           minerr6=1;
           x1 = trin[, featnameallshort[min1]]
           x2 = trin[, featnameallshort[min2]]
           x3 = trin[, featnameallshort[min3]]
           x4 = trin[, featnameallshort[min4]]
           x5 = trin[, featnameallshort[min5]]
           for (i in 1:nfeat) {
                if (!(i %in% c(min1,min2,min3,min4,min5))) {
                    xi = trin[, featnameallshort[i]]
                    fiti = lm(ytr ~ x1 + x2 + x3 + x4 + x5 + xi)
                    ierr = mean(abs(fiti$residuals))
                    if (ierr < minerr6) {
                        minerr6 = ierr 
                        min6 = i 
                    minfiti = fiti
                    }
                }
           }
           if (1) {
           err6.resub[f] = minerr6 ;
           est.resub = minfiti$coefficients %*%  t(cbind(iintr, trin[, featnameallshort[c(min1,min2,min3,min4,min5,min6)]]))
           est.cv    = minfiti$coefficients %*%  t(cbind(iinte, tein[, featnameallshort[c(min1,min2,min3,min4,min5,min6)]]))
           err6.cv[f]    = mean(abs(est.cv - gtte))
           ct6.cv[1:2,f]    = get_sen_spc(est.cv, gtte, thr)
           ct6.resub[1:2,f] = get_sen_spc(est.resub, gttr, thr)
           } else {
           err6.resub[f] = get_resub_err(fiti, ytr) ;
           err6.cv[f]    = get_cv_err(fiti, tein[, featnameallshort[c(min1,min2,min3,min4,min5,min6)]], yte)
           }
           
if (debug) print("6")
if (0) {
           tmp.resub <- sigclust(trin[, featnameallshort[c(min1,min2,min3,min4,min5,min6)]], 100, labflag=1, label=as.vector(est.resub>thr), icovest=2)
           sig6.resub[1:2,f] <- c(tmp.resub@xcindex, -log10(tmp.resub@pvalnorm))
           tmp.cv <- sigclust(tein[, featnameallshort[c(min1,min2,min3,min4,min5,min6)]], 100, labflag=1, label=as.vector(est.cv>thr), icovest=2)
           sig6.cv[1:2,f] <- c(tmp.cv@xcindex, -log10(tmp.cv@pvalnorm))
}
           stepbystep.featname[1:6,f] = featnameallshort[c( min1,min2,min3,min4,min5,min6)]
           
        }
        if (ncol(zin)>=7 && nf>=7) {   
           minerr7=1;
           x1 = trin[, featnameallshort[min1]]
           x2 = trin[, featnameallshort[min2]]
           x3 = trin[, featnameallshort[min3]]
           x4 = trin[, featnameallshort[min4]]
           x5 = trin[, featnameallshort[min5]]
           x6 = trin[, featnameallshort[min6]]
        
           for (i in 1:nfeat) {
                if (!(i %in% c(min1,min2,min3,min4,min5,min6))) {
                    xi = trin[, featnameallshort[i]]
                    fiti = lm(ytr ~ x1 + x2 + x3 + x4 + x5 + x6 + xi)
                    ierr = mean(abs(fiti$residuals))
                    if (ierr < minerr7) {
                        minerr7 = ierr 
                        min7 = i 
                    minfiti = fiti
                    }
                }
           }
           if (1) {
           err7.resub[f] = minerr7 ;
           est.resub = minfiti$coefficients %*%  t(cbind(iintr, trin[, featnameallshort[c(min1,min2,min3,min4,min5,min6,min7)]]))
           est.cv    = minfiti$coefficients %*%  t(cbind(iinte, tein[, featnameallshort[c(min1,min2,min3,min4,min5,min6,min7)]]))
           err7.cv[f]    = mean(abs(est.cv - gtte))
           ct7.cv[1:2,f]    = get_sen_spc(est.cv, gtte, thr)
           ct7.resub[1:2,f] = get_sen_spc(est.resub, gttr, thr)
           } else {
           err7.resub[f] = get_resub_err(fiti, ytr) ;
           err7.cv[f]    = get_cv_err(fiti, tein[, featnameallshort[c(min1,min2,min3,min4,min5,min6,min7)]], yte)
           }
           
if (debug) print("7")
if (0) {
           tmp.resub <- sigclust(trin[, featnameallshort[c(min1,min2,min3,min4,min5,min6,min7)]], 100, labflag=1, label=as.vector(est.resub>thr), icovest=2)
           sig7.resub[1:2,f] <- c(tmp.resub@xcindex, -log10(tmp.resub@pvalnorm))
           tmp.cv <- sigclust(tein[, featnameallshort[c(min1,min2,min3,min4,min5,min6,min7)]], 100, labflag=1, label=as.vector(est.cv>thr), icovest=2)
           sig7.cv[1:2,f] <- c(tmp.cv@xcindex, -log10(tmp.cv@pvalnorm))
}
           stepbystep.featname[1:7,f] = featnameallshort[c( min1,min2,min3,min4,min5,min6,min7)]
           
           
        }
        if (ncol(zin)>=8 && nf>=8) {   
           minerr8=1;
           x1 = trin[, featnameallshort[min1]]
           x2 = trin[, featnameallshort[min2]]
           x3 = trin[, featnameallshort[min3]]
           x4 = trin[, featnameallshort[min4]]
           x5 = trin[, featnameallshort[min5]]
           x6 = trin[, featnameallshort[min6]]
           x7 = trin[, featnameallshort[min7]]
           for (i in 1:nfeat) {
                if (!(i %in% c(min1,min2,min3,min4,min5,min6,min7))) {
                    xi = trin[, featnameallshort[i]]
                    fiti = lm(ytr ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + xi)
                    ierr = mean(abs(fiti$residuals))
                    if (ierr < minerr8) {
                        minerr8 = ierr 
                        min8 = i 
                    minfiti = fiti
                    }
                }
           }
           if (1) {
           err8.resub[f] = minerr8 ;
           est.resub = minfiti$coefficients %*%  t(cbind(iintr, trin[, featnameallshort[c(min1,min2,min3,min4,min5,min6,min7,min8)]]))
           est.cv    = minfiti$coefficients %*%  t(cbind(iinte, tein[, featnameallshort[c(min1,min2,min3,min4,min5,min6,min7,min8)]]))
           err8.cv[f]    = mean(abs(est.cv - gtte))
           ct8.cv[1:2,f]    = get_sen_spc(est.cv, gtte, thr)
           ct8.resub[1:2,f] = get_sen_spc(est.resub, gttr, thr)
           } else {
           err8.resub[f] = get_resub_err(fiti, ytr) ;
           err8.cv[f]    = get_cv_err(fiti, tein[, featnameallshort[c(min1,min2,min3,min4,min5,min6,min7,min8)]], yte)
           } 
           
if (debug) print("8")
if (0) {
           tmp.resub <- sigclust(trin[, featnameallshort[c(min1,min2,min3,min4,min5,min6,min7,min8)]], 100, labflag=1, label=as.vector(est.resub>thr), icovest=2)
           sig8.resub[1:2,f] <- c(tmp.resub@xcindex, -log10(tmp.resub@pvalnorm))
           tmp.cv <- sigclust(tein[, featnameallshort[c(min1,min2,min3,min4,min5,min6,min7,min8)]], 100, labflag=1, label=as.vector(est.cv>thr), icovest=2)
           sig8.cv[1:2,f] <- c(tmp.cv@xcindex, -log10(tmp.cv@pvalnorm))
}
           stepbystep.featname[1:8,f] = featnameallshort[c( min1,min2,min3,min4,min5,min6,min7,min8)]
           
        }
        if (ncol(zin)>=9 && nf>=9) {   
           minerr9=1;
           x1 = trin[, featnameallshort[min1]]
           x2 = trin[, featnameallshort[min2]]
           x3 = trin[, featnameallshort[min3]]
           x4 = trin[, featnameallshort[min4]]
           x5 = trin[, featnameallshort[min5]]
           x6 = trin[, featnameallshort[min6]]
           x7 = trin[, featnameallshort[min7]]
           x8 = trin[, featnameallshort[min8]]
           for (i in 1:nfeat) {
                if (!(i %in% c(min1,min2,min3,min4,min5,min6,min7,min8))) {
                    xi = trin[, featnameallshort[i]]
                    fiti = lm(ytr ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + xi)
                    ierr = mean(abs(fiti$residuals))
                    if (ierr < minerr9) {
                        minerr9 = ierr 
                        min9 = i 
                    minfiti = fiti
                    }
                }
           }
           if (1) {
           err9.resub[f] = minerr9 ;
           est.resub = minfiti$coefficients %*%  t(cbind(iintr, trin[, featnameallshort[c(min1,min2,min3,min4,min5,min6,min7,min8,min9)]]))
           est.cv    = minfiti$coefficients %*%  t(cbind(iinte, tein[, featnameallshort[c(min1,min2,min3,min4,min5,min6,min7,min8,min9)]]))
           err9.cv[f]    = mean(abs(est.cv - gtte))
           ct9.cv[1:2,f]    = get_sen_spc(est.cv, gtte, thr)
           ct9.resub[1:2,f] = get_sen_spc(est.resub, gttr, thr)
           } else {
           err9.resub[f] = get_resub_err(fiti, ytr) ;
           err9.cv[f]    = get_cv_err(fiti, tein[, featnameallshort[c(min1,min2,min3,min4,min5,min6,min7,min8,min9)]], yte)
           } 
        
if (debug) print("9")
if (0) {
           tmp.resub <- sigclust(trin[, featnameallshort[c(min1,min2,min3,min4,min5,min6,min7,min8,min9)]], 100, labflag=1, label=as.vector(est.resub>thr), icovest=2)
           sig9.resub[1:2,f] <- c(tmp.resub@xcindex, -log10(tmp.resub@pvalnorm))
           tmp.cv <- sigclust(tein[, featnameallshort[c(min1,min2,min3,min4,min5,min6,min7,min8,min9)]], 100, labflag=1, label=as.vector(est.cv>thr), icovest=2)
           sig9.cv[1:2,f] <- c(tmp.cv@xcindex, -log10(tmp.cv@pvalnorm))
}
           stepbystep.featname[1:9,f] = featnameallshort[c( min1,min2,min3,min4,min5,min6,min7,min8,min9)]
           
           
        }
        if (ncol(zin)>=10 && nf>=10) {   
           minerr10=1;
           x1 = trin[, featnameallshort[min1]]
           x2 = trin[, featnameallshort[min2]]
           x3 = trin[, featnameallshort[min3]]
           x4 = trin[, featnameallshort[min4]]
           x5 = trin[, featnameallshort[min5]]
           x6 = trin[, featnameallshort[min6]]
           x7 = trin[, featnameallshort[min7]]
           x8 = trin[, featnameallshort[min8]]
           x9 = trin[, featnameallshort[min9]]
           for (i in 1:nfeat) {
                if (!(i %in% c(min1,min2,min3,min4,min5,min6,min7,min8,min9))) {
                    xi = trin[, featnameallshort[i]]
                    fiti = lm(ytr ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + xi)
                    ierr = mean(abs(fiti$residuals))
                    if (ierr < minerr10) {
                        minerr10 = ierr 
                        min10 = i 
                    minfiti = fiti
                    }
                }
           }
        
           if (1) {
           err10.resub[f] = minerr10 ;
           est.resub = minfiti$coefficients %*%  t(cbind(iintr, trin[, featnameallshort[c(min1,min2,min3,min4,min5,min6,min7,min8,min9,min10)]]))
           est.cv    = minfiti$coefficients %*%  t(cbind(iinte, tein[, featnameallshort[c(min1,min2,min3,min4,min5,min6,min7,min8,min9,min10)]]))
           err10.cv[f]    = mean(abs(est.cv - gtte))
           ct10.cv[1:2,f]    = get_sen_spc(est.cv, gtte, thr)
           ct10.resub[1:2,f] = get_sen_spc(est.resub, gttr, thr)
           } else {
           err10.resub[f] = get_resub_err(fiti, ytr) ;
           err10.cv[f]    = get_cv_err(fiti, tein[, featnameallshort[c(min1,min2,min3,min4,min5,min6,min7,min8,min9,min10)]], yte)
           }
        
if (debug) print("10")
if (0) {
           tmp.resub <- sigclust(trin[, featnameallshort[c(min1,min2,min3,min4,min5,min6,min7,min8,min9,min10)]], 100, labflag=1, label=as.vector(est.resub>thr), icovest=2)
           sig10.resub[1:2,f] <- c(tmp.resub@xcindex, -log10(tmp.resub@pvalnorm))
           tmp.cv <- sigclust(tein[, featnameallshort[c(min1,min2,min3,min4,min5,min6,min7,min8,min9,min10)]], 100, labflag=1, label=as.vector(est.cv>thr), icovest=2)
           sig10.cv[1:2,f] <- c(tmp.cv@xcindex, -log10(tmp.cv@pvalnorm))
}

           stepbystep.featname[1:10,f] = featnameallshort[c( min1,min2,min3,min4,min5,min6,min7,min8,min9,min10)]

        } 
      } # if(length(idxte)>0 && length(idxtr)>0) {
        
    }  # nfold
        
        
    avg.sig1.cv <- apply(sig1.cv, 1, mean)
    avg.sig2.cv <- apply(sig2.cv, 1, mean) 
    avg.sig3.cv <- apply(sig3.cv, 1, mean) 
    avg.sig4.cv <- apply(sig4.cv, 1, mean) 
    avg.sig5.cv <- apply(sig5.cv, 1, mean) 
    avg.sig6.cv <- apply(sig6.cv, 1, mean) 
    avg.sig7.cv <- apply(sig7.cv, 1, mean) 
    avg.sig8.cv <- apply(sig8.cv, 1, mean) 
    avg.sig9.cv <- apply(sig9.cv, 1, mean) 
    avg.sig10.cv <- apply(sig10.cv, 1, mean) 
    
    sd.sig1.cv <- apply(sig1.cv, 1, sd) 
    sd.sig2.cv <- apply(sig2.cv, 1, sd) 
    sd.sig3.cv <- apply(sig3.cv, 1, sd) 
    sd.sig4.cv <- apply(sig4.cv, 1, sd) 
    sd.sig5.cv <- apply(sig5.cv, 1, sd) 
    sd.sig6.cv <- apply(sig6.cv, 1, sd) 
    sd.sig7.cv <- apply(sig7.cv, 1, sd) 
    sd.sig8.cv <- apply(sig8.cv, 1, sd) 
    sd.sig9.cv <- apply(sig9.cv, 1, sd) 
    sd.sig10.cv <- apply(sig10.cv, 1, sd) 
   
    avg.sig1.resub <- apply(sig1.resub, 1, mean) 
    avg.sig2.resub <- apply(sig2.resub, 1, mean) 
    avg.sig3.resub <- apply(sig3.resub, 1, mean) 
    avg.sig4.resub <- apply(sig4.resub, 1, mean) 
    avg.sig5.resub <- apply(sig5.resub, 1, mean) 
    avg.sig6.resub <- apply(sig6.resub, 1, mean) 
    avg.sig7.resub <- apply(sig7.resub, 1, mean) 
    avg.sig8.resub <- apply(sig8.resub, 1, mean) 
    avg.sig9.resub <- apply(sig9.resub, 1, mean) 
    avg.sig10.resub <- apply(sig10.resub, 1, mean) 
        
    sd.sig1.resub <- apply(sig1.resub, 1, sd) 
    sd.sig2.resub <- apply(sig2.resub, 1, sd) 
    sd.sig3.resub <- apply(sig3.resub, 1, sd) 
    sd.sig4.resub <- apply(sig4.resub, 1, sd) 
    sd.sig5.resub <- apply(sig5.resub, 1, sd) 
    sd.sig6.resub <- apply(sig6.resub, 1, sd) 
    sd.sig7.resub <- apply(sig7.resub, 1, sd) 
    sd.sig8.resub <- apply(sig8.resub, 1, sd) 
    sd.sig9.resub <- apply(sig9.resub, 1, sd) 
    sd.sig10.resub <- apply(sig10.resub, 1, sd) 
        
    avg.ct1.cv <- apply(ct1.cv, 1, mean)
    avg.ct2.cv <- apply(ct2.cv, 1, mean) 
    avg.ct3.cv <- apply(ct3.cv, 1, mean) 
    avg.ct4.cv <- apply(ct4.cv, 1, mean) 
    avg.ct5.cv <- apply(ct5.cv, 1, mean) 
    avg.ct6.cv <- apply(ct6.cv, 1, mean) 
    avg.ct7.cv <- apply(ct7.cv, 1, mean) 
    avg.ct8.cv <- apply(ct8.cv, 1, mean) 
    avg.ct9.cv <- apply(ct9.cv, 1, mean) 
    avg.ct10.cv <- apply(ct10.cv, 1, mean) 
    
    sd.ct1.cv <- apply(ct1.cv, 1, sd) 
    sd.ct2.cv <- apply(ct2.cv, 1, sd) 
    sd.ct3.cv <- apply(ct3.cv, 1, sd) 
    sd.ct4.cv <- apply(ct4.cv, 1, sd) 
    sd.ct5.cv <- apply(ct5.cv, 1, sd) 
    sd.ct6.cv <- apply(ct6.cv, 1, sd) 
    sd.ct7.cv <- apply(ct7.cv, 1, sd) 
    sd.ct8.cv <- apply(ct8.cv, 1, sd) 
    sd.ct9.cv <- apply(ct9.cv, 1, sd) 
    sd.ct10.cv <- apply(ct10.cv, 1, sd) 
   
    avg.ct1.resub <- apply(ct1.resub, 1, mean) 
    avg.ct2.resub <- apply(ct2.resub, 1, mean) 
    avg.ct3.resub <- apply(ct3.resub, 1, mean) 
    avg.ct4.resub <- apply(ct4.resub, 1, mean) 
    avg.ct5.resub <- apply(ct5.resub, 1, mean) 
    avg.ct6.resub <- apply(ct6.resub, 1, mean) 
    avg.ct7.resub <- apply(ct7.resub, 1, mean) 
    avg.ct8.resub <- apply(ct8.resub, 1, mean) 
    avg.ct9.resub <- apply(ct9.resub, 1, mean) 
    avg.ct10.resub <- apply(ct10.resub, 1, mean) 
        
    sd.ct1.resub <- apply(ct1.resub, 1, sd) 
    sd.ct2.resub <- apply(ct2.resub, 1, sd) 
    sd.ct3.resub <- apply(ct3.resub, 1, sd) 
    sd.ct4.resub <- apply(ct4.resub, 1, sd) 
    sd.ct5.resub <- apply(ct5.resub, 1, sd) 
    sd.ct6.resub <- apply(ct6.resub, 1, sd) 
    sd.ct7.resub <- apply(ct7.resub, 1, sd) 
    sd.ct8.resub <- apply(ct8.resub, 1, sd) 
    sd.ct9.resub <- apply(ct9.resub, 1, sd) 
    sd.ct10.resub <- apply(ct10.resub, 1, sd) 
        
        
    stepbystep.feat.mean.resub <- cbind( avg.ct1.resub , avg.ct2.resub , avg.ct3.resub , avg.ct4.resub , avg.ct5.resub , avg.ct6.resub , avg.ct7.resub , avg.ct8.resub , avg.ct9.resub , avg.ct10.resub)
    stepbystep.feat.mean.cv <- cbind(avg.ct1.cv , avg.ct2.cv, avg.ct3.cv, avg.ct4.cv, avg.ct5.cv, avg.ct6.cv, avg.ct7.cv, avg.ct8.cv, avg.ct9.cv, avg.ct10.cv) 
        
    stepbystep.feat.sd.resub <- cbind(sd.ct1.resub , sd.ct2.resub , sd.ct3.resub , sd.ct4.resub , sd.ct5.resub , sd.ct6.resub , sd.ct7.resub , sd.ct8.resub , sd.ct9.resub , sd.ct10.resub)
    stepbystep.feat.sd.cv <- cbind(sd.ct1.cv , sd.ct2.cv, sd.ct3.cv, sd.ct4.cv, sd.ct5.cv, sd.ct6.cv, sd.ct7.cv, sd.ct8.cv, sd.ct9.cv, sd.ct10.cv)
        
    spcsen.cv <- matrix(0.0,nrow=2,ncol=nf)
    spcsen.resub <- matrix(0.0,nrow=2,ncol=nf)
    sd.spcsen.cv <- matrix(0.0,nrow=2,ncol=nf)
    sd.spcsen.resub <- matrix(0.0,nrow=2,ncol=nf)
 if (nf>=1) {    
    spcsen.cv[,1]    <- c(mean(1-ct1.cv[2,1:nfold]), mean(ct1.cv[1,1:nfold]))
    spcsen.resub[,1]    <- c(mean(1-ct1.resub[2,1:nfold]), mean(ct1.resub[1,1:nfold]))
    sd.spcsen.cv[,1]    <- c(sd(1-ct1.cv[2,1:nfold]), sd(ct1.cv[1,1:nfold]))
    sd.spcsen.resub[,1]    <- c(sd(1-ct1.resub[2,1:nfold]), sd(ct1.resub[1,1:nfold]))
    } 
 if (nf>=2) {    
    spcsen.cv[,2]    <- c(mean(1-ct2.cv[2,1:nfold]), mean(ct2.cv[1,1:nfold]))
    spcsen.resub[,2]    <- c(mean(1-ct2.resub[2,1:nfold]), mean(ct2.resub[1,1:nfold]))
    sd.spcsen.cv[,2]    <- c(sd(1-ct2.cv[2,1:nfold]), sd(ct2.cv[1,1:nfold]))
    sd.spcsen.resub[,2]    <- c(sd(1-ct2.resub[2,1:nfold]), sd(ct2.resub[1,1:nfold]))
    } 
 if (nf>=3) {    
    spcsen.cv[,3]    <- c(mean(1-ct3.cv[2,1:nfold]), mean(ct3.cv[1,1:nfold]))
    spcsen.resub[,3]    <- c(mean(1-ct3.resub[2,1:nfold]), mean(ct3.resub[1,1:nfold]))
    sd.spcsen.cv[,3]    <- c(sd(1-ct3.cv[2,1:nfold]), sd(ct3.cv[1,1:nfold]))
    sd.spcsen.resub[,3]    <- c(sd(1-ct3.resub[2,1:nfold]), sd(ct3.resub[1,1:nfold]))
    } 
 if (nf>=4) {    
    spcsen.cv[,4]    <- c(mean(1-ct4.cv[2,1:nfold]), mean(ct4.cv[1,1:nfold]))
    spcsen.resub[,4]    <- c(mean(1-ct4.resub[2,1:nfold]), mean(ct4.resub[1,1:nfold]))
    sd.spcsen.cv[,4]    <- c(sd(1-ct4.cv[2,1:nfold]), sd(ct4.cv[1,1:nfold]))
    sd.spcsen.resub[,4]    <- c(sd(1-ct4.resub[2,1:nfold]), sd(ct4.resub[1,1:nfold]))
    } 
 if (nf>=5) {    
    spcsen.cv[,5]    <- c(mean(1-ct5.cv[2,1:nfold]), mean(ct5.cv[1,1:nfold]))
    spcsen.resub[,5]    <- c(mean(1-ct5.resub[2,1:nfold]), mean(ct5.resub[1,1:nfold]))
    sd.spcsen.cv[,5]    <- c(sd(1-ct5.cv[2,1:nfold]), sd(ct5.cv[1,1:nfold]))
    sd.spcsen.resub[,5]    <- c(sd(1-ct5.resub[2,1:nfold]), sd(ct5.resub[1,1:nfold]))
    } 
 if (nf>=6) {    
    spcsen.cv[,6]    <- c(mean(1-ct6.cv[2,1:nfold]), mean(ct6.cv[1,1:nfold]))
    spcsen.resub[,6]    <- c(mean(1-ct6.resub[2,1:nfold]), mean(ct6.resub[1,1:nfold]))
    sd.spcsen.cv[,6]    <- c(sd(1-ct6.cv[2,1:nfold]), sd(ct6.cv[1,1:nfold]))
    sd.spcsen.resub[,6]    <- c(sd(1-ct6.resub[2,1:nfold]), sd(ct6.resub[1,1:nfold]))
    } 
 if (nf>=7) {    
    spcsen.cv[,7]    <- c(mean(1-ct7.cv[2,1:nfold]), mean(ct7.cv[1,1:nfold]))
    spcsen.resub[,7]    <- c(mean(1-ct7.resub[2,1:nfold]), mean(ct7.resub[1,1:nfold]))
    sd.spcsen.cv[,7]    <- c(sd(1-ct7.cv[2,1:nfold]), sd(ct7.cv[1,1:nfold]))
    sd.spcsen.resub[,7]    <- c(sd(1-ct7.resub[2,1:nfold]), sd(ct7.resub[1,1:nfold]))
    } 
 if (nf>=8) {    
    spcsen.cv[,8]    <- c(mean(1-ct8.cv[2,1:nfold]), mean(ct8.cv[1,1:nfold]))
    spcsen.resub[,8]    <- c(mean(1-ct8.resub[2,1:nfold]), mean(ct8.resub[1,1:nfold]))
    sd.spcsen.cv[,8]    <- c(sd(1-ct8.cv[2,1:nfold]), sd(ct8.cv[1,1:nfold]))
    sd.spcsen.resub[,8]    <- c(sd(1-ct8.resub[2,1:nfold]), sd(ct8.resub[1,1:nfold]))
    } 
 if (nf>=9) {    
    spcsen.cv[,9]    <- c(mean(1-ct9.cv[2,1:nfold]), mean(ct9.cv[1,1:nfold]))
    spcsen.resub[,9]    <- c(mean(1-ct9.resub[2,1:nfold]), mean(ct9.resub[1,1:nfold]))
    sd.spcsen.cv[,9]    <- c(sd(1-ct9.cv[2,1:nfold]), sd(ct9.cv[1,1:nfold]))
    sd.spcsen.resub[,9]    <- c(sd(1-ct9.resub[2,1:nfold]), sd(ct9.resub[1,1:nfold]))
    } 
 if (nf>=10) {    
    spcsen.cv[,10]    <- c(mean(1-ct10.cv[2,1:nfold]), mean(ct10.cv[1,1:nfold]))
    spcsen.resub[,10]    <- c(mean(1-ct10.resub[2,1:nfold]), mean(ct10.resub[1,1:nfold]))
    sd.spcsen.cv[,10]    <- c(sd(1-ct10.cv[2,1:nfold]), sd(ct10.cv[1,1:nfold]))
    sd.spcsen.resub[,10]    <- c(sd(1-ct10.resub[2,1:nfold]), sd(ct10.resub[1,1:nfold]))
    } 

if (debug) {
    print("spcsen.cv")
    print(spcsen.cv)
    print("sd.spcsen.cv")
    print(sd.spcsen.cv)
    print("spcsen.resub")
    print(spcsen.resub)
    print("sd.spcsen.resub")
    print(sd.spcsen.resub)
}
    pvalcindex.cv <- matrix(0.0,nrow=2,ncol=nf)
    pvalcindex.resub <- matrix(0.0,nrow=2,ncol=nf)
    sd.pvalcindex.cv <- matrix(0.0,nrow=2,ncol=nf)
    sd.pvalcindex.resub <- matrix(0.0,nrow=2,ncol=nf)
    
 if (nf>=1) {    
    pvalcindex.cv[,1]    <- c(mean(sig1.cv[2,1:nfold]), mean(sig1.cv[1,1:nfold]))
    pvalcindex.resub[,1]    <- c(mean(sig1.resub[2,1:nfold]), mean(sig1.resub[1,1:nfold]))
    sd.pvalcindex.cv[,1]    <- c(sd(sig1.cv[2,1:nfold]), sd(sig1.cv[1,1:nfold]))
    sd.pvalcindex.resub[,1]    <- c(sd(sig1.resub[2,1:nfold]), sd(sig1.resub[1,1:nfold]))
    } 
 if (nf>=2) {    
    pvalcindex.cv[,2]    <- c(mean(sig2.cv[2,1:nfold]), mean(sig2.cv[1,1:nfold]))
    pvalcindex.resub[,2]    <- c(mean(sig2.resub[2,1:nfold]), mean(sig2.resub[1,1:nfold]))
    sd.pvalcindex.cv[,2]    <- c(sd(sig2.cv[2,1:nfold]), sd(sig2.cv[1,1:nfold]))
    sd.pvalcindex.resub[,2]    <- c(sd(sig2.resub[2,1:nfold]), sd(sig2.resub[1,1:nfold]))

    } 
 if (nf>=3) {    
    pvalcindex.cv[,3]    <- c(mean(sig3.cv[2,1:nfold]), mean(sig3.cv[1,1:nfold]))
    pvalcindex.resub[,3]    <- c(mean(sig3.resub[2,1:nfold]), mean(sig3.resub[1,1:nfold]))
    sd.pvalcindex.cv[,3]    <- c(sd(sig3.cv[2,1:nfold]), sd(sig3.cv[1,1:nfold]))
    sd.pvalcindex.resub[,3]    <- c(sd(sig3.resub[2,1:nfold]), sd(sig3.resub[1,1:nfold]))

    } 
 if (nf>=4) {    
    pvalcindex.cv[,4]    <- c(mean(sig4.cv[2,1:nfold]), mean(sig4.cv[1,1:nfold]))
    pvalcindex.resub[,4]    <- c(mean(sig4.resub[2,1:nfold]), mean(sig4.resub[1,1:nfold]))
    sd.pvalcindex.cv[,4]    <- c(sd(sig4.cv[2,1:nfold]), sd(sig4.cv[1,1:nfold]))
    sd.pvalcindex.resub[,4]    <- c(sd(sig4.resub[2,1:nfold]), sd(sig4.resub[1,1:nfold]))

    } 
 if (nf>=5) {    
    pvalcindex.cv[,5]    <- c(mean(sig5.cv[2,1:nfold]), mean(sig5.cv[1,1:nfold]))
    pvalcindex.resub[,5]    <- c(mean(sig5.resub[2,1:nfold]), mean(sig5.resub[1,1:nfold]))
    sd.pvalcindex.cv[,5]    <- c(sd(sig5.cv[2,1:nfold]), sd(sig5.cv[1,1:nfold]))
    sd.pvalcindex.resub[,5]    <- c(sd(sig5.resub[2,1:nfold]), sd(sig5.resub[1,1:nfold]))

    } 
 if (nf>=6) {    
    pvalcindex.cv[,6]    <- c(mean(sig6.cv[2,1:nfold]), mean(sig6.cv[1,1:nfold]))
    pvalcindex.resub[,6]    <- c(mean(sig6.resub[2,1:nfold]), mean(sig6.resub[1,1:nfold]))
    sd.pvalcindex.cv[,6]    <- c(sd(sig6.cv[2,1:nfold]), sd(sig6.cv[1,1:nfold]))
    sd.pvalcindex.resub[,6]    <- c(sd(sig6.resub[2,1:nfold]), sd(sig6.resub[1,1:nfold]))

    } 
 if (nf>=7) {    
    pvalcindex.cv[,7]    <- c(mean(sig7.cv[2,1:nfold]), mean(sig7.cv[1,1:nfold]))
    pvalcindex.resub[,7]    <- c(mean(sig7.resub[2,1:nfold]), mean(sig7.resub[1,1:nfold]))
    sd.pvalcindex.cv[,7]    <- c(sd(sig7.cv[2,1:nfold]), sd(sig7.cv[1,1:nfold]))
    sd.pvalcindex.resub[,7]    <- c(sd(sig7.resub[2,1:nfold]), sd(sig7.resub[1,1:nfold]))

    } 
 if (nf>=8) {    
    pvalcindex.cv[,8]    <- c(mean(sig8.cv[2,1:nfold]), mean(sig8.cv[1,1:nfold]))
    pvalcindex.resub[,8]    <- c(mean(sig8.resub[2,1:nfold]), mean(sig8.resub[1,1:nfold]))
    sd.pvalcindex.cv[,8]    <- c(sd(sig8.cv[2,1:nfold]), sd(sig8.cv[1,1:nfold]))
    sd.pvalcindex.resub[,8]    <- c(sd(sig8.resub[2,1:nfold]), sd(sig8.resub[1,1:nfold]))

    } 
 if (nf>=9) {    
    pvalcindex.cv[,9]    <- c(mean(sig9.cv[2,1:nfold]), mean(sig9.cv[1,1:nfold]))
    pvalcindex.resub[,9]    <- c(mean(sig9.resub[2,1:nfold]), mean(sig9.resub[1,1:nfold]))
    sd.pvalcindex.cv[,9]    <- c(sd(sig9.cv[2,1:nfold]), sd(sig9.cv[1,1:nfold]))
    sd.pvalcindex.resub[,9]    <- c(sd(sig9.resub[2,1:nfold]), sd(sig9.resub[1,1:nfold]))

    } 
 if (nf>=10) {    
    pvalcindex.cv[,10]    <- c(mean(sig10.cv[2,1:nfold]), mean(sig10.cv[1,1:nfold]))
    pvalcindex.resub[,10]    <- c(mean(sig10.resub[2,1:nfold]), mean(sig10.resub[1,1:nfold]))
    sd.pvalcindex.cv[,10]    <- c(sd(sig10.cv[2,1:nfold]), sd(sig10.cv[1,1:nfold]))
    sd.pvalcindex.resub[,10]    <- c(sd(sig10.resub[2,1:nfold]), sd(sig10.resub[1,1:nfold]))
 }   
if (debug) {
    print("stepbystep.featname")
    print(stepbystep.featname)
}
    myout <- list()
    myout$stepbystep.featname <- stepbystep.featname[1:nf,]
    myout$stepbystep.feat.mean.resub <- stepbystep.feat.mean.resub
    myout$stepbystep.feat.sd.resub <- stepbystep.feat.sd.resub
    myout$stepbystep.feat.mean.cv <- stepbystep.feat.mean.cv
    myout$stepbystep.feat.sd.cv <- stepbystep.feat.sd.cv
    myout$spcsen.cv <- spcsen.cv
    myout$spcsen.resub <- spcsen.resub
    myout$sd.spcsen.cv <- sd.spcsen.cv
    myout$sd.spcsen.resub <- sd.spcsen.resub
    myout$pvalcindex.cv <- pvalcindex.cv
    myout$pvalcindex.resub <- pvalcindex.resub
    myout$sd.pvalcindex.cv <- sd.pvalcindex.cv
    myout$sd.pvalcindex.resub <- sd.pvalcindex.resub
    myout$nf <- nf
    
if (flagPlot) {print("	making plots ") }
    pdf(rocpdffn, width=8, height=8)
    par(mfrow=c(1,2))
    plot(1:nf, spcsen.cv[1,1:nf], pch=16, col="red", main="false positive rate \n (1 - specificity)", ylab="false positive rate", xlab="number of features", ylim=c(0,0.5), type="o", cex=0.8)
    arrows(1:nf, y0=spcsen.cv[1,1:nf]-sd.spcsen.cv[1,1:nf], y1=spcsen.cv[1,1:nf]+sd.spcsen.cv[1,1:nf], code=3, angle=90, length=.05, col="red")
    lines(1:nf, spcsen.resub[1,1:nf], pch=16, col="blue", type="o")
    arrows(1:nf, y0=spcsen.resub[1,1:nf]-sd.spcsen.resub[1,1:nf], y1=spcsen.resub[1,1:nf]+sd.spcsen.resub[1,1:nf], code=3, angle=90, length=.05, col="blue")
    grid()
    legend("topleft", legend=c(paste(nfold, " fold cv"), "resub"), pch=c(16,16), col=c("red", "blue"))
    plot(1:nf, spcsen.cv[2,1:nf], pch=16, col="red", main="sensitivity", ylab="sensitivity", xlab="number of features", ylim=c(0.5, 1), type="o", cex=0.8)
    arrows(1:nf, y0=spcsen.cv[2,1:nf]-sd.spcsen.cv[2,1:nf], y1=spcsen.cv[2,1:nf]+sd.spcsen.cv[2,1:nf], code=3, angle=90, length=.05, col="red")
    lines(1:nf, spcsen.resub[2,1:nf], pch=16, col="blue", type="o")
    arrows(1:nf, y0=spcsen.resub[2,1:nf]-sd.spcsen.resub[2,1:nf], y1=spcsen.resub[2,1:nf]+sd.spcsen.resub[2,1:nf], code=3, angle=90, length=.05, col="blue")
    grid()
    legend("bottomright", legend=c(paste(nfold, "fold cv"), "resub"), pch=c(16,16), col=c("red", "blue"))

    plot(1:nf, pvalcindex.cv[1,1:nf], pch=16, col="red", main="single cluster : pvalue", ylab="-log10(pval)", xlab="number of features", ylim=c(0,max(c(pvalcindex.cv[1,1:nf]))), type="o", cex=0.8)
    arrows(1:nf, y0=pvalcindex.cv[1,1:nf]-sd.pvalcindex.cv[1,1:nf], y1=pvalcindex.cv[1,1:nf]+sd.pvalcindex.cv[1,1:nf], code=3, angle=90, length=.05, col="red")
    lines(1:nf, pvalcindex.resub[1,1:nf], pch=16, col="blue", type="o")
    arrows(1:nf, y0=pvalcindex.resub[1,1:nf]-sd.pvalcindex.resub[1,1:nf], y1=pvalcindex.resub[1,1:nf]+sd.pvalcindex.resub[1,1:nf], code=3, angle=90, length=.05, col="blue")
    grid()
    legend("topleft", legend=c(paste(nfold, " fold cv"), "resub"), pch=c(16,16), col=c("red", "blue"))
    plot(1:nf, pvalcindex.cv[2,1:nf], pch=16, col="red", main="cluster_index", ylab="cindex", xlab="number of features", ylim=c(0, 1), type="o", cex=0.8)
    arrows(1:nf, y0=pvalcindex.cv[2,1:nf]-sd.pvalcindex.cv[2,1:nf], y1=pvalcindex.cv[2,1:nf]+sd.pvalcindex.cv[2,1:nf], code=3, angle=90, length=.05, col="red")
    lines(1:nf, pvalcindex.resub[2,1:nf], pch=16, col="blue", type="o")
    arrows(1:nf, y0=pvalcindex.resub[2,1:nf]-sd.pvalcindex.resub[2,1:nf], y1=pvalcindex.resub[2,1:nf]+sd.pvalcindex.resub[2,1:nf], code=3, angle=90, length=.05, col="blue")
    grid()
    legend("bottomright", legend=c(paste(nfold, "fold cv"), "resub"), pch=c(16,16), col=c("red", "blue"))
    dev.off()
#
if (0) {
    print("in mystepbystep") 
    print(sss.str)
    print("dim(out[[sss.str]]$spcsen.cv")
    print(dim(out[[sss.str]]$spcsen.cv))
    print("dim(out[[sss.str]]$spcsen.resub")
    print(dim(out[[sss.str]]$spcsen.resub))
}


if (debug) {print("	returning from mystepbystepcv") }
myout


}




rmse <- function (err) {
    sqrt(mean(err^2))
}

cal_anovap <- function(ZZZ, grp) {
    ##############################################################
    print("#####   p-value for each Creline    #####")
    ##############################################################
    feat.sorted <- sort(colnames(ZZZ))
    feat.sorted.binarylast <- c(feat.sorted[-c(9,10,11)], feat.sorted[c(9,10,11)])
    colororder <- order(creC[1:10])
    anovap <- rep(0, ncol(ZZZ))
    for (i in 1:ncol(ZZZ)) {
         iZZZ <- ZZZ[, feat.sorted.binarylast[i]]
         igrp <- grp
         iDF <- data.frame(iZZZ, igrp)
         anovap[i] <- summary(aov(iZZZ ~ igrp, data=iDF))[[1]][1,5]
    }
    tmp <- mt.rawp2adjp(anovap, proc="Bonferroni")
    #tmp <- mt.rawp2adjp(anovap, proc="TSBH")
    anova.adjp <- tmp$adjp[order(tmp$index),2]
    anova.adjp
}

