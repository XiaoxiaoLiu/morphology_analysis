library(WGCNA)
library(flashClust)

#####################################

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
    print(c(length(x), length(tmp)))
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
  #print(bins)
  rf.pred.cl = setNames(rep(NA, length(cl)), names(cl))
  rf.pred.prob = matrix(0, nrow=length(cl), ncol=length(unique(cl)))
  rownames(rf.pred.prob) <- names(cl)
  colnames(rf.pred.prob) <- unique(cl)
  elecmp2<-c()
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


  filename.Id <- AssignedID_FN
  out <- read.csv(filename.Id, header=TRUE)
  Id.name <- gsub("~sim", "",as.character(out[,1]))
  Id.cre <- get.field(as.character(out[,1]), ";",1)
  Id.id <- as.character(out[,2])
  names(Id.id) <- Id.name

  XXX <- Feature
  #rownames(XXX) <- gsub("-197628.06.01.01;NA", ";197628.06.01.01", gsub("-197628.06.02.01;NA", ";197628.06.02.01", gsub("-197628.03.02.01;NA", ";197628.03.02.01", rownames(XXX))))
  nfeat <- ncol(XXX)
  idx <- match(rownames(XXX), names(Id.id))

  XXX.ClusterID <- Id.id[idx]
  unique.ClusterID <- sort(unique(XXX.ClusterID))
  XXX.label <- match(XXX.ClusterID, unique.ClusterID)
  names(XXX.label) <- names(Id.id)[idx]
  unique.label <- unique(XXX.label)

  set.seed(0) ###random number seed
  #excludeclusters<-c(0,11,13,2,4,7)  ###clusters to exclude in the cross-validation, first round
  excludeclusters<-c() ###clusters to exclude in the cross-validation, second round
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

  creline_table = read.csv(paste0(data_dir,"Cre_line_typeC.09102015.csv"), header=TRUE)
  creC <- as.character(creline_table[,"cre_line"])
  mycre <- get.field(rownames(fulltab.reordered), ";",1)
  crecolor <- rainbow(11+1)[match(mycre,creC)]
  leg.pch <- rep(15,11)
  leg.str <- creC
  leg.crecolor <- rainbow(12)[1:11]

  myheat.color <- rev(heat.colors(101)) ; myheat.color[101] <- "black"
  if (flag.plot) {
    pdf(paste(outdir, outfile3, sep="/"), height=8, width=16)
    heatmap.3(t(fulltab.reordered), ColSideColors=crecolor, Colv=TRUE, Rowv=FALSE, trace="none", keysize=0.8, margins=c(10,10), col=myheat.color, margin=c(10,10), dendrogram="none", main=" 10-fold cross validation (100 runs) ")
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

  #colnames(fulltab.reordered.woNA.R) <- paste(class.N2R[,3],colnames(fulltab.reordered.woNA.R),sep="%")
  colnames(fulltab.reordered.woNA.R) <- paste(colnames(fulltab.reordered.woNA.R),' (',class.N2R[,3],'%)')
  if (flag.plot) {
    heatmap.3(t(fulltab.reordered.woNA.R), ColSideColors=crecolor, Colv=TRUE, Rowv=FALSE, trace="none", keysize=0.8, margins=c(10,10), col=myheat.color, margin=c(10,10), dendrogram="none", main="Cluster Validation ; 100 run of 10-fold cross validation")
    legend("bottomleft", col=leg.crecolor, pch=rep(15,11), legend=leg.str,cex=0.75)
    dev.off()
  }
}















