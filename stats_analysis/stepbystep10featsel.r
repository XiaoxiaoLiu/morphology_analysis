source("/home/xiaoxiaol/work/src/cell-type-analysis/stats_analysis/analysis.function.Ephys.r")
options(warn=-1)


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


BuildBinaryClusTree <- function(Nodei, soi.leg, Nshuffled=1000, flagDEC="LDA", flagGRP="SEL", flagPthr=0.05, flagSparse=FALSE, flagPlot=TRUE, flagDIP=FALSE, flagMembership=TRUE, flagPartition="PCA" ) {
  TN.Sample.thr <- 5  # termination condition
  print("##########################################################")
  print(paste("Numbef of samples :", nrow(Nodei$Xsoi),  "&   Number of Features : ", ncol(Nodei$Xsoi)))
  if (!is.matrix(Nodei$Xsoi) || nrow(Nodei$Xsoi) <= TN.Sample.thr) {
    Nodei$terminal <- TRUE
    Assign_NodeID (Nodei)
    print(paste("  ", Nodei$NodeStr, ": terminal node...STOP # Not enough Sample", nrow(Nodei$Xoi))) 
    
  }
 else 
{
     #  NodeStr = empty?
    print(paste("   BuildBinaryClusTree", Nodei$NodeStr))
    Nodei <- clustering2decR(Nodei, Nodei$Xsoi, Nodei$GT_cre, Nodei$ZSel35, Nodei$strin, Nodei$pch, soi.leg$pch, Nodei$crecolor, Nodei$hybrid.crecolorL, 
                             soi.leg$crecolor, soi.leg$crecolorL, soi.leg$str, soi.leg$strL, Nodei$NodeStr, Nshuffled, flagDEC, flagGRP, flagPthr, flagSparse, flagPlot, flagDIP, flagMembership, flagPartition) 
    
    print(paste("   BinaryClusTree is Built with ", Nodei$NodeStr))
    
    
    #        Nodei$Enough <- est_partition_boot(Nodei, seq(10,50,5), 100, 4, Node0.strin) 
    
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
      Nodei$Left <- BuildBinaryClusTree ( Nodei$Left, soi.leg, Nshuffled, flagDEC, flagGRP, flagPthr, flagSparse, flagPlot, flagDIP, flagMembership, flagPartition )
      
      print("Right") 
      NodeRStr <- paste(Nodei$Str, "R", sep="")
      Nodei$Right <- SetUpChild( Nodei, Nodei$idxR, "R" )
      Nodei$Right <- BuildBinaryClusTree ( Nodei$Right, soi.leg, Nshuffled, flagDEC, flagGRP, flagPthr, flagSparse, flagPlot, flagDIP, flagMembership )
      
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

order_feat_byTtest <- function(Zin, Gin, Xin) {
  
  #    print(length(Gin))
  #    print("apply_Ttest") ; print(table(Gin))
  #print(Zin)
  TPval <- apply(Zin, 2, apply_Ttest, Gin)
  #    print("apply_Ttest done") ; print(table(Gin))
  tmp <- mt.rawp2adjp(TPval,'TSBH')
  TadjP <- tmp$adjp[order(tmp$index),2]
  porder <- order(TadjP, decreasing=FALSE)
  
  #    print("apply_TtestM12") ; print(table(Gin))
  Xmean12 <- apply(Xin, 2, apply_TtestM12, Gin)
  #    print("apply_TtestM12 done") ; 
  
  OUT <- list()
  OUT$order <- porder
  OUT$pval.ordered <- TadjP[porder] 
  OUT$pval.lt0.01 <- sum(TadjP[porder]<0.01, na.rm=TRUE) 
  OUT$pval.lt0.05 <- sum(TadjP[porder]<0.05, na.rm=TRUE) 
  OUT$Xmean12 <- Xmean12[,porder] 
  
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
    meth="ward"
  }
  newdistmat<-matrix(0,nrow=nrow(Pin$x),ncol=nrow(Pin$x));
  coordmat <- matrix(0,nrow=nrow(Pin$x),ncol=nPC)
  
  print("====== cell-cell distance matrix ======")
  #    pdf("testnPC1.pdf")
  for (ii in 1:nPC) {
    tempdist<-outer(Pin$x[,ii],Pin$x[,ii],"-");
    if (meth=='ward') {
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
    hhh_PCA <- hclust(as.dist(newdistmat),method="ward")
  }
  print("====== binary partition L R and its significance (binary_partition_PC)======")
  ####### cutree
  hhh_grpLR <- cutree(hhh_PCA,k=2) ;
  hhh_ttt <- table(hhh_grpLR)
  
  OUT <- list()
  if (min(hhh_ttt) <= 1) {
    OUT$cindex <- 1 
    OUT$pval   <- 1
    print(paste("Cindex=", OUT$cindex, "    pval=", OUT$pval))
    
  } else {
    print("calling sigclust")
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


binary_partition_incPC <- function (Zin, Pin, nPC, meth, flagSparse=FALSE) {
  
  debug <- 0
  if (debug==1) {
    Pin <- orig.Pin
    PDFFN <- "Node0.PCA.pval.cindex.pdf" 
    nPC <- 5
    meth="ward"
  }
  if (ncol(Pin$x) < 10) { 
    nPC <- ncol(Pin$x) 
  } else {
    nPC <- 10
  }
  
  newdistmat<-matrix(0,nrow=nrow(Pin$x),ncol=nrow(Pin$x));
  coordmat <- matrix(0,nrow=nrow(Pin$x),ncol=nPC)
  
  #    pdf("testnPC1.pdf")
  print("====== with increasing number of PC, pick most partitioining nPC  ======")
  min.pval <- 1
  min.nPC <- 1
  min.grp12 <- c()
  min.hhh_PCA <- c()
  min.cindex <- 0.0
  for (ii in 1:nPC) {
    tempdist<-outer(Pin$x[,ii],Pin$x[,ii],"-");
    if (meth=='ward') {
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
        ###### use ward by default
        hhh_PCA <- hclust(as.dist(newdistmat),method="ward")
      }
      hhh_grpLR <- cutree(hhh_PCA,k=2) ;
      hhh_ttt <- table(hhh_grpLR)
      
      if (min(hhh_ttt) > 1) {
        #statistical significance of clustering
        print("sigclust")
        sigval <- sigclust(coordmat, 100, labflag=0, icovest=2)  
        
        
        #FEM
#        print("kyle_fem")
       # femval <- kyle_fem(as.matrix(Zin), hhh_grpLR)  # 
        
        
       # print(paste("sigclust", sigval@pvalnorm, sigval@xcindex, ", femval", femval@icl, femval@k1icl, femval@pvalnorm, femval@k))
        if (sigval@pvalnorm < min.pval) {
          min.pval <- sigval@pvalnorm
          min.cindex <- sigval@xcindex
          min.grp12 <- hhh_grpLR 
          min.hhh_PCA <- hhh_PCA 
          min.nPC <- ii
        }
      }
      else
      {
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

part2grp_PCA <- function(Pin, Zin, PDFFN, Nsim=1000, flagSparse, flagPlot) {
  
  debug <- 0
  if (debug==1) {
    Pin <- orig.Pin
    Zin <- orig.ZSel35
    PDFFN <- "Node0.PCA.pval.cindex.pdf" 
  }
  meth <- "ward"
  sigval_min <- 1.0
  xcindex_min <- 1.0
  
  print(paste("      get the number of significant PCs by shuffling", Nsim, " ======") )
  NN1 <- ncol(Pin$x)
  #NN1 <- sum(get_signif_PC(Zin, Pin, Nsim))
  print(paste("      Nfeature", NN1))
  
  if (NN1 > 0) {
    #partLR <- binary_partition_PC (Pin, NN1, "ward", flagSparse) 
    partLR <- binary_partition_incPC (Zin, Pin, NN1, "ward", flagSparse) 
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


binary_partition_DLM <- function (Zin, Pin, nPC, meth, flagSparse=FALSE) {
  
  debug <- 0
  if (debug==1) {
    Pin <- orig.Pin
    PDFFN <- "Node0.PCA.pval.cindex.pdf" 
    nPC <- 5
    meth="ward"
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
    if (meth=='ward') {
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
        hhh_PCA <- hclust(as.dist(newdistmat),method="ward")
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
      print(paste("PC", ii, " eigenvalue=", summary(Pin)$importance[2,ii], "min.pval=",min.pval))
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

part2grp_DLM <- function(Pin, Zin, PDFFN, Nsim=1000, flagSparse, flagPlot) {
  
  debug <- 0
  if (debug==1) {
    Pin <- orig.Pin
    Zin <- orig.ZSel35
    PDFFN <- "Node0.PCA.pval.cindex.pdf" 
  }
  meth <- "ward"
  sigval_min <- 1.0
  xcindex_min <- 1.0
  
  print(paste("      get the number of significant PCs by shuffling", Nsim, " ======") )
  NN1 <- ncol(Pin$x)
  #NN1 <- sum(get_signif_PC(Zin, Pin, Nsim))
  print(paste("      Nfeature", NN1))
  
  if (NN1 > 0) {
    #partLR <- binary_partition_PC (Pin, NN1, "ward", flagSparse) 
    partLR <- binary_partition_DLM (Zin, Pin, NN1, "ward", flagSparse) 
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

get_grp_mean_sd <- function(PCn, grp12, Npc) {
  
  idx1 <- which(grp12==1)
  idx2 <- which(grp12==2)
  
  X1 <- as.matrix(PCn[idx1,1:Npc])
  X2 <- as.matrix(PCn[idx2,1:Npc])
  
  mu1 <- apply(X1, 2, mean)
  sd1 <- apply(X1, 2, sd)
  
  mu2 <- apply(X2, 2, mean)
  sd2 <- apply(X2, 2, sd)
  
  OUT <- list()
  OUT$mu1 <- mu1
  OUT$sd1 <- sd1
  OUT$mu2 <- mu2
  OUT$sd2 <- sd2
  
  OUT
}


cal_12membership <- function (X, param, Nin) {
  
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
  
  c(mem1, mem2)
}


clustering2decR <- function (ORIG, orig.XSel,orig.GT_cre, orig.ZSel35, orig.strin, mypch, leg.pch, crecolor, crecolor2, leg.crecolor, leg.crecolor2, leg.str, leg.str2, nodestr, nsim, flagDEC="LDA", flagGRP="SEL", flagPthr=0.05, flagSparse=FALSE, flagPlot=TRUE, flagDIP=FALSE, flagMembership=FALSE, flagPartition="PCA") { 
  
  set.seed(1)
  OUT <- ORIG
  
  OUT$nPC  <- 0
  #print(dim(orig.ZSel35))
  
  pdffn <- paste(orig.strin, ".pca.pdf", sep="") 
  if (!flagPlot) save(ORIG, file="ORIG.Rdata") 
  if (flagPlot) pdf(pdffn)
 
    ################################ bimodel features are selected 
    print("      selecting features with bimodal distribution")
    nfeat <- ncol(orig.ZSel35) 
    print(paste("           number of cells :", nrow(orig.ZSel35), ",     number of features :", ncol(orig.ZSel35)))
    
    if (!flagDIP) {
      feature_std <- apply(orig.ZSel35, 2, sd)  # col std
      zero.sd <- sum(feature_std==0)
      
      #print(paste("zero.sd =", zero.sd))
      
      if (zero.sd==0) { f0 <- 2} else { f0 <- 1 }
      foi <- c() 
      if (zero.sd!=nfeat) {
        feature_std_order <- order(feature_std, decreasing=FALSE)
        for (f in c((zero.sd+f0):nfeat)) {
          Zf <- orig.ZSel35[,c(feature_std_order[1],feature_std_order[f])]
          ############ sig clust 
          fsig <- sigclust(Zf, 100, labflag=0, label=0, icovest=2)
          ############ sig clust 
          if (fsig@pvalnorm < 0.01) {
            foi <- c(foi, feature_std_order[f]) 
            print(paste('feat=', f, '  sigclust=', fsig@pvalnorm))
          }
        }
      }
    } else {
      foi <- c() 
      for (f in c(1:nfeat)) {
        frun <- dip.test(orig.ZSel35[,f], simulate.p.value=TRUE, B=1000)
        if (frun$p.value < 0.01) foi <- c(foi, f)
      } 
    }
    print(paste("           number of interesting features = ", length(foi), "out of", nfeat))
  

  ######################  PCA
  print("     Do pca")
  
  if (length(foi)>0) {
    
    ZSel35 <- orig.ZSel35[, foi]
    XSel <- orig.XSel[, foi]
    
    #orig.pca <- plotPCA(ZSel35, nodestr, mypch, leg.pch, crecolor, leg.crecolor, leg.str, pdffn, flagPlot)
    orig.pca <- prcomp(ZSel35, cor=TRUE)
    
    if (flagPartition=="PCA") {
      print("     Partition in pca")
      pdffn <- paste(orig.strin, ".pca.pdf", sep="") 
      ###########
      orig.grp <- part2grp_PCA(orig.pca, ZSel35, pdffn, nsim, flagSparse, flagPlot)
     #########    
} else {
      pdffn <- paste(orig.strin, ".DLM.pdf", sep="") 
      orig.grp <- part2grp_DLM(orig.pca, ZSel35, pdffn, nsim, flagSparse, flagPlot)
    }
    if (flagPlot) dev.off()
    
    pdffn <- paste(orig.strin, ".pca.clustering.pdf", sep="")
    
    print(paste("           (orig.grp$nPCsig, orig.grp$pval, orig.grp$pval)=",orig.grp$nPCsig, round(orig.grp$pval,3), round(orig.grp$pval,3))) 
    if ((orig.grp$nPCsig > 0) && !is.na(orig.grp$cindex) && !is.na(orig.grp$pval) && (orig.grp$pval< flagPthr)) {
      clusterPCA (orig.pca, orig.grp$nPCsig, crecolor2, leg.crecolor2, leg.str2, pdffn, flagPlot) 
      #if ((orig.grp$nPCsig > 0) && !is.na(orig.grp$pval) && (orig.grp$pval<1) && ((orig.grp$pval<0.01) || (orig.grp$gap > 0))) {
      print("            any features with discrimination power given binary partition?")
      #print(table(orig.grp$PCA_grp12))
      idx.ordered <- order_feat_byTtest (ZSel35, orig.grp$PCA_grp12, XSel)
      
      ZSel35.pordered <- ZSel35[, idx.ordered$order]
      tmp.pordered <- ZSel35.pordered 
      colnames(tmp.pordered) <- paste(signif(-log10(idx.ordered$pval.ordered),2), colnames(ZSel35)[idx.ordered$order])
      if (flagPthr==0.01) nadjp <- idx.ordered$pval.lt0.01
      if (flagPthr==0.05) nadjp <- idx.ordered$pval.lt0.05 
      print(paste("           ", nadjp, "features with adjp<", flagPthr)) 
      OUT$ttest.feat <- idx.ordered 
      
      if (nadjp > 0) {
        
        write.csv(data.frame(feature=colnames(ZSel35)[idx.ordered$order[1:nadjp]], log10P.neg=signif(-log10(idx.ordered$pval.ordered[1:nadjp]),4), meanL=idx.ordered$Xmean12[1,1:nadjp] , meanR=idx.ordered$Xmean12[2,1:nadjp]), file= paste(orig.strin, '.PCApartition.AllFeatures.adjP', flagPthr, '.sorted.csv', sep=""), row.names=FALSE)
        #write.csv(data.frame(feature=colnames(orig.ZSel35)[idx.ordered$order[1:nadjp]], log10p.neg=signif(-log10(idx.ordered$pval.ordered[1:nadjp]),4)), file= paste(orig.strin, '.pcapartition.allfeatures.adjp0.01.sorted.csv', sep=""), row.names=FALSE)
        
        if (flagPlot) {
          pdf(paste(orig.strin, '.heatmap.pcapartition.pdf', sep=""), height=5+round(ncol(ZSel35)/5), width=8+round(nrow(ZSel35)/10))
          
          hh <- heatmap.3(t(tmp.pordered), hclustfun=hhclust, Colv=as.dendrogram(orig.grp$hhh_PCA), Rowv=TRUE, trace="none",  dendrogram="column",
                          ColSideColors=crecolor2, main=paste("clustering of cells @", nodestr, " with ", orig.grp$nPCsig," pc's \n", 
                                                              " cindex =", signif(orig.grp$cindex,3), " & -log10(p) =", signif(-log10(orig.grp$pval),3), "\n",
                                                              ncol(ZSel35), " features", sep=""), 
                          keysize=0.8, margins=c(10,10), cexcol=0.5,cexrow=1, col=rbg)
          legend("bottomleft", bg="white", legend=leg.str2, pch=15, col=leg.crecolor2, cex=0.75, bty="n")
        }
        
        hhh_allFeat <- hclust(dist(tmp.pordered),method="ward")
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
          
          hhh_SelFeat <- hclust(dist(tmp.adjp.pordered),method="ward")
          hhh_1 <- hclust(dist(t(tmp.adjp.pordered)),method="ward")
          SelFeat_grpLR <- cutree(hhh_SelFeat,k=2)
          OUT$hhh_SEL <- hhh_SelFeat 
          OUT$hhh_SEL_rowInd <- hhh_SelFeat$order 
          OUT$hhh_SEL_colInd <- hhh_1$order 
        } 
        if (flagPlot) dev.off()
        
        
        if (flagMembership) { 
          print("          estimate parameter for two clusters")
          params <- get_grp_mean_sd (orig.pca$x, orig.grp$PCA_grp12, orig.grp$nPCsig)
          print("          calculate memebership score")
          #membershiplr <- cal_12membership(orig.pca$x[1,1:orig.grp$nPCsig], params, orig.grp$nPCsig)
          membershiplrs <- apply(as.matrix(orig.pca$x[,1:orig.grp$nPCsig]), 1, cal_12membership, params, orig.grp$nPCsig)
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
        order_pc1_loading <- order(abs(orig.pca$rotation[,"PC1"]), decreasing=TRUE)
        orig_feat_pcloading <- rownames(orig.pca$rotation)[order_pc1_loading]
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
        sf10.params <- get_grp_mean_sd (sf10.ZSel35, orig.grp$PCA_grp12, orig.pcloading.cv4roc$nf)
        
        OUT$nPC         <- orig.grp$nPCsig 
        OUT$loadingsPC  <- orig.pca$rotation[,1:orig.grp$nPCsig] 
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
            tmp.membershiplrs <- apply(as.matrix(sf10.ZSel35[,1:nsf]), 1, cal_12membership, sf10.params, nsf)
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
          OUT$membershiplrs_sf10 <- sf10.membershiplrs 
          OUT$grouping_pval_cindex_sf10 <- sf10.significance 
        } else {
          OUT$params_sf10 <- c()
          OUT$featname_sf10 <- c()
          OUT$membershiplrs_sf10 <- c() 
          OUT$grouping_pval_cindex_sf10 <- c() 
        } 
        OUT$params_pc   <- params
        OUT$membershiplrs_pc <- membershiplrs 
        
        
        print("calculate creline composition & draw pie chart")
        
        OUT$ct <- table(orig.GT_cre, OUT$grouping)
        print(OUT$ct)
        write.csv(OUT$ct, file=paste(orig.strin, ".cre_line.composition.csv", sep=""))
        
        if (flagPlot) {
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







mystepbystep <- function(gt, y, zin, featnameallshort, nfeat, nfold, thr) {
  nf=10
  if (0) {
    gt <- orig.class
    y <- orig.class
    zin <- orig.zsel35
    featnameallshort <- grp0_feat_pcloading
    nfeat <-38
    nfold <- 10
    thr <- i.thr 
    rocpdffn <- 'cktest.pdf'
  }
  
  y2.idx <- which(y==2)
  y2.n <- length(y2.idx)
  
  y1.idx <- which(y==1)
  y1.n <- length(y1.idx)
  
  subtitlestr <- paste("nfold=", nfold) 
  
  err1.cv <- err2.cv <- err3.cv <- err4.cv <- err5.cv <- err6.cv <- err7.cv <- err8.cv <- err9.cv <- err10.cv <- rep(0,nfold+1) 
  err1.resub <- err2.resub <- err3.resub <- err4.resub <- err5.resub <- err6.resub <- err7.resub <- err8.resub <- err9.resub <- err10.resub <- rep(0,nfold+1) 
  ct1.cv <- ct2.cv <- ct3.cv <- ct4.cv <- ct5.cv <- ct6.cv <- ct7.cv <- ct8.cv <- ct9.cv <- ct10.cv <- matrix(0,2,nfold+1) 
  ct1.resub <- ct2.resub <- ct3.resub <- ct4.resub <- ct5.resub <- ct6.resub <- ct7.resub <- ct8.resub <- ct9.resub <- ct10.resub <- matrix(0,2,nfold+1) 
  
  stepbystep.featname = matrix("hello", 10, nfold+1)
  nsample = nrow(zin)
  
  mymod1 = 1:y1.n %% nfold
  mymod2 = 1:y2.n %% nfold
  for (f in 1:(nfold+1)) {
    #print(paste('thr=', thr, ' ', f, '-th fold /', nfold, 'folds'))
    if (f<=nfold) {
      idxte = c(y1.idx[which(mymod1 == (f-1))], y2.idx[which(mymod2 == (f-1))])
      idxtr = setdiff(1:nsample,idxte)
    } else {
      idxte = 1:nsample
      idxtr = 1:nsample
    }
    trin = zin[idxtr,]
    tein = zin[idxte,]
    iintr  = matrix(1,length(idxtr),1)
    iinte  = matrix(1,length(idxte),1)
    ytr = y[idxtr]
    yte = y[idxte]
    gttr = gt[idxtr]
    gtte = gt[idxte]
    
    if (ncol(zin)>=1) {
      minerr1=1
      print(paste("nfeat=",nfeat))
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
        print("min1 and featnameallshort(min1)")
        print(paste(min1, featnameallshort[c(min1)]))
        err1.resub[f] = minerr1 ;
        est.resub = minfiti$coefficients %*%  t(cbind(iintr, trin[, featnameallshort[c(min1)]]))
        est.cv    = minfiti$coefficients %*%  t(cbind(iinte, tein[, featnameallshort[c(min1)]]))
        err1.cv[f]    = mean(abs(est.cv - yte))
        ct1.cv[1:2,f]    = get_sen_spc(est.cv, gtte, thr)
        ct1.resub[1:2,f] = get_sen_spc(est.resub, gttr, thr)
      } else {
        err1.resub[f] = get_resub_err(fiti, ytr) ;
        err1.cv[f]    = get_cv_err(fiti, tein[, featnameallshort[c(min1)]], yte)
      }
      rm(est.cv, est.resub, minfiti, fiti, xi)
    }
    
    if (ncol(zin)>=2) {   
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
        err2.cv[f]    = mean(abs(est.cv - yte))
        ct2.cv[1:2,f]    = get_sen_spc(est.cv, gtte, thr)
        ct2.resub[1:2,f] = get_sen_spc(est.resub, gttr, thr)
      } else {
        err2.resub[f] = get_resub_err(fiti, ytr) ;
        err2.cv[f]    = get_cv_err(fiti, tein[, featnameallshort[c(min1,min2)]], yte)
      } 
      
      
      rm(est.cv, est.resub, minfiti, fiti, xi)
    }
    if (ncol(zin)>=3) {   
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
        err3.cv[f]    = mean(abs(est.cv - yte))
        ct3.cv[1:2,f]    = get_sen_spc(est.cv, gtte, thr)
        ct3.resub[1:2,f] = get_sen_spc(est.resub, gttr, thr)
      } else {
        err3.resub[f] = get_resub_err(fiti, ytr) ;
        err3.cv[f]    = get_cv_err(fiti, tein[, featnameallshort[c(min1,min2,min3)]], yte)
      } 
      
      rm(est.cv, est.resub, minfiti, fiti, xi)
    }
    if (ncol(zin)>=4) {   
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
        err4.cv[f]    = mean(abs(est.cv - yte))
        ct4.cv[1:2,f]    = get_sen_spc(est.cv, gtte, thr)
        ct4.resub[1:2,f] = get_sen_spc(est.resub, gttr, thr)
      } else {
        err4.resub[f] = get_resub_err(fiti, ytr) ;
        err4.cv[f]    = get_cv_err(fiti, tein[, featnameallshort[c(min1,min2,min3,min4)]], yte)
      }
      
      rm(est.cv, est.resub, minfiti, fiti, xi)
    }
    if (ncol(zin)>=5) {   
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
        err5.cv[f]    = mean(abs(est.cv - yte))
        ct5.cv[1:2,f]    = get_sen_spc(est.cv, gtte, thr)
        ct5.resub[1:2,f] = get_sen_spc(est.resub, gttr, thr)
      } else {
        err5.resub[f] = get_resub_err(fiti, ytr) ;
        err5.cv[f]    = get_cv_err(fiti, tein[, featnameallshort[c(min1,min2,min3,min4,min5)]], yte)
      } 
      
      rm(est.cv, est.resub, minfiti, fiti, xi)
    }
    if (ncol(zin)>=6) {   
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
        err6.cv[f]    = mean(abs(est.cv - yte))
        ct6.cv[1:2,f]    = get_sen_spc(est.cv, gtte, thr)
        ct6.resub[1:2,f] = get_sen_spc(est.resub, gttr, thr)
      } else {
        err6.resub[f] = get_resub_err(fiti, ytr) ;
        err6.cv[f]    = get_cv_err(fiti, tein[, featnameallshort[c(min1,min2,min3,min4,min5,min6)]], yte)
      }
      
      rm(est.cv, est.resub, minfiti, fiti, xi)
    }
    if (ncol(zin)>=7) {   
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
        err7.cv[f]    = mean(abs(est.cv - yte))
        ct7.cv[1:2,f]    = get_sen_spc(est.cv, gtte, thr)
        ct7.resub[1:2,f] = get_sen_spc(est.resub, gttr, thr)
      } else {
        err7.resub[f] = get_resub_err(fiti, ytr) ;
        err7.cv[f]    = get_cv_err(fiti, tein[, featnameallshort[c(min1,min2,min3,min4,min5,min6,min7)]], yte)
      }
      
      rm(est.cv, est.resub, minfiti, fiti, xi)
      
    }
    if (ncol(zin)>=8) {   
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
        err8.cv[f]    = mean(abs(est.cv - yte))
        ct8.cv[1:2,f]    = get_sen_spc(est.cv, gtte, thr)
        ct8.resub[1:2,f] = get_sen_spc(est.resub, gttr, thr)
      } else {
        err8.resub[f] = get_resub_err(fiti, ytr) ;
        err8.cv[f]    = get_cv_err(fiti, tein[, featnameallshort[c(min1,min2,min3,min4,min5,min6,min7,min8)]], yte)
      } 
      
      rm(est.cv, est.resub, minfiti, fiti, xi)
    }
    if (ncol(zin)>=9) {   
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
        err9.cv[f]    = mean(abs(est.cv - yte))
        ct9.cv[1:2,f]    = get_sen_spc(est.cv, gtte, thr)
        ct9.resub[1:2,f] = get_sen_spc(est.resub, gttr, thr)
      } else {
        err9.resub[f] = get_resub_err(fiti, ytr) ;
        err9.cv[f]    = get_cv_err(fiti, tein[, featnameallshort[c(min1,min2,min3,min4,min5,min6,min7,min8,min9)]], yte)
      } 
      
      rm(est.cv, est.resub, minfiti, fiti, xi)
      
    }
    if (ncol(zin)>=10) {   
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
        err10.cv[f]    = mean(abs(est.cv - yte))
        ct10.cv[1:2,f]    = get_sen_spc(est.cv, gtte, thr)
        ct10.resub[1:2,f] = get_sen_spc(est.resub, gttr, thr)
      } else {
        err10.resub[f] = get_resub_err(fiti, ytr) ;
        err10.cv[f]    = get_cv_err(fiti, tein[, featnameallshort[c(min1,min2,min3,min4,min5,min6,min7,min8,min9,min10)]], yte)
      }
      rm(est.cv, est.resub, minfiti, fiti, xi)
      
      stepbystep.featname[1:10,f] = featnameallshort[c( min1,min2,min3,min4,min5,min6,min7,min8,min9,min10)]
    } #10
    
  } # fold
  
  
  spcsen.cv <- matrix(0.0,nrow=2,ncol=nf)
  spcsen.resub <- matrix(0.0,nrow=2,ncol=nf)
  sd.spcsen.cv <- matrix(0.0,nrow=2,ncol=nf)
  sd.spcsen.resub <- matrix(0.0,nrow=2,ncol=nf)
  
  spcsen.cv[,1]    <- c(mean(1-ct1.cv[1,1:nfold]), mean(ct1.cv[2,1:nfold]))
  spcsen.cv[,2]    <- c(mean(1-ct2.cv[1,1:nfold]), mean(ct2.cv[2,1:nfold]))
  spcsen.cv[,3]    <- c(mean(1-ct3.cv[1,1:nfold]), mean(ct3.cv[2,1:nfold]))
  spcsen.cv[,4]    <- c(mean(1-ct4.cv[1,1:nfold]), mean(ct4.cv[2,1:nfold]))
  spcsen.cv[,5]    <- c(mean(1-ct5.cv[1,1:nfold]), mean(ct5.cv[2,1:nfold]))
  spcsen.cv[,6]    <- c(mean(1-ct6.cv[1,1:nfold]), mean(ct6.cv[2,1:nfold]))
  spcsen.cv[,7]    <- c(mean(1-ct7.cv[1,1:nfold]), mean(ct7.cv[2,1:nfold]))
  spcsen.cv[,8]    <- c(mean(1-ct8.cv[1,1:nfold]), mean(ct8.cv[2,1:nfold]))
  spcsen.cv[,9]    <- c(mean(1-ct9.cv[1,1:nfold]), mean(ct9.cv[2,1:nfold]))
  spcsen.cv[,10]    <- c(mean(1-ct10.cv[1,1:nfold]), mean(ct10.cv[2,1:nfold]))
  
  
  spcsen.resub[,1]    <- c(mean(1-ct1.resub[1,1:nfold]), mean(ct1.resub[2,1:nfold]))
  spcsen.resub[,2]    <- c(mean(1-ct2.resub[1,1:nfold]), mean(ct2.resub[2,1:nfold]))
  spcsen.resub[,3]    <- c(mean(1-ct3.resub[1,1:nfold]), mean(ct3.resub[2,1:nfold]))
  spcsen.resub[,4]    <- c(mean(1-ct4.resub[1,1:nfold]), mean(ct4.resub[2,1:nfold]))
  spcsen.resub[,5]    <- c(mean(1-ct5.resub[1,1:nfold]), mean(ct5.resub[2,1:nfold]))
  spcsen.resub[,6]    <- c(mean(1-ct6.resub[1,1:nfold]), mean(ct6.resub[2,1:nfold]))
  spcsen.resub[,7]    <- c(mean(1-ct7.resub[1,1:nfold]), mean(ct7.resub[2,1:nfold]))
  spcsen.resub[,8]    <- c(mean(1-ct8.resub[1,1:nfold]), mean(ct8.resub[2,1:nfold]))
  spcsen.resub[,9]    <- c(mean(1-ct9.resub[1,1:nfold]), mean(ct9.resub[2,1:nfold]))
  spcsen.resub[,10]    <- c(mean(1-ct10.resub[1,1:nfold]), mean(ct10.resub[2,1:nfold]))
  
  
  sd.spcsen.cv[,1]    <- c(sd(1-ct1.cv[1,1:nfold]), sd(ct1.cv[2,1:nfold]))
  sd.spcsen.cv[,2]    <- c(sd(1-ct2.cv[1,1:nfold]), sd(ct2.cv[2,1:nfold]))
  sd.spcsen.cv[,3]    <- c(sd(1-ct3.cv[1,1:nfold]), sd(ct3.cv[2,1:nfold]))
  sd.spcsen.cv[,4]    <- c(sd(1-ct4.cv[1,1:nfold]), sd(ct4.cv[2,1:nfold]))
  sd.spcsen.cv[,5]    <- c(sd(1-ct5.cv[1,1:nfold]), sd(ct5.cv[2,1:nfold]))
  sd.spcsen.cv[,6]    <- c(sd(1-ct6.cv[1,1:nfold]), sd(ct6.cv[2,1:nfold]))
  sd.spcsen.cv[,7]    <- c(sd(1-ct7.cv[1,1:nfold]), sd(ct7.cv[2,1:nfold]))
  sd.spcsen.cv[,8]    <- c(sd(1-ct8.cv[1,1:nfold]), sd(ct8.cv[2,1:nfold]))
  sd.spcsen.cv[,9]    <- c(sd(1-ct9.cv[1,1:nfold]), sd(ct9.cv[2,1:nfold]))
  sd.spcsen.cv[,10]    <- c(sd(1-ct10.cv[1,1:nfold]), sd(ct10.cv[2,1:nfold]))
  
  
  sd.spcsen.resub[,1]    <- c(sd(1-ct1.resub[1,1:nfold]), sd(ct1.resub[2,1:nfold]))
  sd.spcsen.resub[,2]    <- c(sd(1-ct2.resub[1,1:nfold]), sd(ct2.resub[2,1:nfold]))
  sd.spcsen.resub[,3]    <- c(sd(1-ct3.resub[1,1:nfold]), sd(ct3.resub[2,1:nfold]))
  sd.spcsen.resub[,4]    <- c(sd(1-ct4.resub[1,1:nfold]), sd(ct4.resub[2,1:nfold]))
  sd.spcsen.resub[,5]    <- c(sd(1-ct5.resub[1,1:nfold]), sd(ct5.resub[2,1:nfold]))
  sd.spcsen.resub[,6]    <- c(sd(1-ct6.resub[1,1:nfold]), sd(ct6.resub[2,1:nfold]))
  sd.spcsen.resub[,7]    <- c(sd(1-ct7.resub[1,1:nfold]), sd(ct7.resub[2,1:nfold]))
  sd.spcsen.resub[,8]    <- c(sd(1-ct8.resub[1,1:nfold]), sd(ct8.resub[2,1:nfold]))
  sd.spcsen.resub[,9]    <- c(sd(1-ct9.resub[1,1:nfold]), sd(ct9.resub[2,1:nfold]))
  sd.spcsen.resub[,10]    <- c(sd(1-ct10.resub[1,1:nfold]), sd(ct10.resub[2,1:nfold]))
  
  pdf("roc.final.test.in.correct.pdf", width=8, height=8)
  par(mfrow=c(1,2))
  plot(1:nf, spcsen.cv[1,], pch=16, col="red", main="false positive rate \n (1 - specificity)", ylab="false positive rate(%)", xlab="number of features", ylim=c(0,5), type="o", cex=0.8)
  lines(1:nf, spcsen.resub[1,], pch=16, col="darkgreen", type="o")
  grid()
  legend("topleft", legend=c(paste(nfold, " fold cv"), "resub"), pch=c(16,16), col=c("red", "darkgreen"))
  
  plot(1:nf, spcsen.cv[2,], pch=16, col="red", main="sensitivity", ylab="sensitivity", xlab="number of features", ylim=c(50,100), type="o", cex=0.8)
  lines(1:nf, spcsen.resub[2,], pch=16, col="darkgreen", type="o")
  grid()
  legend("bottomright", legend=c(paste(nfold, "fold cv"), "resub"), pch=c(16,16), col=c("red", "darkgreen"))
  dev.off()
  
  myout <- list()
  myout$stepbystep.featname <- stepbystep.featname
  myout$spcsen.cv <- spcsen.cv
  myout$spcsen.resub <- spcsen.resub
  myout$sd.spcsen.cv <- sd.spcsen.cv
  myout$sd.spcsen.resub <- sd.spcsen.resub
  myout
  
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
  
  if (0) {
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
          #      print("min1 and featnameallshort(min1)")
          #      print(paste(min1, featnameallshort[c(min1)]))
          #      print("dim(tein)")
          #      print(dim(tein))
          #      print("minfiti$coefficients")
          #      print(minfiti$coefficients)
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


mystepbystepboot <- function(gt, y, zin, featnameallshort, nfeat, samplesizes, nboot, nfold, thr, rocpdfn) {
  
  
  #if (1) {
  #        y <- h0_class
  #        zin <- h0_zsel35
  #        featnameallshort <- grp0_feat_pcloading
  #        nfeat <-39
  #        nfold <- 10
  #        thr <- 1 
  #        rocpdffn <- 'cktest.pdf'
  #}
  nf=10
  
  y2.idx <- which(y==2)
  y2.n <- length(y2.idx)
  
  y1.idx <- which(y==1)
  y1.n <- length(y1.idx)
  
  
  nsamplesize = length(samplesizes)
  subtitlestr <- paste("nsamplesize=", nsamplesize) 
  
  err1.cv <- err2.cv <- err3.cv <- err4.cv <- err5.cv <- err6.cv <- err7.cv <- err8.cv <- err9.cv <- err10.cv <- rep(0,nboot+1) 
  err1.resub <- err2.resub <- err3.resub <- err4.resub <- err5.resub <- err6.resub <- err7.resub <- err8.resub <- err9.resub <- err10.resub <- rep(0,nboot+1) 
  ct1.cv <- ct2.cv <- ct3.cv <- ct4.cv <- ct5.cv <- ct6.cv <- ct7.cv <- ct8.cv <- ct9.cv <- ct10.cv <- matrix(0,ncol=nboot,nrow=2)
  ct1.resub <- ct2.resub <- ct3.resub <- ct4.resub <- ct5.resub <- ct6.resub <- ct7.resub <- ct8.resub <- ct9.resub <- ct10.resub <- matrix(0,ncol=nboot,nrow=2) 
  
  stepbystep.featname = matrix("hello", 10, nsamplesize+1)
  nsample = nrow(zin)
  
  mymod1 = 1:y1.n %% nfold
  mymod2 = 1:y2.n %% nfold
  
  out <- list()
  for (sss in 1:nsamplesize) {
    
    f.samplesize = samplesizes[sss]
    
    print(paste('samplesize=', f.samplesize))
    sss.str <- paste('samplesize', f.samplesize, sep="")
    for (f in 1:nboot) {
      idxtr = c(sample(y2.idx,f.samplesize), y1.idx) 
      idxte = c(setdiff(1:nsample,idxtr), y1.idx)
      
      trin = zin[idxtr,]
      tein = zin[idxte,]
      iintr  = matrix(1,length(idxtr),1)
      iinte  = matrix(1,length(idxte),1)
      ytr = y[idxtr]
      yte = y[idxte]
      gttr = gt[idxtr]
      gtte = gt[idxte]
      
      if (ncol(zin)>=1) {
        minerr1=1
        for (i in 1:nfeat) {
          xi = trin[, featnameallshort[i]]
          fiti = lm(gttr ~ xi)
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
        rm(est.cv, est.resub, minfiti, fiti, xi)
      }
      if (ncol(zin)>=2) {   
        minerr2=1;
        x1 = trin[, featnameallshort[min1]]
        for (i in 1:nfeat) {
          if (!(i %in% min1)) {
            xi = trin[, featnameallshort[i]]
            fiti = lm(gttr ~ x1 + xi)
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
        
        
        rm(est.cv, est.resub, minfiti, fiti, xi)
      }
      if (ncol(zin)>=3) {   
        minerr3=1;
        x1 = trin[, featnameallshort[min1]]
        x2 = trin[, featnameallshort[min2]]
        for (i in 1:nfeat) {
          if (!(i %in% c(min1,min2))) {
            xi = trin[, featnameallshort[i]]
            fiti = lm(gttr ~ x1 + x2 + xi)
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
        
        rm(est.cv, est.resub, minfiti, fiti, xi)
      }
      if (ncol(zin)>=4) {   
        minerr4=1;
        x1 = trin[, featnameallshort[min1]]
        x2 = trin[, featnameallshort[min2]]
        x3 = trin[, featnameallshort[min3]]
        for (i in 1:nfeat) {
          if (!(i %in% c(min1,min2,min3))) {
            xi = trin[, featnameallshort[i]]
            fiti = lm(gttr ~ x1 + x2 + x3 + xi)
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
        
        rm(est.cv, est.resub, minfiti, fiti, xi)
      }
      if (ncol(zin)>=5) {   
        minerr5=1;
        x1 = trin[, featnameallshort[min1]]
        x2 = trin[, featnameallshort[min2]]
        x3 = trin[, featnameallshort[min3]]
        x4 = trin[, featnameallshort[min4]]
        for (i in 1:nfeat) {
          if (!(i %in% c(min1,min2,min3,min4))) {
            xi = trin[, featnameallshort[i]]
            fiti = lm(gttr ~ x1 + x2 + x3 + x4 + xi)
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
        
        rm(est.cv, est.resub, minfiti, fiti, xi)
      }
      if (ncol(zin)>=6) {   
        minerr6=1;
        x1 = trin[, featnameallshort[min1]]
        x2 = trin[, featnameallshort[min2]]
        x3 = trin[, featnameallshort[min3]]
        x4 = trin[, featnameallshort[min4]]
        x5 = trin[, featnameallshort[min5]]
        for (i in 1:nfeat) {
          if (!(i %in% c(min1,min2,min3,min4,min5))) {
            xi = trin[, featnameallshort[i]]
            fiti = lm(gttr ~ x1 + x2 + x3 + x4 + x5 + xi)
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
        
        rm(est.cv, est.resub, minfiti, fiti, xi)
      }
      if (ncol(zin)>=7) {   
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
            fiti = lm(gttr ~ x1 + x2 + x3 + x4 + x5 + x6 + xi)
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
        
        rm(est.cv, est.resub, minfiti, fiti, xi)
        
      }
      if (ncol(zin)>=8) {   
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
            fiti = lm(gttr ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + xi)
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
        
        rm(est.cv, est.resub, minfiti, fiti, xi)
      }
      if (ncol(zin)>=9) {   
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
            fiti = lm(gttr ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + xi)
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
        
        rm(est.cv, est.resub, minfiti, fiti, xi)
        
      }
      if (ncol(zin)>=10) {   
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
            fiti = lm(gttr ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + xi)
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
        rm(est.cv, est.resub, minfiti, fiti, xi)
        
        stepbystep.featname[1:10,sss] = featnameallshort[c( min1,min2,min3,min4,min5,min6,min7,min8,min9,min10)]
      } 
      
    }  # nboot
    
    
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
    
    spcsen.cv[,1]    <- c(mean(1-ct1.cv[2,1:nboot]), mean(ct1.cv[1,1:nboot]))
    spcsen.cv[,2]    <- c(mean(1-ct2.cv[2,1:nboot]), mean(ct2.cv[1,1:nboot]))
    spcsen.cv[,3]    <- c(mean(1-ct3.cv[2,1:nboot]), mean(ct3.cv[1,1:nboot]))
    spcsen.cv[,4]    <- c(mean(1-ct4.cv[2,1:nboot]), mean(ct4.cv[1,1:nboot]))
    spcsen.cv[,5]    <- c(mean(1-ct5.cv[2,1:nboot]), mean(ct5.cv[1,1:nboot]))
    spcsen.cv[,6]    <- c(mean(1-ct6.cv[2,1:nboot]), mean(ct6.cv[1,1:nboot]))
    spcsen.cv[,7]    <- c(mean(1-ct7.cv[2,1:nboot]), mean(ct7.cv[1,1:nboot]))
    spcsen.cv[,8]    <- c(mean(1-ct8.cv[2,1:nboot]), mean(ct8.cv[1,1:nboot]))
    spcsen.cv[,9]    <- c(mean(1-ct9.cv[2,1:nboot]), mean(ct9.cv[1,1:nboot]))
    spcsen.cv[,10]    <- c(mean(1-ct10.cv[2,1:nboot]), mean(ct10.cv[1,1:nboot]))
    
    
    spcsen.resub[,1]    <- c(mean(1-ct1.resub[2,1:nboot]), mean(ct1.resub[1,1:nboot]))
    spcsen.resub[,2]    <- c(mean(1-ct2.resub[2,1:nboot]), mean(ct2.resub[1,1:nboot]))
    spcsen.resub[,3]    <- c(mean(1-ct3.resub[2,1:nboot]), mean(ct3.resub[1,1:nboot]))
    spcsen.resub[,4]    <- c(mean(1-ct4.resub[2,1:nboot]), mean(ct4.resub[1,1:nboot]))
    spcsen.resub[,5]    <- c(mean(1-ct5.resub[2,1:nboot]), mean(ct5.resub[1,1:nboot]))
    spcsen.resub[,6]    <- c(mean(1-ct6.resub[2,1:nboot]), mean(ct6.resub[1,1:nboot]))
    spcsen.resub[,7]    <- c(mean(1-ct7.resub[2,1:nboot]), mean(ct7.resub[1,1:nboot]))
    spcsen.resub[,8]    <- c(mean(1-ct8.resub[2,1:nboot]), mean(ct8.resub[1,1:nboot]))
    spcsen.resub[,9]    <- c(mean(1-ct9.resub[2,1:nboot]), mean(ct9.resub[1,1:nboot]))
    spcsen.resub[,10]    <- c(mean(1-ct10.resub[2,1:nboot]), mean(ct10.resub[1,1:nboot]))
    
    
    sd.spcsen.cv[,1]    <- c(sd(1-ct1.cv[2,1:nboot]), sd(ct1.cv[1,1:nboot]))
    sd.spcsen.cv[,2]    <- c(sd(1-ct2.cv[2,1:nboot]), sd(ct2.cv[1,1:nboot]))
    sd.spcsen.cv[,3]    <- c(sd(1-ct3.cv[2,1:nboot]), sd(ct3.cv[1,1:nboot]))
    sd.spcsen.cv[,4]    <- c(sd(1-ct4.cv[2,1:nboot]), sd(ct4.cv[1,1:nboot]))
    sd.spcsen.cv[,5]    <- c(sd(1-ct5.cv[2,1:nboot]), sd(ct5.cv[1,1:nboot]))
    sd.spcsen.cv[,6]    <- c(sd(1-ct6.cv[2,1:nboot]), sd(ct6.cv[1,1:nboot]))
    sd.spcsen.cv[,7]    <- c(sd(1-ct7.cv[2,1:nboot]), sd(ct7.cv[1,1:nboot]))
    sd.spcsen.cv[,8]    <- c(sd(1-ct8.cv[2,1:nboot]), sd(ct8.cv[1,1:nboot]))
    sd.spcsen.cv[,9]    <- c(sd(1-ct9.cv[2,1:nboot]), sd(ct9.cv[1,1:nboot]))
    sd.spcsen.cv[,10]    <- c(sd(1-ct10.cv[2,1:nboot]), sd(ct10.cv[1,1:nboot]))
    
    
    sd.spcsen.resub[,1]    <- c(sd(1-ct1.resub[2,1:nboot]), sd(ct1.resub[1,1:nboot]))
    sd.spcsen.resub[,2]    <- c(sd(1-ct2.resub[2,1:nboot]), sd(ct2.resub[1,1:nboot]))
    sd.spcsen.resub[,3]    <- c(sd(1-ct3.resub[2,1:nboot]), sd(ct3.resub[1,1:nboot]))
    sd.spcsen.resub[,4]    <- c(sd(1-ct4.resub[2,1:nboot]), sd(ct4.resub[1,1:nboot]))
    sd.spcsen.resub[,5]    <- c(sd(1-ct5.resub[2,1:nboot]), sd(ct5.resub[1,1:nboot]))
    sd.spcsen.resub[,6]    <- c(sd(1-ct6.resub[2,1:nboot]), sd(ct6.resub[1,1:nboot]))
    sd.spcsen.resub[,7]    <- c(sd(1-ct7.resub[2,1:nboot]), sd(ct7.resub[1,1:nboot]))
    sd.spcsen.resub[,8]    <- c(sd(1-ct8.resub[2,1:nboot]), sd(ct8.resub[1,1:nboot]))
    sd.spcsen.resub[,9]    <- c(sd(1-ct9.resub[2,1:nboot]), sd(ct9.resub[1,1:nboot]))
    sd.spcsen.resub[,10]    <- c(sd(1-ct10.resub[2,1:nboot]), sd(ct10.resub[1,1:nboot]))
    
    myout <- list()
    myout$stepbystep.featname <- stepbystep.featname
    myout$stepbystep.feat.mean.resub <- stepbystep.feat.mean.resub
    myout$stepbystep.feat.sd.resub <- stepbystep.feat.sd.resub
    myout$stepbystep.feat.mean.cv <- stepbystep.feat.mean.cv
    myout$stepbystep.feat.sd.cv <- stepbystep.feat.sd.cv
    myout$spcsen.cv <- spcsen.cv
    myout$spcsen.resub <- spcsen.resub
    myout$sd.spcsen.cv <- sd.spcsen.cv
    myout$sd.spcsen.resub <- sd.spcsen.resub
    
    par(mfrow=c(1,2))
    plot(1:nf, spcsen.cv[1,], pch=16, col="red", main="false positive rate \n (1 - specificity)", ylab="false positive rate(%)", xlab="number of features", ylim=c(0,5), type="o", cex=0.8)
    lines(1:nf, spcsen.resub[1,], pch=16, col="darkgreen", type="o")
    grid()
    legend("topleft", legend=c(paste(nfold, " fold cv"), "resub"), pch=c(16,16), col=c("red", "darkgreen"))
    plot(1:nf, spcsen.cv[2,], pch=16, col="red", main="sensitivity", ylab="sensitivity", xlab="number of features", ylim=c(50,nf), type="o", cex=0.8)
    lines(1:nf, spcsen.resub[2,], pch=16, col="darkgreen", type="o")
    grid()
    legend("bottomright", legend=c(paste(nfold, "fold cv"), "resub"), pch=c(16,16), col=c("red", "darkgreen"))
    #
    out[[sss.str]] <- myout
    if (0) {
      print("in mystepbystep") 
      print(sss.str)
      print("dim(out[[sss.str]]$spcsen.cv")
      print(dim(out[[sss.str]]$spcsen.cv))
      print("dim(out[[sss.str]]$spcsen.resub")
      print(dim(out[[sss.str]]$spcsen.resub))
    }
  } # ach samplesize
  
  out
  
}


mystepbystepboot10 <- function(gt, y, zin, featnameallshort, nfeat, samplesizes, nboot, nfold, thr, rocpdfn) {
  
  
  #if (1) {
  #        y <- h0_class
  #        zin <- h0_zsel35
  #        featnameallshort <- grp0_feat_pcloading
  #        nfeat <-39
  #        nfold <- 10
  #        thr <- 1 
  #        rocpdffn <- 'cktest.pdf'
  #}
  nf=10
  
  y2.idx <- which(y==2)
  y2.n <- length(y2.idx)
  
  y1.idx <- which(y==1)
  y1.n <- length(y1.idx)
  
  nsamplesize = length(samplesizes)
  subtitlestr <- paste("nsamplesize=", nsamplesize) 
  
  err1.cv <- err2.cv <- err3.cv <- err4.cv <- err5.cv <- err6.cv <- err7.cv <- err8.cv <- err9.cv <- err10.cv <- rep(0,nboot+1) 
  err1.resub <- err2.resub <- err3.resub <- err4.resub <- err5.resub <- err6.resub <- err7.resub <- err8.resub <- err9.resub <- err10.resub <- rep(0,nboot+1) 
  ct1.cv <- ct2.cv <- ct3.cv <- ct4.cv <- ct5.cv <- ct6.cv <- ct7.cv <- ct8.cv <- ct9.cv <- ct10.cv <- matrix(0,ncol=nboot,nrow=2)
  ct1.resub <- ct2.resub <- ct3.resub <- ct4.resub <- ct5.resub <- ct6.resub <- ct7.resub <- ct8.resub <- ct9.resub <- ct10.resub <- matrix(0,ncol=nboot,nrow=2) 
  
  stepbystep.featname = matrix("hello", 10, nsamplesize+1)
  nsample = nrow(zin)
  
  y2size <- round(y2.n/10)
  y1size <- round(y1.n/10)
  
  out <- list()
  ntest = 10
  for (fff in 1:ntest) {
    fff.idx2 =c(sample(y2.idx,y2size)) #pvalb
    fff.idx1 =c(sample(y1.idx,y1size))
    
    fff.idxte = c(fff.idx2,fff.idx1) 
    fff.idxtr2 = setdiff(y2.idx, fff.idx2)
    fff.idxtr1 = setdiff(y1.idx, fff.idx1)
    
    for (sss in 1:nsamplesize) {
      
      sss.samplesize = samplesizes[sss]
      print(paste('samplesize=', sss.samplesize))
      sss.str <- paste('test', fff, '.samplesize', sss.samplesize, sep="")
      
      for (f in 1:nboot) {
        idxtr = c(sample(fff.idxtr2,sss.samplesize), fff.idxtr1) 
        idxte = fff.idxte
        
        trin = zin[idxtr,]
        tein = zin[idxte,]
        iintr  = matrix(1,length(idxtr),1)
        iinte  = matrix(1,length(idxte),1)
        ytr = y[idxtr]
        yte = y[idxte]
        gttr = gt[idxtr]
        gtte = gt[idxte]
        
        if (ncol(zin)>=1) {
          minerr1=1
          for (i in 1:nfeat) {
            xi = trin[, featnameallshort[i]]
            fiti = lm(gttr ~ xi)
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
          rm(est.cv, est.resub, minfiti, fiti, xi)
        }
        if (ncol(zin)>=2) {   
          minerr2=1;
          x1 = trin[, featnameallshort[min1]]
          for (i in 1:nfeat) {
            if (!(i %in% min1)) {
              xi = trin[, featnameallshort[i]]
              fiti = lm(gttr ~ x1 + xi)
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
          
          
          rm(est.cv, est.resub, minfiti, fiti, xi)
        }
        if (ncol(zin)>=3) {   
          minerr3=1;
          x1 = trin[, featnameallshort[min1]]
          x2 = trin[, featnameallshort[min2]]
          for (i in 1:nfeat) {
            if (!(i %in% c(min1,min2))) {
              xi = trin[, featnameallshort[i]]
              fiti = lm(gttr ~ x1 + x2 + xi)
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
          
          rm(est.cv, est.resub, minfiti, fiti, xi)
        }
        if (ncol(zin)>=4) {   
          minerr4=1;
          x1 = trin[, featnameallshort[min1]]
          x2 = trin[, featnameallshort[min2]]
          x3 = trin[, featnameallshort[min3]]
          for (i in 1:nfeat) {
            if (!(i %in% c(min1,min2,min3))) {
              xi = trin[, featnameallshort[i]]
              fiti = lm(gttr ~ x1 + x2 + x3 + xi)
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
          
          rm(est.cv, est.resub, minfiti, fiti, xi)
        }
        if (ncol(zin)>=5) {   
          minerr5=1;
          x1 = trin[, featnameallshort[min1]]
          x2 = trin[, featnameallshort[min2]]
          x3 = trin[, featnameallshort[min3]]
          x4 = trin[, featnameallshort[min4]]
          for (i in 1:nfeat) {
            if (!(i %in% c(min1,min2,min3,min4))) {
              xi = trin[, featnameallshort[i]]
              fiti = lm(gttr ~ x1 + x2 + x3 + x4 + xi)
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
          
          rm(est.cv, est.resub, minfiti, fiti, xi)
        }
        if (ncol(zin)>=6) {   
          minerr6=1;
          x1 = trin[, featnameallshort[min1]]
          x2 = trin[, featnameallshort[min2]]
          x3 = trin[, featnameallshort[min3]]
          x4 = trin[, featnameallshort[min4]]
          x5 = trin[, featnameallshort[min5]]
          for (i in 1:nfeat) {
            if (!(i %in% c(min1,min2,min3,min4,min5))) {
              xi = trin[, featnameallshort[i]]
              fiti = lm(gttr ~ x1 + x2 + x3 + x4 + x5 + xi)
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
          
          rm(est.cv, est.resub, minfiti, fiti, xi)
        }
        if (ncol(zin)>=7) {   
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
              fiti = lm(gttr ~ x1 + x2 + x3 + x4 + x5 + x6 + xi)
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
          
          rm(est.cv, est.resub, minfiti, fiti, xi)
          
        }
        if (ncol(zin)>=8) {   
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
              fiti = lm(gttr ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + xi)
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
          
          rm(est.cv, est.resub, minfiti, fiti, xi)
        }
        if (ncol(zin)>=9) {   
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
              fiti = lm(gttr ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + xi)
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
          
          rm(est.cv, est.resub, minfiti, fiti, xi)
          
        }
        if (ncol(zin)>=10) {   
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
              fiti = lm(gttr ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + xi)
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
          rm(est.cv, est.resub, minfiti, fiti, xi)
          
          stepbystep.featname[1:10,sss] = featnameallshort[c( min1,min2,min3,min4,min5,min6,min7,min8,min9,min10)]
        } 
        
      }  # nboot
      
      
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
      
      spcsen.cv[,1]    <- c(mean(1-ct1.cv[2,1:nboot]), mean(ct1.cv[1,1:nboot]))
      spcsen.cv[,2]    <- c(mean(1-ct2.cv[2,1:nboot]), mean(ct2.cv[1,1:nboot]))
      spcsen.cv[,3]    <- c(mean(1-ct3.cv[2,1:nboot]), mean(ct3.cv[1,1:nboot]))
      spcsen.cv[,4]    <- c(mean(1-ct4.cv[2,1:nboot]), mean(ct4.cv[1,1:nboot]))
      spcsen.cv[,5]    <- c(mean(1-ct5.cv[2,1:nboot]), mean(ct5.cv[1,1:nboot]))
      spcsen.cv[,6]    <- c(mean(1-ct6.cv[2,1:nboot]), mean(ct6.cv[1,1:nboot]))
      spcsen.cv[,7]    <- c(mean(1-ct7.cv[2,1:nboot]), mean(ct7.cv[1,1:nboot]))
      spcsen.cv[,8]    <- c(mean(1-ct8.cv[2,1:nboot]), mean(ct8.cv[1,1:nboot]))
      spcsen.cv[,9]    <- c(mean(1-ct9.cv[2,1:nboot]), mean(ct9.cv[1,1:nboot]))
      spcsen.cv[,10]    <- c(mean(1-ct10.cv[2,1:nboot]), mean(ct10.cv[1,1:nboot]))
      
      
      spcsen.resub[,1]    <- c(mean(1-ct1.resub[2,1:nboot]), mean(ct1.resub[1,1:nboot]))
      spcsen.resub[,2]    <- c(mean(1-ct2.resub[2,1:nboot]), mean(ct2.resub[1,1:nboot]))
      spcsen.resub[,3]    <- c(mean(1-ct3.resub[2,1:nboot]), mean(ct3.resub[1,1:nboot]))
      spcsen.resub[,4]    <- c(mean(1-ct4.resub[2,1:nboot]), mean(ct4.resub[1,1:nboot]))
      spcsen.resub[,5]    <- c(mean(1-ct5.resub[2,1:nboot]), mean(ct5.resub[1,1:nboot]))
      spcsen.resub[,6]    <- c(mean(1-ct6.resub[2,1:nboot]), mean(ct6.resub[1,1:nboot]))
      spcsen.resub[,7]    <- c(mean(1-ct7.resub[2,1:nboot]), mean(ct7.resub[1,1:nboot]))
      spcsen.resub[,8]    <- c(mean(1-ct8.resub[2,1:nboot]), mean(ct8.resub[1,1:nboot]))
      spcsen.resub[,9]    <- c(mean(1-ct9.resub[2,1:nboot]), mean(ct9.resub[1,1:nboot]))
      spcsen.resub[,10]    <- c(mean(1-ct10.resub[2,1:nboot]), mean(ct10.resub[1,1:nboot]))
      
      
      sd.spcsen.cv[,1]    <- c(sd(1-ct1.cv[2,1:nboot]), sd(ct1.cv[1,1:nboot]))
      sd.spcsen.cv[,2]    <- c(sd(1-ct2.cv[2,1:nboot]), sd(ct2.cv[1,1:nboot]))
      sd.spcsen.cv[,3]    <- c(sd(1-ct3.cv[2,1:nboot]), sd(ct3.cv[1,1:nboot]))
      sd.spcsen.cv[,4]    <- c(sd(1-ct4.cv[2,1:nboot]), sd(ct4.cv[1,1:nboot]))
      sd.spcsen.cv[,5]    <- c(sd(1-ct5.cv[2,1:nboot]), sd(ct5.cv[1,1:nboot]))
      sd.spcsen.cv[,6]    <- c(sd(1-ct6.cv[2,1:nboot]), sd(ct6.cv[1,1:nboot]))
      sd.spcsen.cv[,7]    <- c(sd(1-ct7.cv[2,1:nboot]), sd(ct7.cv[1,1:nboot]))
      sd.spcsen.cv[,8]    <- c(sd(1-ct8.cv[2,1:nboot]), sd(ct8.cv[1,1:nboot]))
      sd.spcsen.cv[,9]    <- c(sd(1-ct9.cv[2,1:nboot]), sd(ct9.cv[1,1:nboot]))
      sd.spcsen.cv[,10]    <- c(sd(1-ct10.cv[2,1:nboot]), sd(ct10.cv[1,1:nboot]))
      
      
      sd.spcsen.resub[,1]    <- c(sd(1-ct1.resub[2,1:nboot]), sd(ct1.resub[1,1:nboot]))
      sd.spcsen.resub[,2]    <- c(sd(1-ct2.resub[2,1:nboot]), sd(ct2.resub[1,1:nboot]))
      sd.spcsen.resub[,3]    <- c(sd(1-ct3.resub[2,1:nboot]), sd(ct3.resub[1,1:nboot]))
      sd.spcsen.resub[,4]    <- c(sd(1-ct4.resub[2,1:nboot]), sd(ct4.resub[1,1:nboot]))
      sd.spcsen.resub[,5]    <- c(sd(1-ct5.resub[2,1:nboot]), sd(ct5.resub[1,1:nboot]))
      sd.spcsen.resub[,6]    <- c(sd(1-ct6.resub[2,1:nboot]), sd(ct6.resub[1,1:nboot]))
      sd.spcsen.resub[,7]    <- c(sd(1-ct7.resub[2,1:nboot]), sd(ct7.resub[1,1:nboot]))
      sd.spcsen.resub[,8]    <- c(sd(1-ct8.resub[2,1:nboot]), sd(ct8.resub[1,1:nboot]))
      sd.spcsen.resub[,9]    <- c(sd(1-ct9.resub[2,1:nboot]), sd(ct9.resub[1,1:nboot]))
      sd.spcsen.resub[,10]    <- c(sd(1-ct10.resub[2,1:nboot]), sd(ct10.resub[1,1:nboot]))
      
      myout <- list()
      myout$stepbystep.featname <- stepbystep.featname
      myout$stepbystep.feat.mean.resub <- stepbystep.feat.mean.resub
      myout$stepbystep.feat.sd.resub <- stepbystep.feat.sd.resub
      myout$stepbystep.feat.mean.cv <- stepbystep.feat.mean.cv
      myout$stepbystep.feat.sd.cv <- stepbystep.feat.sd.cv
      myout$spcsen.cv <- spcsen.cv
      myout$spcsen.resub <- spcsen.resub
      myout$sd.spcsen.cv <- sd.spcsen.cv
      myout$sd.spcsen.resub <- sd.spcsen.resub
      
      par(mfrow=c(1,2))
      plot(1:nf, spcsen.cv[1,], pch=16, col="red", main="false positive rate \n (1 - specificity)", ylab="false positive rate(%)", xlab="number of features", ylim=c(0,5), type="o", cex=0.8)
      lines(1:nf, spcsen.resub[1,], pch=16, col="darkgreen", type="o")
      grid()
      legend("topleft", legend=c(paste(nfold, " fold cv"), "resub"), pch=c(16,16), col=c("red", "darkgreen"))
      plot(1:nf, spcsen.cv[2,], pch=16, col="red", main="sensitivity", ylab="sensitivity", xlab="number of features", ylim=c(50,nf), type="o", cex=0.8)
      lines(1:nf, spcsen.resub[2,], pch=16, col="darkgreen", type="o")
      grid()
      legend("bottomright", legend=c(paste(nfold, "fold cv"), "resub"), pch=c(16,16), col=c("red", "darkgreen"))
      #
      out[[sss.str]] <- myout
      
      if (0) {
        print("in mystepbystep") 
        print(sss.str)
        print("dim(out[[sss.str]]$spcsen.cv")
        print(dim(out[[sss.str]]$spcsen.cv))
        print("dim(out[[sss.str]]$spcsen.resub")
        print(dim(out[[sss.str]]$spcsen.resub))
      }
    } # samplesize
    
  } # nfold
  
  out
  
}


est_partition <- function (h0_gt, h0_class, h0_zsel35, grp0_feat_pcloading, nallfeat, nfold, outfn) {
  nf = 10
  thrs = seq(1.4,1.6,0.05)
  
  sensitivity.cv <- fpr.cv <- matrix(0,nf,length(thrs))
  sensitivity.resub <- fpr.resub <- matrix(0,nf,length(thrs))
  sd.sensitivity.cv <- sd.fpr.cv <- matrix(0,nf,length(thrs))
  sd.sensitivity.resub <- sd.fpr.resub <- matrix(0,nf,length(thrs))
  
  for (iiii in 1:length(thrs)) {
    i.thr <- thrs[iiii]
    h0_pcloading.cv4roc <- mystepbystep(h0_gt, h0_class, h0_zsel35, grp0_feat_pcloading, nallfeat, nfold, i.thr)
    
    sensitivity.cv[1:nf,iiii] <- h0_pcloading.cv4roc$spcsen.cv[2,1:nf]
    fpr.cv[1:nf,iiii] <- h0_pcloading.cv4roc$spcsen.cv[1,1:nf]
    sensitivity.resub[1:nf,iiii] <- h0_pcloading.cv4roc$spcsen.resub[2,1:nf]
    fpr.resub[1:nf,iiii] <- h0_pcloading.cv4roc$spcsen.resub[1,1:nf]
    
    sd.sensitivity.cv[1:nf,iiii] <- h0_pcloading.cv4roc$sd.spcsen.cv[2,1:nf]
    sd.fpr.cv[1:nf,iiii] <- h0_pcloading.cv4roc$sd.spcsen.cv[1,1:nf]
    sd.sensitivity.resub[1:nf,iiii] <- h0_pcloading.cv4roc$sd.spcsen.resub[2,1:nf]
    sd.fpr.resub[1:nf,iiii] <- h0_pcloading.cv4roc$sd.spcsen.resub[1,1:nf]
  }
  
  aroc <- rep(0,nf)
  
  for (i in 1:nf) {
    aroc[i] <- cal_aroc(fpr.cv[i,], sensitivity.cv[i,])
  }
  
  pdf(outfn, width=8, height=8)
  plot(fpr.cv[1,], sensitivity.cv[1,], col="grey", pch=".", ylab="sensitivity", xlab="false positive rate = 1-specificity", main="roc of pvalb call", xlim=c(0,1),ylim=c(0,1), type="o", lty=1, lwd=2)
  for (i in 2:8) {
    lines(fpr.cv[i,], sensitivity.cv[i,], col=i, pch=".", type="o", lty=1, lwd=2)
  }
  lines(fpr.cv[9,], sensitivity.cv[9,], col=9, pch=".", type="o", lty=2,lwd=1)
  lines(fpr.cv[10,], sensitivity.cv[10,], col=10, pch=".", type="o", lty=2,lwd=1)
  thisstr <- paste(signif(aroc,3),":", 1:nf, "features")
  thiscol <- 1:nf
  legend("bottomright", legend=thisstr, col=thiscol, pch=rep(15,10) )
  dev.off()
  
  out <- list()
  out$aroc <- aroc
  out$sensitivity.cv <- sensitivity.cv
  out$frp.cv <- fpr.cv
  out$sensitivity.resub <- sensitivity.resub
  out$frp.resub <- fpr.resub
  out
}



est_partition_boot10 <- function (h0_gt, h0_class, h0_zsel35, grp0_feat_pcloading, nallfeat, samplesizes, nboot, nfold, h0_strin) {
  
  
  #print("###############################")
  mycase=""
  if(mycase=="node0") {
    h0_gt <- node0$groupingpc
    h0_class <- node0$groupingpc
    h0_zsel35 <- node0.zsel35
    grp0_feat_pcloading <- colnames(node0.zsel35)
    nallfeat <- ncol(node0.zsel35)
    samplesizes <- seq(10,50,5)
    nboot <- 100
    nfold <- 4
    h0_strin <- 'vera.node0' 
  }
  if (mycase== "noder") {
    h0_gt <- noder$groupingpc
    h0_class <- noder$groupingpc
    h0_zsel35 <- noder.zsel35
    grp0_feat_pcloading <- colnames(noder.zsel35)
    nallfeat <- ncol(noder.zsel35)
    samplesizes <- seq(10,50,5)
    nboot <- 100
    nfold <- 4
    h0_strin <- 'vera.noder'
  }
  if (mycase== "noderl") {
    h0_gt <- noderl$groupingpc
    h0_class <- noderl$groupingpc
    h0_zsel35 <- noderl.zsel35
    grp0_feat_pcloading <- colnames(noderl.zsel35)
    nallfeat <- ncol(noderl.zsel35)
    samplesizes <- seq(10,50,5)
    nboot <- 100
    nfold <- 4
    h0_strin <- 'vera.noderl'
  }
  
  outfn <- paste(h0_strin, ".run.morecells.pdf", sep="")
  pdf(outfn)
  
  nf = 10
  ns = length(samplesizes)
  thrs = seq(1.4,1.6,0.05)
  
  sensitivity.cv <- fpr.cv <- array(0, dim=c(length(thrs), nf, ns))
  sensitivity.resub <- fpr.resub <- array(0, dim=c(length(thrs), nf, ns))
  sd.sensitivity.cv <- sd.fpr.cv <- array(0, dim=c(length(thrs), nf, ns))
  sd.sensitivity.resub <- sd.fpr.resub <- array(0, dim=c(length(thrs), nf, ns))
  aroc <- sd.aroc <- matrix(0,nrow=ns,ncol=nf)
  
  
  for (iiii in 1:length(thrs)) {
    i.thr <- thrs[iiii]
    print(paste("ithr=",i.thr))
    h0_pcloading.cv4roc <- mystepbystepboot10(h0_gt, h0_class, h0_zsel35, grp0_feat_pcloading, nallfeat, samplesizes, nboot, nfold, i.thr)
    
    for (sss in 1:ns) {
      sss.samplesize = samplesizes[sss]
      sss.str <- paste('samplesize', sss.samplesize, sep="")
      
      sensitivity.cv[iiii, 1:nf,sss] <- h0_pcloading.cv4roc[[sss.str]]$spcsen.cv[2,1:nf]
      fpr.cv[iiii, 1:nf,sss] <- h0_pcloading.cv4roc[[sss.str]]$spcsen.cv[1,1:nf]
      sensitivity.resub[iiii, 1:nf,sss] <- h0_pcloading.cv4roc[[sss.str]]$spcsen.resub[2,1:nf]
      fpr.resub[iiii, 1:nf,sss] <- h0_pcloading.cv4roc[[sss.str]]$spcsen.resub[1,1:nf]
      
      sd.sensitivity.cv[iiii, 1:nf,sss] <- h0_pcloading.cv4roc[[sss.str]]$sd.spcsen.cv[2,1:nf]
      sd.fpr.cv[iiii, 1:nf,sss] <- h0_pcloading.cv4roc[[sss.str]]$sd.spcsen.cv[1,1:nf]
      sd.sensitivity.resub[iiii, 1:nf,sss] <- h0_pcloading.cv4roc[[sss.str]]$sd.spcsen.resub[2,1:nf]
      sd.fpr.resub[iiii, 1:nf,sss] <- h0_pcloading.cv4roc[[sss.str]]$sd.spcsen.resub[1,1:nf]
    }
  }
  dev.off()
  
  pdf(paste(h0.strin, ".bootstrap.summary.nb100.nf4.relgain.0.pdf", sep=""), width=8, height=8)
  for (s in 1:ns) {
    for (f in 1:nf) {
      aroc[s,f] <- cal_aroc(fpr.cv[,f,s], sensitivity.cv[,f,s])
      sd.aroc[s,f] <- mean(c(sd.fpr.cv[,f,s], sd.sensitivity.cv[,f,s]))
    }
    
    plot(fpr.cv[,1,s], sensitivity.cv[,1,s], col="grey", pch=15, ylab="sensitivity", xlab="false positive rate = 1-specificity", main=paste("roc of pvalb call :", samplesizes[s]), xlim=c(0,0.5),ylim=c(0.5,1), type="o", lty=1, lwd=3)
    for (f in 2:8) {
      lines(fpr.cv[,f,s], sensitivity.cv[,f,s], col=f, pch=15, type="o", lty=1, lwd=3)
    }
    lines(fpr.cv[,9,s], sensitivity.cv[,9,s], col=9, pch=22, type="o", lty=2,lwd=2)
    lines(fpr.cv[,10,s], sensitivity.cv[,10,s], col=10, pch=22, type="o", lty=2,lwd=2)
    thisstr <- paste(signif(aroc[s,],3),":", 1:nf, "features")
    thiscol <- 1:nf
    grid()
    legend("bottomright", legend=thisstr, col=thiscol, pch=c(rep(15,8),rep(22,2)) )
  }
  
  for (f in 1:nf) {
    if (f==1) {
      plot(samplesizes, aroc[,f], col=f, pch=15, ylab="area of roc", xlab="pvalb sample size",
           main="area of roc vs. sample size : pvalb", type="o", lwd=2, ylim=c(0.5,1.1))
      arrows(samplesizes, y0=aroc[,f]-sd.aroc[,f], y1=aroc[,f]+sd.aroc[,f], code=3, angle=90, length=.05, col=f)
      
    } else {
      if (f > 8) {
        lines(samplesizes, aroc[,f], col=f, pch=22, type="o", lwd=2,lty=2)
        arrows(samplesizes, y0=aroc[,f]-sd.aroc[,f], y1=aroc[,f]+sd.aroc[,f], code=3, angle=90, length=0.1, col=f)
      } else {
        lines(samplesizes, aroc[,f], col=f, pch=15, type="o", lwd=2,lty=1)
        arrows(samplesizes, y0=aroc[,f]-sd.aroc[,f], y1=aroc[,f]+sd.aroc[,f], code=3, angle=90, length=0.1, col=f)
      }
    }
    grid()
    thisstr <- paste(1:nf, "features")
    thiscol <- 1:nf
    legend("bottomright", legend=thisstr, col=thiscol, pch=c(rep(15,8),rep(22,2)) )
  }
  
  tmp=1
  if (tmp==1) {
    
    sen.cv.mean <-sensitivity.cv[2,,]
    sen.resub.mean <-sensitivity.resub[2,,]
    sen.cv.sd <-sd.sensitivity.cv[2,,]
    sen.resub.sd <-sd.sensitivity.resub[2,,]
    
    mean.sens.cv <-sensitivity.cv[2,,]
    mean.sens.resub <-sensitivity.resub[2,,]
    sd.sens.cv <-sd.sensitivity.cv[2,,]
    sd.sens.resub <-sd.sensitivity.resub[2,,]
    
    par(mfrow=c(2,2))
    for (f in c(1,3,7,9)) {
      plot(samplesizes, mean.sens.resub[f,], pch=15,col="black", lty=1, lwd=2,type="o", main=paste("nfeat=",f), ylim=c(0.5,1), ylab="sensitivity")
      arrows(samplesizes,y0=mean.sens.resub[f,]-(0.316*sd.sens.resub[f,]), y1=mean.sens.resub[f,]+(0.316*sd.sens.resub[f,]),
             code=3,angle=90, length=0.03, col="black")
      
      lines(samplesizes, mean.sens.cv[f,], pch=15, col="blue", lty=2, lwd=2,type="o")
      arrows(samplesizes,y0=mean.sens.cv[f,]-(0.316*sd.sens.cv[f,]), y1=mean.sens.cv[f,]+(0.316*sd.sens.cv[f,]),
             code=3,angle=90, length=0.03, col="blue")
      legend("bottomright", legend=c("training", "testing"), pch=15, col=c("black","blue"))
      grid()
    }
    
    gain.sens.resub <- t(apply(mean.sens.resub, 1, diff))
    gain.cum.sens.resub <- t(apply(gain.sens.resub, 1, cumsum))
    gain.percent.sens.resub <- signif(100*(gain.sens.resub / mean.sens.resub[,1:7]),2)
    
    gain.sens.cv <- t(apply(mean.sens.cv, 1, diff))
    gain.cum.sens.cv <- t(apply(gain.sens.cv, 1, cumsum))
    gain.percent.sens.cv <- signif(100*(gain.sens.cv / mean.sens.cv[,1:7]),2)
    
    par(mfrow=c(2,2))
    for (f in c(1,3,7,9)) {
      plot(samplesizes, c(na,gain.sens.resub[f,]), pch=15,col="black", lty=1, lwd=2,type="o", main=paste("nfeat=",f), 
           ylim=c(-0.01,0.35), ylab="sensitivity gain")
      lines(samplesizes, c(na,gain.sens.cv[f,]), pch=15, col="blue", lty=2, lwd=2,type="o")
      legend("topright", legend=c("training", "testing"), pch=15, col=c("black","blue"))
      grid()
    }
    
    for (f in c(1,3,7,9)) {
      plot(samplesizes, c(na,gain.percent.sens.resub[f,]), pch=15,col="black", lty=1, lwd=2,type="o", main=paste("nfeat=",f), 
           ylim=c(-5,50), ylab="sensitivity relative gain(%)")
      lines(samplesizes, c(na,gain.percent.sens.cv[f,]), pch=15, col="blue", lty=2, lwd=2,type="o")
      legend("topright", legend=c("training", "testing"), pch=15, col=c("black","blue"))
      grid()
    }
    
    
  }
  dev.off()
  save.image(file=paste(h0_strin, ".all.data.rdata", sep=""))
  
  
  
  
  
  out <- list()
  out$aroc <- aroc
  out$sensitivity.cv <- sensitivity.cv
  out$frp.cv <- fpr.cv
  out$sensitivity.resub <- sensitivity.resub
  out$frp.resub <- fpr.resub
  out
  
}


est_partition_boot <- function (node, samplesizes, nboot, nfold, strin) {
  
  h0_gt <- node$groupingpc
  h0_class <- node$groupingpc
  h0_zsel35 <- node$zsel35
  grp0_feat_pcloading <- colnames(node$zsel35)
  nallfeat <- ncol(node$zsel35)
  h0_strin <- strin 
  
  
  print("###############################")
  mycase=""
  if(mycase=="node0") {
    h0_gt <- node0$groupingpc
    h0_class <- node0$groupingpc
    h0_zsel35 <- node0.zsel35
    grp0_feat_pcloading <- colnames(node0.zsel35)
    nallfeat <- ncol(node0.zsel35)
    samplesizes <- seq(10,50,5)
    nboot <- 100
    nfold <- 4
    h0_strin <- 'vera.node0' 
  }
  if (mycase== "noder") {
    h0_gt <- noder$groupingpc
    h0_class <- noder$groupingpc
    h0_zsel35 <- noder.zsel35
    grp0_feat_pcloading <- colnames(noder.zsel35)
    nallfeat <- ncol(noder.zsel35)
    samplesizes <- seq(10,50,5)
    nboot <- 100
    nfold <- 4
    h0_strin <- 'vera.noder'
  }
  if (mycase== "noderl") {
    h0_gt <- noderl$groupingpc
    h0_class <- noderl$groupingpc
    h0_zsel35 <- noderl.zsel35
    grp0_feat_pcloading <- colnames(noderl.zsel35)
    nallfeat <- ncol(noderl.zsel35)
    samplesizes <- seq(10,50,5)
    nboot <- 100
    nfold <- 4
    h0_strin <- 'vera.noderl'
  }
  
  outfn <- paste(h0_strin, ".run.morecells.pdf", sep="")
  pdf(outfn)
  
  nf = 10
  ns = length(samplesizes)
  thrs = seq(1.4,1.6,0.05)
  
  sensitivity.cv <- fpr.cv <- array(0, dim=c(length(thrs), nf, ns))
  sensitivity.resub <- fpr.resub <- array(0, dim=c(length(thrs), nf, ns))
  sd.sensitivity.cv <- sd.fpr.cv <- array(0, dim=c(length(thrs), nf, ns))
  sd.sensitivity.resub <- sd.fpr.resub <- array(0, dim=c(length(thrs), nf, ns))
  aroc <- sd.aroc <- matrix(0,nrow=ns,ncol=nf)
  
  
  for (iiii in 1:length(thrs)) {
    i.thr <- thrs[iiii]
    print(paste("ithr=",i.thr))
    h0_pcloading.cv4roc <- mystepbystepboot(h0_gt, h0_class, h0_zsel35, grp0_feat_pcloading, nallfeat, samplesizes, nboot, nfold, i.thr)
    
    for (sss in 1:ns) {
      sss.samplesize = samplesizes[sss]
      sss.str <- paste('samplesize', sss.samplesize, sep="")
      
      sensitivity.cv[iiii, 1:nf,sss] <- h0_pcloading.cv4roc[[sss.str]]$spcsen.cv[2,1:nf]
      fpr.cv[iiii, 1:nf,sss] <- h0_pcloading.cv4roc[[sss.str]]$spcsen.cv[1,1:nf]
      sensitivity.resub[iiii, 1:nf,sss] <- h0_pcloading.cv4roc[[sss.str]]$spcsen.resub[2,1:nf]
      fpr.resub[iiii, 1:nf,sss] <- h0_pcloading.cv4roc[[sss.str]]$spcsen.resub[1,1:nf]
      
      sd.sensitivity.cv[iiii, 1:nf,sss] <- h0_pcloading.cv4roc[[sss.str]]$sd.spcsen.cv[2,1:nf]
      sd.fpr.cv[iiii, 1:nf,sss] <- h0_pcloading.cv4roc[[sss.str]]$sd.spcsen.cv[1,1:nf]
      sd.sensitivity.resub[iiii, 1:nf,sss] <- h0_pcloading.cv4roc[[sss.str]]$sd.spcsen.resub[2,1:nf]
      sd.fpr.resub[iiii, 1:nf,sss] <- h0_pcloading.cv4roc[[sss.str]]$sd.spcsen.resub[1,1:nf]
    }
  }
  dev.off()
  
  pdf(paste(h0.strin, ".bootstrap.summary.nb100.nf4.relgain.0.pdf", sep=""), width=8, height=8)
  for (s in 1:ns) {
    for (f in 1:nf) {
      aroc[s,f] <- cal_aroc(fpr.cv[,f,s], sensitivity.cv[,f,s])
      sd.aroc[s,f] <- mean(c(sd.fpr.cv[,f,s], sd.sensitivity.cv[,f,s]))
    }
    
    plot(fpr.cv[,1,s], sensitivity.cv[,1,s], col="grey", pch=15, ylab="sensitivity", xlab="false positive rate = 1-specificity", main=paste("roc of pvalb call :", samplesizes[s]), xlim=c(0,0.5),ylim=c(0.5,1), type="o", lty=1, lwd=3)
    for (f in 2:8) {
      lines(fpr.cv[,f,s], sensitivity.cv[,f,s], col=f, pch=15, type="o", lty=1, lwd=3)
    }
    lines(fpr.cv[,9,s], sensitivity.cv[,9,s], col=9, pch=22, type="o", lty=2,lwd=2)
    lines(fpr.cv[,10,s], sensitivity.cv[,10,s], col=10, pch=22, type="o", lty=2,lwd=2)
    thisstr <- paste(signif(aroc[s,],3),":", 1:nf, "features")
    thiscol <- 1:nf
    grid()
    legend("bottomright", legend=thisstr, col=thiscol, pch=c(rep(15,8),rep(22,2)) )
  }
  
  for (f in 1:nf) {
    if (f==1) {
      plot(samplesizes, aroc[,f], col=f, pch=15, ylab="area of roc", xlab="pvalb sample size",
           main="area of roc vs. sample size : pvalb", type="o", lwd=2, ylim=c(0.5,1.1))
      arrows(samplesizes, y0=aroc[,f]-sd.aroc[,f], y1=aroc[,f]+sd.aroc[,f], code=3, angle=90, length=.05, col=f)
      
    } else {
      if (f > 8) {
        lines(samplesizes, aroc[,f], col=f, pch=22, type="o", lwd=2,lty=2)
        arrows(samplesizes, y0=aroc[,f]-sd.aroc[,f], y1=aroc[,f]+sd.aroc[,f], code=3, angle=90, length=0.1, col=f)
      } else {
        lines(samplesizes, aroc[,f], col=f, pch=15, type="o", lwd=2,lty=1)
        arrows(samplesizes, y0=aroc[,f]-sd.aroc[,f], y1=aroc[,f]+sd.aroc[,f], code=3, angle=90, length=0.1, col=f)
      }
    }
    grid()
    thisstr <- paste(1:nf, "features")
    thiscol <- 1:nf
    legend("bottomright", legend=thisstr, col=thiscol, pch=c(rep(15,8),rep(22,2)) )
  }
  
  tmp=1
  if (tmp==1) {
    
    sen.cv.mean <-sensitivity.cv[2,,]
    sen.resub.mean <-sensitivity.resub[2,,]
    sen.cv.sd <-sd.sensitivity.cv[2,,]
    sen.resub.sd <-sd.sensitivity.resub[2,,]
    
    mean.sens.cv <-sensitivity.cv[2,,]
    mean.sens.resub <-sensitivity.resub[2,,]
    sd.sens.cv <-sd.sensitivity.cv[2,,]
    sd.sens.resub <-sd.sensitivity.resub[2,,]
    
    par(mfrow=c(2,2))
    for (f in c(1,3,7,9)) {
      plot(samplesizes, mean.sens.resub[f,], pch=15,col="black", lty=1, lwd=2,type="o", main=paste("nfeat=",f), ylim=c(0.5,1), ylab="sensitivity")
      arrows(samplesizes,y0=mean.sens.resub[f,]-(0.316*sd.sens.resub[f,]), y1=mean.sens.resub[f,]+(0.316*sd.sens.resub[f,]),
             code=3,angle=90, length=0.03, col="black")
      
      lines(samplesizes, mean.sens.cv[f,], pch=15, col="blue", lty=2, lwd=2,type="o")
      arrows(samplesizes,y0=mean.sens.cv[f,]-(0.316*sd.sens.cv[f,]), y1=mean.sens.cv[f,]+(0.316*sd.sens.cv[f,]),
             code=3,angle=90, length=0.03, col="blue")
      legend("bottomright", legend=c("training", "testing"), pch=15, col=c("black","blue"))
      grid()
    }
    
    gain.sens.resub <- t(apply(mean.sens.resub, 1, diff))
    gain.cum.sens.resub <- t(apply(gain.sens.resub, 1, cumsum))
    gain.percent.sens.resub <- signif(100*(gain.sens.resub / mean.sens.resub[,1:7]),2)
    
    gain.sens.cv <- t(apply(mean.sens.cv, 1, diff))
    gain.cum.sens.cv <- t(apply(gain.sens.cv, 1, cumsum))
    gain.percent.sens.cv <- signif(100*(gain.sens.cv / mean.sens.cv[,1:7]),2)
    
    par(mfrow=c(2,2))
    for (f in c(1,3,7,9)) {
      plot(samplesizes, c(na,gain.sens.resub[f,]), pch=15,col="black", lty=1, lwd=2,type="o", main=paste("nfeat=",f), 
           ylim=c(-0.01,0.35), ylab="sensitivity gain")
      lines(samplesizes, c(na,gain.sens.cv[f,]), pch=15, col="blue", lty=2, lwd=2,type="o")
      legend("topright", legend=c("training", "testing"), pch=15, col=c("black","blue"))
      grid()
    }
    
    for (f in c(1,3,7,9)) {
      plot(samplesizes, c(na,gain.percent.sens.resub[f,]), pch=15,col="black", lty=1, lwd=2,type="o", main=paste("nfeat=",f), 
           ylim=c(-5,50), ylab="sensitivity relative gain(%)")
      lines(samplesizes, c(na,gain.percent.sens.cv[f,]), pch=15, col="blue", lty=2, lwd=2,type="o")
      legend("topright", legend=c("training", "testing"), pch=15, col=c("black","blue"))
      grid()
    }
    
    
  }
  dev.off()
  save.image(file=paste(h0_strin, ".all.data.rdata", sep=""))
  
  out <- list()
  out$aroc <- aroc
  out$sensitivity.cv <- sensitivity.cv
  out$frp.cv <- fpr.cv
  out$sensitivity.resub <- sensitivity.resub
  out$frp.resub <- fpr.resub
  out
  
}


