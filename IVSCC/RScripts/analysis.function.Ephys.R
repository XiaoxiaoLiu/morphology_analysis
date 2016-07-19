#Sys.setenv(R_SESSION_TMPDIR="/local2/tmp")
#library(multtest)
source("/data/informatics/changkyul/Ephys/Ephys_Unix/Hotelling/R/hotelling.test.r")
#source("/data/informatics/changkyul/Ephys/Ephys_Unix/clValid/R/clValid-Classes.R")
#source("/data/informatics/changkyul/Ephys/Ephys_Unix/clValid/R/clValid-functions.R")
#source("/data/informatics/changkyul/Ephys/Ephys_Unix/clValid/R/clValid-internal.R")
#source("/data/informatics/changkyul/Ephys/Ephys_Unix/clValid/R/clValid-Methods.R")

library( "FisherEM" )
library( gtools )
library( matrixcalc )
library( WGCNA )
library( multtest)



rbg = colorRampPalette(c("green", "black", "red"))(100)
bwr = colorRampPalette(c("black", "white", "red"))(100)
ryb = colorRampPalette(c("blue", "yellow", "red"))(100)
iryb = colorRampPalette(c("red", "yellow", "blue"))(100)
iryw = colorRampPalette(c("white", "yellow", "red"))(100)
bwy = colorRampPalette(c("darkslateblue", "white", "orange2"))(100)
ryw = heat.colors(100)
wyr = rev(ryw)
w2b = colorRampPalette(c("white", "blue"))(100)

add_sampleid <- function( FN, FNOUT, MappingFN="/data/informatics/changkyul/Ephys/Data/810SampleNameSampleID.table.csv" ) {

    datain <- read.csv(FN, header=TRUE)
    samplename <- datain[,1]

    mapping <- read.csv(MappingFN, header=TRUE)
    tmp1 <- get.field(as.character(mapping[,"samplename"]), ";", 1)
    tmp2 <- get.field(as.character(mapping[,"samplename"]), ";", 2)
    tmp22 <- get.field(tmp2, "-", 2)
    mapping.samplename <- paste(tmp1, tmp22, sep=";")

    specimenid <- as.character(mapping[match(samplename, mapping.samplename), "sampleid"])
    sampleid <- paste(get.field(samplename, ";", 1), specimenid, sep=";")
    write.csv(data.frame(sampleid=sampleid, datain), file=FNOUT, row.names=FALSE)

}

add_samplename <- function( FN, FNOUT, MappingFN="/data/informatics/changkyul/Ephys/Data/810SampleNameSampleID.table.csv" ) {

    datain <- read.csv(FN, header=TRUE)
    sampleid <- get.field(as.character(datain[,1]), ";", 2)
    
    mapping <- read.csv(MappingFN, header=TRUE)
    samplename <- as.character(mapping[match(sampleid, mapping[,"sampleid"]), "samplename"])
    write.csv(data.frame(samplename=samplename, datain), file=FNOUT, row.names=FALSE)
}

get_crecolor_from_samplename <- function(samplename, datastr="09102015") {

   Sample.creall <- get.field(samplename, ";", 1)
   creall <- gsub("-Cre-D", "-Cre", Sample.creall)
    
   DATADir = "/data/informatics/changkyul/Ephys/Data/"
   creline_table = read.csv(paste(DATADir, "Cre_line_typeC.", datastr, ".csv", sep=""))
   creC <- as.character(creline_table[,"cre_line"])
   typeC <- substr(as.character(creline_table[,"type"]),1,3)
   pchC <- as.numeric(creline_table[,"pch1"])
   colorC <- as.character(creline_table[,"color"])

print(creC)

   print("Setting up Color & Symbols")
   mycolor <- rainbow(12)
   outcolor <- matrix(mycolor[match(creall,creC)], ncol=1)
   rownames(outcolor) <- samplename

   outcolor

}


add_meta_color <- function (mycrecolor) {

   meta.table.all <- read.csv("/data/informatics/changkyul/Ephys/Data/metadata.09252015.csv", header=TRUE, colClasses="character")
   #meta.table <- read.csv("\\\\Aibsdata/informatics/changkyul/Ephys/Data/metadata.09102015.csv", header=TRUE, colClasses="character")
   selected <- c(5:8,11:17)
   meta.names <- colnames(meta.table.all)[selected]
   tmpsamplename <- as.character(meta.table.all[,"name"])
   id <- as.character(meta.table.all[,"id"])
   cre <- as.character(meta.table.all[,"cre_line"])
   meta.samplename.all <- paste(cre, id, sep=";")
   
   idx <- match(rownames(mycrecolor), meta.samplename.all)
   meta.table <- meta.table.all[idx,]
   meta.samplename <- meta.samplename.all[idx]

   meta.value <- matrix(0, nrow=nrow(meta.table), ncol=length(selected))
   colnames(meta.value) <- meta.names
   rownames(meta.value) <- meta.samplename
   
   print(dim(mycrecolor))
   meta.color <- mycrecolor

   for (i in 1:length(selected)) {
        if (i==2) {
            itmp <- as.numeric(meta.table[, selected[i]])
            mask.missing <- which(is.na(itmp))
        } else {
            itmp <- as.character(meta.table[, selected[i]])
            mask.missing <- which(sapply(itmp, nchar)==0)
        }
        itmp[mask.missing] <- NA
        meta.value[,i] <- itmp
        
        jtmp <- match(itmp, sort(unique(itmp)))
        jN <- length(unique(jtmp))
        print(paste(meta.names[i], jN))
   
        if (i%%2==0) { 
            mycolormap <- topo.colors(jN+1)[sample(jN+1)] 
        } else { 
            #mycolormap <- cm.colors(jN+1)  
            mycolormap <- colorRampPalette(c("black", "white", "red"))(jN+1)[sample(jN+1)]
        }
        if (i==8) { mycolormap <- grey(c(0,1)) }
        mycolor <- mycolormap[jtmp] 
        mycolor[is.na(meta.value[,i])] <- "grey"
        mycolor1 <- matrix(mycolor, nrow=length(mycolor), ncol=1) 
   
        meta.color <- cbind(meta.color, mycolor1) 
   }
   colnames(meta.color) <- c("Assigned_Label", "Cre-line", meta.names)
   rownames(meta.color) <- meta.samplename

   meta.color
}

get_cluster_specific_features <- function (OUTDIR0, XXX, Id.id, pthr=0.01) { 

    idx <- match(rownames(XXX), names(Id.id))
    XXX.ClusterID <- Id.id[idx]
    unique.ClusterID <- sort(unique(XXX.ClusterID))
    XXX.label <- match(XXX.ClusterID, unique.ClusterID)
    unique.label <- unique(XXX.label)
    
    print("One Cluster vs. All Others")
    
    foi <- matrix("", nrow=50, ncol=length(unique.label))
    colnames(foi) <- unique.ClusterID
    allfoi <- c()
    allsoi <- c()
    allcol <- c()
    for (lbl in c(unique.label)) {

       lll <- as.numeric(XXX.label == lbl) 
       allsoi <- c(allsoi, which(XXX.label==lbl))

       myt.p <- standardScreeningBinaryTrait(XXX, lll)[, "pvalueStudent.0.vs.1"]
       tmp <- mt.rawp2adjp(myt.p,proc='TSBH')
       myt.adjp <- tmp$adjp[order(tmp$index),2]
    
       myt.adjp.order <- order(myt.adjp, decreasing=FALSE)
       myt.adjp.idx <- myt.adjp.order[which(myt.adjp[myt.adjp.order] < pthr)]
    
       DDD <- XXX[, myt.adjp.idx]
       fiti <- lm(lll ~ DDD - 1) 
       summary.fiti <- summary(fiti)
       t.pval <- summary.fiti[[4]][,4]
       t.pval.order <- order(t.pval, decreasing=FALSE)
       t.featname <- gsub("DDD", "", names(t.pval[t.pval.order])[t.pval[t.pval.order] < pthr])
       if (sum(t.pval<pthr) > 0) {
           foi[1:sum(t.pval<pthr),lbl] <- t.featname
           allfoi <- c(allfoi, t.featname)
           print(paste(lbl, sum(lll), length(myt.adjp.idx), sum(t.pval<pthr)))
       }
       allcol <- c(allcol, rep(lbl, sum(t.pval<pthr)))
       
       
    }

    out <- data.frame(label=allcol, featname=allfoi)
    write.csv(data.frame(label=allcol, featname=allfoi), file=paste(OUTDIR0, "/LDA.oneLabel.vs.others.p", pthr, ".csv", sep=""))

    save(allsoi, allfoi, foi, allcol, file=paste(OUTDIR0, "/LDA.oneLabel.vs.others.p", pthr, ".Rdata", sep=""))
    out   

}


gather_AssignedID_plotHeatmap_tree <- function(CLSDir, infolder, pthr, Node0, leg, filekey, flagPlot=TRUE) {

   print(names(leg))
   print(names(Node0))
   print('###########################################################################')
   print('######## PUT TOGETHER RESULTS  ###########') 
   print('###########################################################################')


   CLSfolder = paste(CLSDir, infolder, sep="/")
   print(CLSfolder)
   filelist <- dir(CLSfolder)   
   foi <- filelist[grep("AssignedID.csv", filelist)]

   ncntr <- 0
   outSampleName <-c()
   outID <-c()
   for (n in c(1:length(foi))) {
       nfname <- paste(CLSfolder, foi[n], sep="/")
       nresult <- read.csv(paste(CLSfolder, foi[n], sep="/"), header=TRUE)
       if (nrow(nresult)>0) {
           nsample <- as.character(nresult[,"samplename"])
           nid <- as.character(nresult[,"clusterID"])

           outSampleName <- c(outSampleName, nsample)
           outID <- c(outID, nid)
       }   
   }
   print(paste("number of cell ", length(outID))) 
   print(paste("number of unique cluster id's", length(unique(outID))) )
   print("table of cluster id's")
   ttt=table(outID)
   out <- data.frame(samplename=outSampleName, clusterID=outID)

   names(outID) <- outSampleName
   IDcolor <- grey(seq(0,1,length(unique(outID))))[match(outID, sample(unique(outID), length(unique(outID))))]
   outFN <- paste(CLSfolder, '/Final.', infolder, '.', pthr, length(ttt), '.csv', sep="")
   write.csv(data.frame(samplename=outSampleName, clusterID=outID), file=paste(CLSfolder, '/Final.', infolder, '.', pthr, length(ttt), '.csv', sep=""), row.names=FALSE)
   write.csv(ttt, file=paste(CLSfolder, '/Final.', infolder, '.', length(ttt), '.table.csv', sep=""))
   print("Cluster ID's are assigned by Final Nodes") 
   print(paste("Output File :", CLSfolder, '/Final.', infolder, '.', pthr, length(ttt), '.csv', sep=""))

   if (flagPlot) {
       print("Heatmap with Terminal Node Cluster Assignement")
       X35 <- Node0$ZSel35
       idx <- match(out$samplename, rownames(X35))

       cluster_feat <- get_cluster_specific_features(CLSfolder, X35, outID, pthr)  
       gdx <- match(unique(cluster_feat$featname), colnames(X35))
    
       pdf(paste(CLSfolder, '/FinalNode.Heatmap.ClusterSpecificFeature.', filekey, '.', pthr, '.0.pdf', sep=""), height=10, width=20)
       color.2row <- Node0$hybrid.crecolor2[, idx] 
       color.11row <- rbind(matrix(IDcolor, nrow=1),Node0$hybrid.crecolor2[2, idx]) 

       heatmap.3(t(X35[idx,gdx]), hclustfun=hclust, trace="none", Colv=FALSE, ColSideColors=color.2row, 
                 Rowv=TRUE, main=paste("Iterative PCA based Clustering \n with Cluster Specific Features, adjp<", pthr), 
                 keysize=0.8, margins=c(10,13), cexCol=0.5,cexRow=1, col=rbg)
       legend("bottomleft", bg="white", legend=leg$str2, pch=15, col=leg$crecolor2, cex=0.75) #, bty="n")

#   pdf(paste(CLSfolder, '/FinalNode.Heatmap.ClusterSpecificFeature.wMetadata.', filekey, '.', pthr, '.1.pdf', sep=""), height=10, width=20)
       print("MetaData")

       assigned.cls <- match(outID, sort(unique(outID)))

       Ncls <- length(unique(outID))
       #   label.color <- colorRampPalette(c("black", "white", "red"))(Ncls)[c(seq(1,Ncls,2), seq(2,Ncls,2))]
       label.color <- topo.colors(Ncls)[c(seq(1,Ncls,2), seq(2,Ncls,2))][sample(Ncls)]

       assigned.color.cre <- get_crecolor_from_samplename(outSampleName)
       print(dim(assigned.color.cre))
       assigned.color.label <- matrix(label.color[assigned.cls], ncol=1)
       print(dim(assigned.color.label))
    
       mycolor2 <- matrix(cbind(assigned.color.label, assigned.color.cre), ncol=2)
       print(dim(mycolor2))
       rownames(mycolor2) <- outSampleName
       mycolor3 <- add_meta_color(mycolor2) 

       meta.samplename <- rownames(mycolor3)
       #plotDendroAndColors( hh, meta.color[idx,], main="Clustering by Ephys Feature : any grouping by meta.data?" )
       plotOrderedColors( match(meta.samplename, outSampleName), mycolor3, main="Clustering by Ephys Feature : any grouping by meta.data?" )
       dev.off()
   }
   outFN

}

meanLR <- function(X) {
     grp1 <- X[1,] >= X[2,]
     grp2 <- !grp1
     thr <- sum(grp1)/(sum(grp1)+sum(grp2))
     meanL <- mean(abs(X[1,grp1]-thr))
     meanR <- mean(abs(X[2,grp2]-thr))
     c(meanL/(meanL+meanR),meanR/(meanL+meanR))
}

get_meanMembership2 <- function( Node, LRstr ) {

   if (Node$terminal) {
       out <- paste(LRstr,'[', nrow(Node$Xsoi),']', sep="")
   } else {
       mean.membership <- meanLR(Node$membershiplrs_pc)
       print(mean.membership)
       Node_Left <- Node$Left
       Node_Left$mean.membership <- round(mean.membership[1],3)
       Node_Right <- Node$Right
       Node_Right$mean.membership <- round(mean.membership[2],3)
       Lstr <- paste(LRstr, "L", sep="")
       Rstr <- paste(LRstr, "R", sep="")
       #print(paste(Lstr, Node_Left$mean.membership), sep="=")
       #print(paste(Rstr, Node_Right$mean.membership), sep="=")
       out <- paste("(", get_meanMembership2(Node_Left, Lstr), ":", Node_Left$mean.membership,",",
                         get_meanMembership2(Node_Right, Rstr), ":", Node_Right$mean.membership,")", sep="")
   }
   out
}


find_idx_p_byTtest <- function(Zin, Gin, pthr) {

    if (1) {
        TPval <- apply(Zin, 2, apply_Ttest, Gin)
        tmp <- mt.rawp2adjp(TPval,'TSBH')
        TadjP <- tmp$adjp[order(tmp$index),2]
        out <- which(TadjP < pthr)
    } else {
        feat.chosen <- c( "ap_2_peak_to_tr_time_long_sq_3",    "ap_2_peak_to_tr_time_long_sq_4",   
 		"ap_2_width_long_sq_4",              "ap_2_width_ramp",                  
 		"ap_fast_tr_change_long_sq_1",       "ap_fast_tr_change_long_sq_2",      
 		"ap_fast_tr_change_long_sq_3",       "ap_fast_tr_change_long_sq_4",      
 		"ap_fast_tr_change_long_sq_5",       "ap_fast_tr_change_ramp",           
 		"ap_height_change_long_sq_5",        "ap_updown_stroke_change_long_sq_1",
 		"ap_updown_stroke_change_long_sq_2", "ap_updown_stroke_change_ramp",     
 		"ap_width_change_long_sq_2",         "ap_width_change_long_sq_3",        
 		"ap_width_change_long_sq_4",         "ap_width_change_ramp",             
 		"ap_1_fast_tr_depth_ramp",           "ap_1_latency_ramp",                
 		"ap_1_peak_ramp",                    "ap_1_peak_to_tr_time_ramp",        
 		"ap_1_thresh_to_peak_time_ramp",     "ap_1_thr_ramp",                    
 		"ap_1_updown_stroke_ratio_ramp",     "ap_1_width_ramp",                  
 		"ap_2_fast_tr_depth_ramp",           "ap_2_thr_ramp",                    
 		"ap_height_change_ramp",             "ap_thr_change_ramp",               
 		"first_isi_ramp")                   
        out <- match(feat.chosen, colnames(Zin)) 
    }
    out
}


kyle_fem <- function ( featall, grpLR, Ks=c(2,3,4), modelin="all", methodin="svd", initin="kmeans", startin=10, flag.disp=FALSE, flag.debug=FALSE, ...) {

if (flag.debug) {
  featall <- Zin
  Ks <- c(2,3,4)
  modelin <- "DkBk"
  methodin <- "svd"
  initin <- "kmeans" 
  startin <-10 
  flag.disp <- FALSE
}
  if (0) {
  idx.adjp0.01 <- find_idx_p_byTtest (featall, grpLR, 0.01)
  print(idx.adjp0.01)
  feat <- featall[, idx.adjp0.01]
  }

  feat <- featall
  outK234 <-  fem(Y = feat, K=Ks , model=modelin, method=methodin, init=initin, nstart=startin, disp=flag.disp)
  Nparam  <- 2
  Nsample <- dim(feat)[1] 
  Nfeat <- dim(feat)[2] 

  u <- outK234$U[,1]
  qq <-  feat %*% u
  v <- c(var(qq))
  mu <- mean(qq)
  uut <- u %*% t(u)

  MU <- matrix(u * mu, nrow=length(u),ncol=1)

  u_ <- diag(length(u)) - uut
  qq_ <-  feat %*% u_
  v_ <- var(c(qq_)) #(115*115)
  #v_ <- var(qq_) #(115*115)
  u_u_t <- u_ %*% t(u_)

  # covarance matrix
  S <- (uut * v) + (u_u_t * v_)

  setClass("FEM", representation(k="numeric", pvalnorm="numeric", icl="numeric", k1icl="numeric"))
  out <- new("FEM") 
  if (is.singular.matrix(S)) {
      print("S is singular")
      slot(out,"icl") <- outK234$icl 
      slot(out,"k1icl") <- -1000000000
      slot(out,"k") <- 1
      slot(out,"pvalnorm") <- 1 
  } else {
      Sinv <- qr.solve(S)
  
      X <- t(feat)

      log_lik <- -((Nfeat/2)*log(2*pi)) - log(det(S))*(1/2)
      for (i in c(1:Nsample)) {
           dx <- X[,i] - MU 
           log_lik <- log_lik - (((t(dx) %*% Sinv) %*% dx)*(1/2))
      }
    
      aic <- (-2*log_lik) + (2*Nparam)
      icl <- log_lik - log(Nsample)

      if (outK234$icl > icl) { tmpicl <- outK234$icl } else { tmpicl <- icl }
      slot(out, "icl") <- c(tmpicl) 
      slot(out, "k1icl") <- c(icl)
      if (outK234$icl > icl) { slot(out, "k") <- outK234$K; slot(out,"pvalnorm") <- 0.001 }
      else { slot(out, "k") <- 1; slot(out,"pvalnorm") <- 1 }
      print(paste("k= 1", "    ", icl))
      print(paste("k=",outK234$K, "    ", outK234$icl))
  }
  print("======================================")
  
  out
}







gather_ID <- function (CLSDir, infolder, Nsample=70) {

   print('###########################################################################')
   print('######## PUT TOGETHER RESULTS  ###########')
   print('###########################################################################')


   CLSfolder = paste(CLSDir, infolder, sep="/")
   print(CLSfolder)
   filelist <- dir(CLSfolder)
   foi <- filelist[grep("AssignedID.csv", filelist)]

   ncntr <- 0
   outSampleName <-c()
   outID <-c()
   for (n in c(1:length(foi))) {
       nfname <- paste(CLSfolder, foi[n], sep="/")
       nresult <- read.csv(paste(CLSfolder, foi[n], sep="/"), header=TRUE)
       if (nrow(nresult)>0) {
           nsample <- as.character(nresult[,"samplename"])
           nid <- as.character(nresult[,"clusterID"])

           outSampleName <- c(outSampleName, nsample)
           outID <- c(outID, nid)
       }
   }
   if (length(outID) != Nsample) {print(paste(infolder, "number of cell does not match", length(outID))) }
   write.csv(data.frame(samplename=outSampleName, clusterID=outID), file=paste(CLSDir, '/', infolder, '.', length(outID), '.csv', sep=""), row.names=FALSE)
}




get.field <- function (tissueid, key, col) {
        slab <- rep(0, length(tissueid))
        for (i in c(1:length(tissueid))) {
	    tmp <- unlist(strsplit(as.character(tissueid[i]), key, fixed=TRUE))
            slab[i] <- tmp[col]
        }
        slab
}
   
#if (version$major==3) {
#library(fastcluster)
#hhclust <- function (d, method="ward", members=NULL) {
#    if (method == "ward") {
#        message("The \"ward\" method has been renamed to \"ward.D\"; note new \"ward.D2\"")
#        method <- "ward.D"
#    }
#    METHODS <- c("single", "complete", "average", "mcquitty",
#        "ward.D", "centroid", "median", "ward.D2")
#    method <- pmatch(method, METHODS)
#    if (is.na(method))
#        stop("Invalid clustering method.")
#    if (method == -1)
#        stop("Ambiguous clustering method.")
#    dendrogram <- c(.Call(fastcluster, attr(d, "Size"), method,
#        d, members), list(labels = attr(d, "Labels"), method = METHODS[method],
#        call = match.call(), dist.method = attr(d, "method")))
#    class(dendrogram) <- "hclust"
#    return(dendrogram)
#}
#}
#
#
#if (version$major==2) {
hhclust <- function(d, ...) {
    out <- hclust(d, method="ward")
}
#}




##############################################################
##############################################################
cal_entropy <- function (grpLR) {

    tt <- table(grpLR)
    lbl <- names(tt)
    Nlbl <- length(lbl)
    Plbl <- as.vector(tt/sum(tt))     
 
    out <- 0
    for (i in 1:Nlbl) {  
         out <- out -(Plbl[i] * log2(Plbl[i])) 
    }
    out
}



calcSeparationN <- function(variables,groupvariable) {
     # find out how many variables we have
   Vw <- calcWithinGroupsCovarianceN(variables, groupvariable)
   Vb <- calcBetweenGroupsCovarianceN(variables, groupvariable)
   sep <- Vb/Vw
   print(paste("Vw=",signif(Vw,3),"Vb=",signif(Vb,3),"separation=Vb/Ww",signif(sep,3)))
   c(Vw, Vb, sep)
}

calcBetweenGroupsCovarianceN <- function(variableN, groupvariable) {
     # find out how many values the group variable can take
     groupvariable2 <- as.factor(groupvariable)
     levels <- levels(groupvariable2)
     numlevels <- length(levels)
     # calculate the grand means
     nfeat <- ncol(variableN)
     variableNmean <- apply(variableN,2,mean)
     # calculate the between-groups covariance
     Covb <- 0
     for (i in 1:numlevels)
     {
        leveli <- levels[i]
        levelidataN <- variableN[groupvariable==leveli,]
        meanN <- apply(levelidataN,2, mean)
        levelilength <- nrow(levelidataN)
        term1 <- levelilength
        for (j in 1:nfeat) {
             term1 <- term1 * (meanN[j] - variableNmean[j])
        }
        Covb <- Covb + term1
     }
     Covb <- Covb / (numlevels - 1)
     return(Covb[[1]])
}

calcWithinGroupsCovarianceN <- function(variableN,groupvariable) {
     # find out how many values the group variable can take
     groupvariable2 <- as.factor(groupvariable)
     levels <- levels(groupvariable2)
     numlevels <- length(levels)
     # get the covariance of variable 1 and variable 2 for each group:

     nfeat <- ncol(variableN)
     variableNmean <- apply(variableN,2,mean)

     Covw <- 0
     for (i in 1:numlevels)
     {
        leveli <- levels[i]
        levelidataN <- variableN[groupvariable==leveli,]
        meanN <- apply(levelidataN,2, mean)
        levelilength <- nrow(levelidataN)

        # get the covariance for this group:
        term1 <- 0
        for (j in 1:levelilength)
        {
             tmp <- 1
             for (k in 1:nfeat) {
                 tmp <- tmp * (levelidataN[j,k] - meanN[k])
             }
             term1 <- term1 + tmp
        }
        Cov_groupi <- term1 # covariance for this group
        Covw <- Covw + Cov_groupi
     }
     totallength <- nrow(variableN)
     Covw <- Covw / (totallength - numlevels)
     return(Covw[[1]])
}



cal_p_norm <- function (x, muin, sdin, N=1) {

   # print(c(x, muin, sdin, N))
    z=(x-muin)/(sdin/sqrt(N))
   # print(z)
   # print(1-pnorm(abs(z)))
   # print(2*(1-pnorm(x, mean=muin, sd=sdin/sqrt(N))))
    out <- 2*(1-pnorm(abs(z)))
    out
}

calcBetweenGroupsCovariance2 <- function(variable1,variable2,groupvariable) {
     # find out how many values the group variable can take
if (0) {
variable1 <- aa1
variable2 <- aa2
groupvariable <- grp
}
     groupvariable2 <- as.factor(groupvariable)
     levels <- levels(groupvariable2)
     numlevels <- length(levels)
     # calculate the grand means
     variable1mean <- mean(variable1)
     variable2mean <- mean(variable2)
     # calculate the between-groups covariance
     Covb <- 0
     for (i in 1:numlevels)
     {
        leveli <- levels[i]
        levelidata1 <- variable1[groupvariable==leveli]
        levelidata2 <- variable2[groupvariable==leveli]
        mean1 <- mean(levelidata1)
        mean2 <- mean(levelidata2)
        levelilength <- length(levelidata1)
        term1 <- (mean1 - variable1mean)*(mean2 - variable2mean)*(levelilength)
        Covb <- Covb + term1
     }
     Covb <- Covb / (numlevels - 1)
     return(Covb)
}

calcWithinGroupsCovariance2 <- function(variable1,variable2,groupvariable) {
     # find out how many values the group variable can take
if (0) {
variable1 <- aa1
variable2 <- aa2
groupvariable <- grp
}
     groupvariable2 <- as.factor(groupvariable)
     levels <- levels(groupvariable2)
     numlevels <- length(levels)
     # get the covariance of variable 1 and variable 2 for each group:
     Covw <- 0
     for (i in 1:numlevels)
     {
        leveli <- levels[i]
        levelidata1 <- variable1[groupvariable==leveli]
        levelidata2 <- variable2[groupvariable==leveli]
        mean1 <- mean(levelidata1)
        mean2 <- mean(levelidata2)
        levelilength <- length(levelidata1)
        # get the covariance for this group:
        term1 <- 0
        for (j in 1:levelilength)
        {
           term1 <- term1 + ((levelidata1[j] - mean1)*(levelidata2[j] - mean2))
        }
        Cov_groupi <- term1 # covariance for this group
        Covw <- Covw + Cov_groupi
        #print(c(i, term1, Covw))
     }
     totallength <- length(variable1)
     Covw <- Covw / (totallength - numlevels)
     return(Covw)
}

calcWithinGroupsVariance <- function(variable,groupvariable) {
     # find out how many values the group variable can take
     #groupvariable2 <- as.factor(groupvariable[[1]])
     groupvariable2 <- as.factor(groupvariable)
     levels <- levels(groupvariable2)
     numlevels <- length(levels)
     # get the mean and standard deviation for each group:
     numtotal <- 0
     denomtotal <- 0
     for (i in 1:numlevels)
     {
        leveli <- levels[i]
        levelidata <- variable[groupvariable==leveli]
        levelilength <- length(levelidata)
        # get the standard deviation for group i:
        sdi <- sd(levelidata)
        numi <- (levelilength - 1)*(sdi * sdi)
        denomi <- levelilength
        numtotal <- numtotal + numi
        denomtotal <- denomtotal + denomi
     }
     # calculate the within-groups variance
     Vw <- numtotal / (denomtotal - numlevels)
     return(Vw)
}

calcBetweenGroupsVariance <- function(variable,groupvariable) {
     # find out how many values the group variable can take
     #groupvariable2 <- as.factor(groupvariable[[1]])
     groupvariable2 <- as.factor(groupvariable)
     levels <- levels(groupvariable2)
     numlevels <- length(levels)
     # calculate the overall grand mean:
     grandmean <- mean(variable)
     # get the mean and standard deviation for each group:
     numtotal <- 0
     denomtotal <- 0
     for (i in 1:numlevels)
     {
        leveli <- levels[i]
        levelidata <- variable[groupvariable==leveli]
        levelilength <- length(levelidata)
        # get the mean and standard deviation for group i:
        meani <- mean(levelidata)
        sdi <- sd(levelidata)
        numi <- levelilength * ((meani - grandmean)^2)
        denomi <- levelilength
        numtotal <- numtotal + numi
        denomtotal <- denomtotal + denomi
     }
     # calculate the between-groups variance
     Vb <- numtotal / (numlevels - 1)
     Vb <- Vb[[1]]
     return(Vb)
}

calcSeparations <- function(variables,groupvariable) {
     # find out how many variables we have
     variables <- as.data.frame(variables)
     numvariables <- length(variables)
     # find the variable names
     variablenames <- colnames(variables)
     # calculate the separation for each variable
     Vwithin <- rep(0, numvariables)
     Vbetween <- rep(0, numvariables)
     Separation <- rep(0, numvariables)

     for (i in 1:numvariables)
     {
        variablei <- variables[i]
        variablename <- variablenames[i]
        Vw <- calcWithinGroupsVariance(variablei, groupvariable)
        Vb <- calcBetweenGroupsVariance(variablei, groupvariable)
        sep <- Vb/Vw
        print(paste("variable",variablename,"Vw=",signif(Vw,3),"Vb=",signif(Vb,3),"separation=Vb/Ww",signif(sep,3)))
        Vwithin[i] <- Vw
        Vbetween[i] <- Vb
        Separation[i] <- sep
     }
     OUT <- matrix(0, nrow=numvariables, ncol=3)
     rownames(OUT) <- variablenames
     colnames(OUT) <- c("withinV", "betweenV", "between_within")
     Vbetween
}

printMeanAndSdByGroup <- function(variables,groupvariable) {
     # find the names of the variables
     variablenames <- c(names(groupvariable),names(as.data.frame(variables)))
     # within each group, find the mean of each variable
     #groupvariable <- groupvariable[,1] # ensures groupvariable is not a list
     means <- aggregate(as.matrix(variables) ~ groupvariable, FUN = mean)
     names(means) <- variablenames
     #print(paste("Means:"))
     #print(means)
     # within each group, find the standard deviation of each variable:
     sds <- aggregate(as.matrix(variables) ~ groupvariable, FUN = sd)
     names(sds) <- variablenames
     #print(paste("Standard deviations:"))
     #print(sds)
     # within each group, find the number of samples:
     samplesizes <- aggregate(as.matrix(variables) ~ groupvariable, FUN = length)
     names(samplesizes) <- variablenames
     #print(paste("Sample sizes:"))
     #print(samplesizes)

     sds
}




impute.my <- function(data) {
    if (is.matrix(data)) {
        mean.allsample <- apply(data, 1, function(x) { tmpmean <- mean(x, na.rm=TRUE) ; tmpR <- runif(1,0,0.01) ; tmpmean <- tmpmean * tmpR} ) 
        mean.allsample[is.na(mean.allsample)] = 1
        my.imputed <- data
        for (i in c(1:nrow(data))) {
            i.mask <- which(is.na(data[i,]))
            if (length(i.mask) > 0) {
               my.imputed[i, i.mask] <- mean.allsample[i]
            }
        }
    } else {
        my.imputed <- data
        #mean.allsample <- mean(data, na.rm=TRUE)
        tmpmean <- mean(data, na.rm=TRUE) ; tmpR <- runif(1,0,0.01) ; tmpmean <- tmpmean * tmpR 
        mean.allsample <- tmpmean ; 
        my.imputed[is.na(data)] = mean.allsample
        my.imputed = matrix(my.imputed, ncol=1)
         
    }
    my.imputed
}

impute.mymy <- function(data) {
    if (is.matrix(data)) {
        mean.allsample <- apply(data, 1, function(x) { mean(x, na.rm=TRUE)} )
        my.imputed <- data
        for (i in c(1:nrow(data))) {
            i.mask <- which(is.na(data[i,]))
            if (length(i.mask) > 0) {
               my.imputed[i, i.mask] <- mean.allsample[i]
            }
        }
    } else {
        my.imputed <- data
        mean.allsample <- mean(data, na.rm=TRUE)
        my.imputed[is.na(data)] = mean.allsample
        my.imputed = matrix(my.imputed, ncol=1)
         
    }
    my.imputed
}

Cneighbor <- function(x, stdev) {
     OUT <- norm(as.matrix(x), type="F") < stdev
     OUT
}

mynorm <- function(x) {
     OUT <- norm(as.matrix(x), type="F")
     OUT
}

CVIdensity <- function(Xi,Xj,Mui, Muj, stdev) {
     Ni <- nrow(Xi)
     Nj <- nrow(Xj)
     Xij <- rbind(Xi,Xj)
     Muij <- apply(Xij, 2, mean)
     Nij <- nrow(Xij)

     Zi <- t(matrix(rep(Mui, Ni), ncol=Ni))
     Yi <- Xi  - Zi
     neighi  <- sum(apply(Xi  - t(matrix(rep(Mui, Ni), ncol=Ni)), 1, Cneighbor, stdev))
     neighj  <- sum(apply(Xj  - t(matrix(rep(Muj, Nj), ncol=Nj)), 1, Cneighbor, stdev))
     neighij <- sum(apply(Xij - t(matrix(rep(Muij, Nij), ncol=Nij)), 1, Cneighbor, stdev))
     OUT <- (neighij / max(neighi, neighj))
     print(c(neighi, neighj, neighij, OUT))
     OUT
}

calCVI <- function (Xin, CLS) {
     print("calculating Cluster_Index for X[nsample X nfeat]")
     unique.CLS <- unique(CLS)
     nCLS <- length(unique.CLS)
     sdCLS <- rep(0, nCLS)

     if (!is.matrix(Xin)) {
        X <- matrix(Xin, nrow=length(Xin), ncol=1)
     } else {
        X <- as.matrix(Xin)
     }
     print(dim(X))
     
     meanCLS <- matrix(0, nrow=nCLS, ncol=ncol(X))

     sdALL <- norm(as.matrix(apply(X, 2, sd)), type="F")
     print(paste("sdALL=", sdALL))
     for (c in 1:nCLS) {
          sdCLS[c] <- norm(as.matrix(apply(X[CLS==unique.CLS[c],], 2, sd)), type="F")
          meanCLS[c,] <- apply(X[CLS==unique.CLS[c],], 2, mean)
          print(paste("cluster", c, "sdCLS[c]=",sdCLS[c]))
     }
     ICV <- mean(sdCLS)/sdALL
     print(paste("ICV=", ICV))

     stdev <- mean(sdCLS)
     print(paste("stdev=", stdev))

     cvi <- 0
     for (i in 1:nCLS) {
          Xi <- X[CLS==unique.CLS[i],]  
          Mui <- meanCLS[i,]
          for (j in setdiff(c(1:nCLS),i)) {
               Xj <- X[CLS==unique.CLS[j],]  
               Muj <- meanCLS[j,]
               Muij <- (Mui + Muj)/2
               cvi <- cvi + CVIdensity(Xi,Xj,Mui,Muj,stdev)
         }
     }
     ICD <- cvi/(nCLS*(nCLS-1))
     ICV + ICD
}


get_resub_err <- function(fiti, y) {
   est01 = rep(0, length(y))
   est01[fiti$fitted.values>0.5] <- 1
   clserr=100*sum(abs(y-est01))/length(y)
   clserr
}
get_cv_err <- function(fiti, TEin, yTE) {
   Iin  = matrix(1,length(yTE),1)
   est01 = rep(0, length(yTE))
   est.cv    = fiti$coefficients %*%  t(cbind(Iin, TEin))
   est01[est.cv>0.5] <- 1
   clserr=100*sum(abs(yTE-est01))/length(yTE)
   clserr
}

get_sen_spc <- function(est, gt, thr) {

   est01 <- rep(1,length(gt))
   mask <- est >= thr
   if (sum(mask, na.rm=TRUE)>0) {
       est01[est >= thr] <- 2
   }

   myct <- matrix(0,2,2)
   myct[1,1] <- sum(gt==1 & est01==1)
   myct[1,2] <- sum(gt==1 & est01==2)
   myct[2,1] <- sum(gt==2 & est01==1)
   myct[2,2] <- sum(gt==2 & est01==2)
#   print(myct)

   if (sum(myct[1,])==0) {
      spc <- NA
   } else {
      spc <- myct[1,1]/sum(myct[1,])
   }
   if (sum(myct[2,])==0) {
      sen <- NA
   } else {
      sen <- myct[2,2]/sum(myct[2,])
   }
   myOUT <- c(sen, spc)
   myOUT
}


getthrPCA1 <- function (pca0, grp12) {
    pcin <- pca0$x[,1] 
    pcin1 <- pcin[grp12==1]
    pcin2 <- pcin[grp12==2]
    if (mean(pcin1) < mean(pcin2)) {
        pcthr <- (max(pcin1) + min(pcin2))/2
    } else {
        pcthr <- (max(pcin2) + min(pcin1))/2
    }
    pcthr
}


getthrPCAn <- function (pca0, grp12) {
    pcin <- pca0$x[,1] 
    pcin1 <- pcin[grp12==1]
    pcin2 <- pcin[grp12==2]
    if (mean(pcin1) < mean(pcin2)) {
        pcthr <- (max(pcin1) + min(pcin2))/2
    } else {
        pcthr <- (max(pcin2) + min(pcin1))/2
    }
    pcthr
}

applyPCA1 <- function(pca0, pcthr0, Xin, GTin) {
   PCin <- Xin %*% pca0$rotation[,1]
   estin <- PCin > pcthr0

#   if (sum(estin==1) > sum(estin==0)) {
#       estin <- 1-estin
#   }
#   if (sum(GTin==2) > sum(GTin==1)) {
#       GTin <- 3 - GTin
#   }

   gt <- GTin
   est01 <- estin + 1 
   myct <- matrix(0,2,2)
   myct[1,1] <- sum(gt==1 & est01==1)
   myct[1,2] <- sum(gt==1 & est01==2)
   myct[2,1] <- sum(gt==2 & est01==1)
   myct[2,2] <- sum(gt==2 & est01==2)

   if (sum(myct[1,])==0) {
      spc <- NA
   } else {
      spc <- myct[1,1]/sum(myct[1,])
   }
   if (sum(myct[2,])==0) {
      sen <- NA
   } else {
      sen <- myct[2,2]/sum(myct[2,])
   }
   myOUT <- c(sen, spc)
   print(myOUT)
   myOUT

}


cal_aroc <- function (fpr_in, sen_in) {
   fpr <- c(0, rev(fpr_in), 1) 
   sen <- c(0, rev(sen_in), 1)
   aroc <- 0
   NN <- length(sen)
   for (i in (2:NN)) {
      dx <- abs(fpr[i] - fpr[i-1])
      dy <- abs(sen[i] + sen[i-1]) / 2
      darea <- dx*dy
      aroc <- aroc + darea 
      #print(signif(c(fpr[i], sen[i],dx,dy,darea,aroc),3))
   }
   signif(aroc,4)
}



plot_mean_std <- function (xval, Xin, ...) {
    if (length(xval) == nrow(Xin)) key=1 
    if (length(xval) == ncol(Xin)) key=2 
    Xmean <- apply(Xin, key, mean)
    Xsd <- apply(Xin, key, sd)
    Xplus <- Xmean + Xsd
    Xminus <- Xmean - Xsd
    plot(xval, Xmean, ...) 
    arrows(xval, y0=Xplus, y1=Xminus, length=0.025)
    grid()
}





plotMDS <- function (Xin, titlestr, mypch, leg.pch, crecolor, leg.crecolor, leg.str) {

    d <- as.dist(1-cor(t(Xin), use="pairwise.complete.obs"))
    fit <- cmdscale(d,eig=TRUE, k=3) # k is the number of dim
    par(fig=c(0,1,0,1),las=2, cex=0.75)
    par(mfrow=c(2,2))
    x <- fit$points[,1]
    y <- fit$points[,2]
    z <- fit$points[,3]
    
    par(mfrow=c(2,2))
    plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
         main=titlestr, pch=mypch, col=crecolor) #, xlim=c(-1,1), ylim=c(-1,1))
    grid()
    plot(x, z, xlab="Coordinate 1", ylab="Coordinate 3",
         main=titlestr, pch=mypch, col=crecolor) #, xlim=c(-1,1), ylim=c(-1,1))
    grid()
    plot(y, z, xlab="Coordinate 2", ylab="Coordinate 3",
         main=titlestr, pch=mypch, col=crecolor) #, xlim=c(-1,1), ylim=c(-1,1))
    grid()
    plot.new()
    legend('top', leg.str, pch=leg.pch, col=leg.crecolor, cex=1, bty="n" )
    
    par(mfrow=c(1,1))
    plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
         main=titlestr, pch=mypch, col=crecolor) #, xlim=c(-1,1), ylim=c(-1,1))
    grid()
    legend('topleft', leg.str, pch=leg.pch, col=leg.crecolor, cex=0.6 )
    rm(d,x,y,z,fit)
}

plotPCA <- function (Xin, titlestr, mypch, leg.pch, crecolor, leg.crecolor, leg.str, PDFFN, flagPlot=TRUE) {

    myPCA <- prcomp(Xin, cor=TRUE)
    mycumvar <- cumsum(myPCA$sdev^2)/sum(myPCA$sdev^2)
    myvarpercent <- signif(100*(myPCA$sdev^2/sum(myPCA$sdev^2)),2)
 
    x <- myPCA$x[,1]
    y <- myPCA$x[,2]
    z <- myPCA$x[,3]
    
    if (flagPlot) {
        par(fig=c(0,1,0,1),las=2, cex=0.75)
        par(mfrow=c(2,2), las=2)
        plot(x, y, xlab=paste("PC 1:", myvarpercent[1], "%"), ylab=paste("PC 2:", myvarpercent[2], "%"),
             main=titlestr, pch=mypch, col=crecolor) #, xlim=c(-1,1), ylim=c(-1,1))
        grid()
        plot(x, z, xlab=paste("PC 1:", myvarpercent[1], "%"),ylab=paste("PC 3:", myvarpercent[3], "%"),
             main=titlestr, pch=mypch, col=crecolor) #, xlim=c(-1,1), ylim=c(-1,1))
        grid()
    
        barplot(summary(myPCA)$importance[2,], xlab="Principal Component", ylab="Variance (%)", main='Variance Contribution')
        #plot(100*mycumvar, xlab="Principal Component", ylab="Cumulative Variance (%)", main='Variance Explained', pch=15, type="o", lty=2)
        grid()
    
        plot.new()
        legend('top', leg.str, pch=leg.pch, col=leg.crecolor, cex=0.80, bty="n" )
    }

    rm(x,y,z,mycumvar, myvarpercent)

    myPCA
}


clusterPCA <- function (Pin, NN1, crecolorL, leg.crecolorL, leg.strL, PDFFN, flagPlot=TRUE) {

    print(paste("====== clustering cells using", NN1, " PCs ======")) 
    #print(dim(crecolorL))
    dist.Pin <- dist(Pin$x[, 1:NN1])
    #print(dim(as.matrix(dist.Pin)))
    hpca.d <- hclust(dist.Pin, method="ward.D")

    if (NN1 > 1) {
        cor.Pin <- as.dist((1 + cor(t(Pin$x[, 1:NN1])))/2) 
        hpca.c <- hclust(cor.Pin, method="ward.D")
    }

    if (flagPlot) {
        pdf(PDFFN, width=12, height=12)
        heatmap.3(as.matrix(dist.Pin), Rowv=as.dendrogram(hpca.d), Colv=as.dendrogram(hpca.d), RowSideColors=t(crecolorL), ColSideColors=crecolorL, trace="none",
                  main=paste("Clustering of Cells using", NN1, "PCs, dist=Euclideian"), col=rbg, keysize=0.8)
        legend("bottomleft", legend=leg.strL, pch=15, col=leg.crecolorL, bg="white", cex=0.85)
    
        if (NN1>1) {
            heatmap.3(as.matrix(cor.Pin), Rowv=as.dendrogram(hpca.c), Colv=as.dendrogram(hpca.c), RowSideColors=t(crecolorL), ColSideColors=crecolorL, trace="none",
                      main=paste("Clustering of Cells using", NN1, "PCs, dist=1-correlation"), col=rbg, keysize=0.8)
            legend("bottomleft", legend=leg.strL, pch=15, col=leg.crecolorL, bg="white", cex=0.85)
        }
        dev.off()
    }
}


partPCA <- function(Pin, Zin, OUTFN, flagPlot) {
    meth <- "ward"
    sigval_min <- 1.0
    wsigval_min <- 1.0
   
    NN=5 
    if (ncol(Zin) < NN) { NN=ncol(Zin)-1}
    NN1=NN+1
    grp12 <- wgrp12 <- matrix(0, nrow(Zin), NN1)
    cindex <- wcindex <- rep(1, NN1)
    pnorm <- wpnorm <- rep(1,NN1)
    for (Npc in c(1:NN1)) {

         if (Npc<NN1) { 
             newdistmat<-matrix(0,nrow=nrow(Pin$x),ncol=nrow(Pin$x));
             coordmat <- matrix(0,nrow=nrow(Pin$x),ncol=Npc)
             wnewpin <- matrix(0,nrow=nrow(Pin$x),ncol=Npc)
             newpin <- matrix(0,nrow=nrow(Pin$x),ncol=Npc)

             for (ii in 1:Npc) {
                  tempdist<-outer(Pin$x[,ii],Pin$x[,ii],"-");
                  if (meth=='ward') {
                      newdistmat<-newdistmat+(abs(tempdist)*(summary(Pin)$importance[2,ii]))^2;  
                  } else { 
                      newdistmat<-newdistmat+abs(tempdist)*(summary(Pin)$importance[2,ii]);
                  }
    
                  coordmat[,ii]<-Pin$x[,ii]*summary(Pin)$importance[2,ii];
                  newpin[,ii] <- Pin$x[,ii] 
                  wnewpin[,ii] <- scale(Pin$x[,ii])
             }
         } else {
             newdistmat <- dist(Zin) ;
             coordmat <- Zin ;
             newpin <- Zin ;
             wnewpin <- Zin ;
         }

         #hhhw=hclust(dist(wnewpin),method="ward")
         #hhhw_2grp <- cutree(hhhw,k=2)
         #hhh <- hclust(dist(newpin), method="ward")
         #hhh_2grp <- cutree(hhh, k=2)

         hhhw=hclust(as.dist(newdistmat),method="ward")
         hhhw_2grp <- cutree(hhhw,k=2)
         hhh <- hclust(dist(coordmat), method="ward")
         hhh_2grp <- cutree(hhh, k=2)

         ntry <- 2*nrow(Zin)
         sigval0 <- sigclust(coordmat, ntry, icovest=2) 
         sigval <- sigclust(coordmat, ntry, label=as.vector(hhh_2grp), icovest=2) 
         sigvalw <- sigclust(coordmat, ntry, label=as.vector(hhhw_2grp), icovest=2) 

         if (sum(hhhw_2grp==1) < sum(hhhw_2grp==2)) {
             grp12[,Npc] <- 3 - hhh_2grp
             wgrp12[,Npc] <- 3 - hhhw_2grp
         } else {
print("0")
             grp12[,Npc] <- hhh_2grp
             wgrp12[,Npc] <- hhhw_2grp
         }
         cindex[Npc] <- sigval@xcindex
         wcindex[Npc] <- sigvalw@xcindex
         pnorm[Npc] <- -log10(sigval@pvalnorm)
         wpnorm[Npc] <- -log10(sigvalw@pvalnorm)
         #print(paste("sigclust=", c(-log10(sigval0@pvalnorm), -log10(sigval@pvalnorm), -log10(sigvalw@pvalnorm))))

          
         if( sigval@pvalnorm < sigval_min ) {
                  sigval_min <- sigval@pvalnorm
                  min_Npc <- Npc
                  min_2grp <- hhh_2grp
         }
         if( sigvalw@pvalnorm < wsigval_min ) {
                  wsigval_min <- sigvalw@pvalnorm
                  wmin_Npc <- Npc
                  wmin_2grp <- hhhw_2grp
         }
    }
    print(c(min_Npc, wmin_Npc))

    if (flagPlot) {
        pdf(OUTFN, width=10, height=6)
        par(mfrow=c(1,2))
        plot(1:NN, pnorm[1:NN], pch=1, main="pvalue", xlab="number of PC''s", ylab="-log10(Pval)"); grid() 
        plot(1:NN, cindex[1:NN], pch=1, main="C-index", xlab="number of PC''s", ylab="C-index" ); grid() 
        legend("bottomright", legend="C-index=sum(WithinCls.RSS)/Total.RSS")
        plot(1:NN, wpnorm[1:NN], pch=1, main="pvalue", xlab="number of PC''s", ylab="-log10(Pval)"); grid() 
        plot(1:NN, wcindex[1:NN], pch=1, main="C-index", xlab="number of PC''s", ylab="C-index"); grid() 
        legend("bottomright", legend="C-index=sum(WithinCls.RSS)/Total.RSS")
        dev.off()
    }

    Out <- list()

    Out$grp12 <- grp12
    Out$wgrp12 <- wgrp12

    Out$min_grp12 <- min_2grp
    Out$wmin_grp12 <- wmin_2grp
    Out$pnorm <- pnorm
    Out$wpnorm <- wpnorm
    Out$cindex <- cindex
    Out$wcindex <- wcindex

    Out$Npc <- min_Npc
    Out$wNpc <- wmin_Npc
    Out
}

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

plot_those <- function(SSS, idx0, idx1, titlestr) {
    S0 = SSS[idx0,]
    C0=rep(colnames(S0),length(idx0))
    Profile0=c(t(S0))
    S1 = SSS[idx1,]
    C1=rep(colnames(S1),length(idx1))
    Profile1=c(t(S1))

    par(mfrow=c(2,2), las=2)
    verboseBarplot.fixedorder(Profile0, C0,
    ylab="scaled",xlab="",main=paste(titlestr, 'subgroup L', sep=" : "),las=2, 
    cex=0.85, cex.main=0.75, cex.lab=0.75, cex.axis=0.75, color=3, 
    numberStandardErrors = 1, KruskalTest = FALSE, AnovaTest = FALSE, two.sided = FALSE, 
    ylim=c(0,1))
    grid()
    verboseBarplot.fixedorder(Profile1, C1,
    ylab="scaled",xlab="",main=paste(titlestr, 'subgroup R', sep=" : "),las=2, 
    cex=0.85, cex.main=0.75, cex.lab=0.75, cex.axis=0.75, color=2, 
    numberStandardErrors = 1, KruskalTest = FALSE, AnovaTest = FALSE, two.sided = FALSE, 
    ylim=c(0,1))
    grid()
}

plot_these <- function( inresubcv.mean, inresubcv.feat, ZZZ, idx0, idx1, legstr, titlestr, crecolor2, leg.str, leg.crecolorIE, nselfeat=4, XSel ) {

    par(mfrow=c(1,1))
    plot(1:10, inresubcv.mean[,1],pch=20,col='blue', type="o", main=paste("Selected Features with LDA : \n data-driven grouping", titlestr),
        ylim=c(0.0,max(inresubcv.mean[,2])+0.05), xlab="Number of Features", ylab='Fitting.Error')
    lines(1:10,inresubcv.mean[,2],pch=20,col='red', type="o")
    grid()
    legend("topright", legend=c("Resubstitution", paste(10, "fold Cross-validation")), pch=20, col=c("blue", "red"))

    selected.featnameshort <- inresubcv.feat[,ncol(inresubcv.feat)]
    par(mfrow=c(ceiling(nselfeat/2),2))
    for (i in 1:nselfeat) {
       i.feat = selected.featnameshort[i]
       iZ.0 = ZZZ[idx0, i.feat]
       iZ.1 = ZZZ[idx1, i.feat]
       ihc.0 = hist(iZ.0, seq(-3.5,3.5, 0.2), plot=FALSE)
       ihc.1 = hist(iZ.1, seq(-3.5,3.5, 0.2), plot=FALSE)
       ihc = cbind(ihc.1$counts, ihc.0$counts)
       rownames(ihc) <- as.character(ihc.0$mids)
       barplot(t(ihc), beside=TRUE, main=i.feat, col=c(1,2), border=NA)
       legend('topright', legend=legstr, cex=0.8, col=c(1,2), pch=15, bty="n")
    }

    par(mfrow=c(1,1))
    h = heatmap.3(t(ZZZ[, selected.featnameshort[1:nselfeat]]), Rowv=FALSE, hclustfun=hclust, trace="none", ColSideColors=crecolor2, #RowSideColors=crecolorIE, dendrogram="col", keysize=1,
              main=paste("Partition by", nselfeat, "Ephys Features : \n data-driven", titlestr), 
              margins=c(10,10), cexCol=0.5,cexRow=0.9, col=rbg)
    legend("bottomleft", bg="white", legend=leg.str, pch=15, col=leg.crecolorIE, cex=0.75) #, bty="n")

    pairs(XSel[, selected.featnameshort[1:nselfeat]], col=crecolor2)
    legend("bottomleft", bg="white", legend=leg.str, pch=15, col=leg.crecolorIE, cex=0.75) #, bty="n")
    pairs(ZZZ[, selected.featnameshort[1:nselfeat]], col=crecolor2)

}


zscore.mask <- function(x, mask) {
     x.na <- x
     nfeat <- ncol(x.na)
     for (i in 1:nfeat) {
          i.name <- colnames(x.na)[i]
          x.na[which(is.na(mask[,i.name])), i] <- NA
     }
     tx <- t(x)
     tx.na <- t(x.na)

     nc <- ncol(tx.na)
     mymean <- matrix(rep(apply(tx.na, 1, mean, na.rm=TRUE), nc), ncol=nc)
     mystd <- matrix(rep(apply(tx.na, 1, sd, na.rm=TRUE),nc), ncol=nc)
     myz <- (tx - mymean)/mystd
}

zscore <- function (x) {
     nc <- dim(x)[2]
     mymean <- matrix(rep(apply(x, 1, mean, na.rm=TRUE), nc), ncol=nc)
     mystd <- matrix(rep(apply(x, 1, sd, na.rm=TRUE),nc), ncol=nc)
     myz <- (x - mymean)/mystd

     myz[is.na(myz)] <- 0
     return(myz)
}


zscore1d <- function (x) {
   if (var(x)==0) { myz <- x } else {
   myz <- (x -  mean(x)) / sqrt(var(x))
   }
   myz[is.na(myz)] <- 0
   return(myz)
}


heatmap.2 <- function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
    distfun = dist, hclustfun = hclust, dendrogram = c("both",
        "row", "column", "none"), symm = FALSE, scale = c("none",
        "row", "column"), na.rm = TRUE, revC = identical(Colv,
        "Rowv"), add.expr, breaks, symbreaks = min(x < 0, na.rm = TRUE) ||
        scale != "none", col = "heat.colors", colsep, rowsep,
    sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote, notecex = 1,
    notecol = "cyan", na.color = par("bg"), trace = c("column",
        "row", "both", "none"), tracecol = "cyan", hline = median(breaks),
    vline = median(breaks), linecol = tracecol, margins = c(5,
        5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr),
    cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL,
    key = TRUE, keysize = 1.5, density.info = c("histogram",
        "density", "none"), denscol = tracecol, symkey = min(x <
        0, na.rm = TRUE) || symbreaks, densadj = 0.25, main = NULL,
    xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL, ...)
{
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale))
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col))
        col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none"))
        warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv))
        Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv))
        Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv))
        Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x))
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1)
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2)
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote))
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
            c("both", "row"))) {
            if (is.logical(Colv) && (Colv))
                dendrogram <- "column"
            else dedrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
            c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv))
                dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc)
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm)
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm)
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow))
        labRow <- if (is.null(rownames(x)))
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol))
        labCol <- if (is.null(colnames(x)))
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) <
        1) {
        if (missing(col) || is.function(col))
            breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks)
            breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function")
        col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei))
        lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid))
        lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)
        if (!missing(ColSideColors)) {
            if (!is.character(ColSideColors) || length(ColSideColors) !=
                nc)
                stop("'ColSideColors' must be a character vector of length ncol(x)")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] +
                1)
            lhei <- c(lhei[1], 0.2, lhei[2])
        }
        if (!missing(RowSideColors)) {
            if (!is.character(RowSideColors) || length(RowSideColors) !=
                nr)
                stop("'RowSideColors' must be a character vector of length nrow(x)")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) -
                1), 1), lmat[, 2] + 1)
            lwid <- c(lwid[1], 0.2, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }
    if (length(lhei) != nrow(lmat))
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat))
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    if (!missing(RowSideColors)) {
        par(mar = c(margins[1], 0, 0, 0.5))
        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    }
    if (!missing(ColSideColors)) {
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    }
    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr"))
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 +
        c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col,
        breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr"))
        retval$rowDendrogram <- ddr
    if (exists("ddc"))
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) {
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
            col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
        cex.axis = cexCol)
    if (!is.null(xlab))
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
        cex.axis = cexRow)
    if (!is.null(ylab))
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr))
        eval(substitute(add.expr))
    if (!missing(colsep))
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0,
            length(csep)), xright = csep + 0.5 + sepwidth[1],
            ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1,
            col = sepcolor, border = sepcolor)
    if (!missing(rowsep))
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) +
            1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) +
            1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1,
            col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol,
                  lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote))
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
            col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main))
        title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        par(mar = c(5, 4, 2, 1), cex = 0.75)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
        }
        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
            xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        axis(1, at = xv, labels = lv)
        if (scale == "row")
            mtext(side = 1, "Row Z-Score", line = 2)
        else if (scale == "column")
            mtext(side = 1, "Column Z-Score", line = 2)
        else mtext(side = 1, "Value", line = 2)
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
                lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("Color Key\nand Density Plot")
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
                col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title("Color Key\nand Histogram")
            par(cex = 0.5)
            mtext(side = 2, "Count", line = 2)
        }
        else title("Color Key")
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
        high = retval$breaks[-1], color = retval$col)
    invisible(retval)
}



heatmap.3 <- function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
    distfun = dist, hclustfun = hclust, dendrogram = c("both", 
        "row", "column", "none"), symm = FALSE, scale = c("none", 
        "row", "column"), na.rm = TRUE, revC = identical(Colv, 
        "Rowv"), add.expr, breaks, symbreaks = min(x < 0, na.rm = TRUE) || 
        scale != "none", col = "heat.colors", colsep, rowsep, 
    sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote, notecex = 1, 
    notecol = "cyan", na.color = par("bg"), trace = c("column", 
        "row", "both", "none"), tracecol = "cyan", hline = median(breaks), 
    vline = median(breaks), linecol = tracecol, margins = c(5, 
        5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr), 
    cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL, 
    key = TRUE, keysize = 1.5, density.info = c("histogram", 
        "density", "none"), denscol = tracecol, symkey = min(x < 
        0, na.rm = TRUE) || symbreaks, densadj = 0.25, main = NULL, 
    xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL, 
    ...) 
{
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale)) 
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col)) 
        col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none")) 
        warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv)) 
        Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv)) 
        Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv)) 
        Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2) 
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote)) 
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% 
            c("both", "row"))) {
            if (is.logical(Colv) && (Colv)) 
                dendrogram <- "column"
            else dedrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% 
            c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv)) 
                dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc) 
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }

    if (exists("hcc")) {retval$hcc <- hcc}
    if (exists("hcr")) {retval$hcr <- hcr}
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow)) 
        labRow <- if (is.null(rownames(x))) 
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol)) 
        labCol <- if (is.null(colnames(x))) 
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 
        1) {
        if (missing(col) || is.function(col)) 
            breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks) 
            breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
                length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function") 
        col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei)) 
        lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid)) 
        lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)
        if (!missing(ColSideColors)) {
            if (!is.character(ColSideColors)) #|| ncol(ColSideColors) != nc) 
                stop("'ColSideColors' must be a character ") #vector of length ncol(x)")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
            if (is.vector(ColSideColors)) nnn=1 
            else nnn=nrow(ColSideColors)
            lhei <- c(lhei[1], nnn*0.1, lhei[2])
        }
        if (!missing(RowSideColors)) {
            if (!is.character(RowSideColors)) #|| length(RowSideColors) != nr) 
                stop("'RowSideColors' must be a character ")
            if (is.vector(RowSideColors)) nnn=1 
            else nnn=ncol(RowSideColors)
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[, 2] + 1)
            lwid <- c(lwid[1], nnn*0.1, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }
    if (length(lhei) != nrow(lmat)) 
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat)) 
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    if (!missing(RowSideColors)) {
        par(mar = c(margins[1], 0, 0, 0.5))

        if (is.vector(RowSideColors)) {
            image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
        } 
        if (is.matrix(RowSideColors)) {
            jk.row = RowSideColors
            jk.xy = matrix(which(jk.row != "0"), dim(jk.row))
colnames(jk.xy) <- colnames(jk.row)
            #image(t(jk.xy), col = jk.row[rowInd, ], xaxt="n", yaxt="n")
            #grid(nx=ncol(jk.row), ny=nrow(jk.row), lty=1, col="black")
            image(t(jk.xy), col = jk.row[rowInd, ], xaxt="n", yaxt="n")
#            axis(3, at=seq(0,1,1/(ncol(jk.xy)-1)),labels=colnames(jk.xy), las=2, cex.axis = cexCol, tick=0)
            axis(1, at=seq(0,1,1/(ncol(jk.xy)-1)),labels=colnames(jk.xy), las=2, cex.axis = cexCol, tick=0)
         #   grid(nx=ncol(jk.row), ny=nrow(jk.row), lty=1, col="black")
        }

    }
    if (!missing(ColSideColors)) {
        par(mar = c(0.5, 0, 0, margins[2]))
        if (is.vector(ColSideColors)) {
            image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
        } 
        if (is.matrix(ColSideColors)) {
            jk.col = ColSideColors
            jk.xy = matrix(which(jk.col != "0"), dim(jk.col))
            image(t(jk.xy), col = jk.col[, colInd], axes = FALSE)
        }

    }
    par(mar = c(margins[1], 0, 0, margins[2]))
    if (!symm || scale != "none") {
        x <- t(x)
        cellnote <- t(cellnote)
    }
    if (revC) {
        iy <- nr:1
        if (exists("ddr")) 
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
        c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
        breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr")) 
        retval$rowDendrogram <- ddr
    if (exists("ddc")) 
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
   # if (!invalid(na.color) & any(is.na(x))) {
   #     mmat <- ifelse(is.na(x), 1, NA)
   #    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
   #         col = na.color, add = TRUE)
   # }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexCol)
    if (!is.null(xlab)) 
        mtext(xlab, side = 1, line = margins[1] - 1.25)
#    axis(3, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
#        cex.axis = cexCol)
    if (!is.null(xlab)) 
        mtext(xlab, side = 3, line = margins[1] - 1.25)

    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexRow)
    if (!is.null(ylab)) 
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr)) 
        eval(substitute(add.expr))
    if (!missing(colsep)) 
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, 
            length(csep)), xright = csep + 0.5 + sepwidth[1], 
            ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, 
            col = sepcolor, border = sepcolor)
    if (!missing(rowsep)) 
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
            1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
            1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
            col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol, 
                  lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote)) 
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
            col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main)) 
        title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        par(mar = c(5, 4, 2, 1), cex = 0.75)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x))
            tmpbreaks[length(tmpbreaks)] <- max(abs(x))
        }
        else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
        }
        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, 
            xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        axis(1, at = xv, labels = lv)
        if (scale == "row") 
            mtext(side = 1, "Row Z-Score", line = 2)
        else if (scale == "column") 
            mtext(side = 1, "Column Z-Score", line = 2)
        else { 
            #mtext(side = 1, "Value", line = 2)
            mtext(side = 1, "", line = 2)
        }
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
                lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("") ;#Color Key and Density Plot", cex=0.25)
            #title("Color Key and Density Plot", cex=0.25)
            par(cex = 0.25)
            mtext(side = 2, "", line = 2)
            #mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
                col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            #title("Color Key and Histogram", cex=0.25)
            title("", cex=0.25)
            par(cex = 0.25)
            mtext(side = 2, "", line = 2)
            #mtext(side = 2, "Count", line = 2)
        }
        else { title("", cex=0.25)
               #title("Color Key", cex=0.25)
        }
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], 
        high = retval$breaks[-1], color = retval$col)
    invisible(retval)
}
################################################################################

