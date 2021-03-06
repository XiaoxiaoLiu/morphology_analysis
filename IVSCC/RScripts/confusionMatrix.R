library(WGCNA)

plot_cm2 <- function(out1_id,out2_id,title_str,ystr,xstr, output_image) {
  png(filename=output_image)
  a = out1_id$cluster_id
  names(a) = out1_id$specimen_name

  b = out2_id$cluster_id
  names(b) = out2_id$specimen_name


  #drop the specimen name that does not match?

  idx <- match(names(a), names(b))

  tt <- table(a,b[idx])
  tt.rowsum <- apply(tt, 1, sum)
  tt.colsum <- apply(tt, 2, sum)
  rownames(tt) <- paste(rownames(tt),tt.rowsum,sep=":")
  colnames(tt) <- paste(unique(out2_id$types),tt.colsum,sep=":")
  tt.norm <- tt / matrix(rep(apply(tt, 1, sum), ncol(tt)), ncol=ncol(tt))
  tt.norm.pct <- round(100 * tt.norm)
  # reorder confusion matrix
  d<-tt.norm.pct

  for (col in 1: dim(d)[2]  )
  {
    #find the max value and record the col index
    row_order<-1:dim(d)[1]
    max_row_id <- which.max(d[col:dim(d)[1],col] )

    max_row_id <- (max_row_id + col -1)
    row_order[col]<-max_row_id
    row_order[max_row_id]<- col

    d<- d[row_order,]
    tt<-tt[row_order,]
    tt.norm<-tt.norm[row_order,]
  }

  for (row in 1: dim(d)[1]  )
  {

     #find the max value and record the col index
     col_order<-1:dim(d)[2]
     max_col_id <- which.max(d[row,row:dim(d)[2]] )

     max_col_id <- (max_col_id + row -1)
     col_order[row]<-max_col_id
     col_order[max_col_id]<- row

     d<- d[,col_order]
     tt<-tt[,col_order]
     tt.norm<-tt.norm[,col_order]
  }






  mycolors =
  labeledHeatmap(d, textMatrix=tt, colors=rev(heat.colors(100)), xLabels=colnames(tt.norm), yLabels=rownames(tt.norm),
                 main=title_str, cex.lab=0.95, cex.text=0.75, xlab=xstr, ylab=ystr)

  dev.off()
}


data_dir = "/Users/xiaoxiaoliu/work/data/lims2/0923_pw_aligned"


types_id = read.csv(paste(data_dir ,"staci_tag_id.csv",sep="/"))
unique.type <- as.character(unique(types_id$types))
type <- as.character(types_id$types)
groupId <- match(type,unique.type)
types_id$cluster_id = groupId


out1_id =  read.csv(paste(data_dir ,"clustering_results/ward_all_ol_clipped/cluster_id.csv",sep="/"))
plot_cm2(out1_id, types_id,"confusion matrix","ward","", paste(data_dir, "confusion_ward_types.png",sep="/"))


# out2_id =  read.csv(paste(data_dir ,"clustering_results/ward_mrmr_ol_clipped/cluster_id.csv",sep="/"))
# plot_cm2(out2_id, types_id,"confusion matrix","ward mrmr","staci's types", paste(data_dir, "confusion_wardmrmr_types.png",sep="/"))
#plot_cm2(out2_id, types_id,"confusion matrix","ward mrmr","staci's types", paste(data_dir, "confusion_wardmrmr_types.png",sep="/"))


out3_id =  read.csv(paste(data_dir ,"clustering_results/ap_all_ol_clipped/cluster_id.csv",sep="/"))
plot_cm2(out3_id, types_id,"confusion matrix","affinity propagation","", paste(data_dir, "confusion_ap_types.png",sep="/") )


# out4_id =  read.csv(paste(data_dir ,"clustering_results/ap_mrmr_ol_clipped/cluster_id.csv",sep="/"))
# plot_cm2(out4_id, types_id,"confusion matrix","affinity propagation mrmr","staci's types", paste(data_dir, "confusion_apmrmr_types.png",sep="/"))
#


# out5_id =  read.csv(paste(data_dir ,"morph.interactions.label.csv",sep="/"))
# plot_cm2(out5_id, types_id,"confusion matrix","mean_shift_interactions","staci's types", paste(data_dir, "confusion_meanshiftinteract_types.png",sep="/"))


# out6_id =  read.csv(paste(data_dir ,"Final.Ephys.seed0.Regular.P0.01.GrpbyPC.FSbyNA.IncPC.PCA.0.058.csv",sep="/"))
# #plot_cm2(out6_id, types_id,"confusion matrix","iterative PCA","staci's types", paste(data_dir, "confusion_iterativepca_types.png",sep="/"))


































#plot_cm2(out1_id, out2_id,"confusion matrix","ward","affinity propagation")


