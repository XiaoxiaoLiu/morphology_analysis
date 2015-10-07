

library(WGCNA)




plot_cm2 <- function(out1_id,out2_id,title_str,ystr,xstr) {

  idx <- match(out1_id$specimen_name, out2_id$specimen_name)

  tt <- table(out1_id$cluster_id,out2_id$cluster_id[idx])
  tt.rowsum <- apply(tt, 1, sum)
  tt.colsum <- apply(tt, 2, sum)
  rownames(tt) <- paste(rownames(tt),tt.rowsum,sep=":")
  colnames(tt) <- paste(colnames(tt),tt.colsum,sep=":")
  tt.norm <- tt / matrix(rep(apply(tt, 1, sum), ncol(tt)), ncol=ncol(tt))
  tt.norm.pct <- round(100 * tt.norm)
  labeledHeatmap(tt.norm.pct, textMatrix=tt, colors=rev(heat.colors(100)), xLabels=colnames(tt.norm), yLabels=rownames(tt.norm),
                 main=title_str, cex.lab=0.75, cex.text=0.75, xlab=xstr, ylab=ystr)
}



data_dir = "/Users/xiaoxiaoliu/work/data"
out1_id =  read.csv(paste(data_dir ,"lims2/0923_pw_aligned/clustering_results/ward_all_ol_clipped/cluster_id.csv",sep="/"))
out2_id =  read.csv(paste(data_dir ,"lims2/0923_pw_aligned/clustering_results/ap_all_ol_clipped/cluster_id.csv",sep="/"))

types_id = read.csv(paste(data_dir ,"lims2/0923_pw_aligned/staci_tag_id.csv",sep="/"))
types_id


plot_cm2(out1_id, out2_id,"confusion matrix","ward","affinity propagation")

plot_cm2(out1_id, types_id,"confusion matrix","ward","staci's types")


plot_cm2(out2_id, types_id,"confusion matrix","affinity propagation","staci's types")



