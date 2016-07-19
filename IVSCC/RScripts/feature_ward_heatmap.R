
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}


#Load necessary packages
library("gplots")
library("devtools")

#Load latest version of heatmap.3 function
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

data_dir<-"/data/mat/xiaoxiaol/data/lims2/0923_pw_aligned"
merged<-read.csv(paste(data_dir, "meta_merged_allFeatures.csv", sep="/"))
zscores<-read.csv(paste(data_dir, "zscore.csv", sep="/"))


dims=dim(merged)



# ZSCORE_THRESHOLD<-3.5
# zscores[zscores<-ZSCORE_THRESHOLD] <- -ZSCORE_THRESHOLD
# zscores[zscores>ZSCORE_THRESHOLD] <- ZSCORE_THRESHOLD


################ remove outliers, normalize data
num_features = dim(zscores)[2]
zscore_mat=data.matrix(zscores) #-1:remove avg_r
rownames(zscore_mat) = merged$specimen_name




unique.dt <- as.character(unique(merged$dendrite_type))
type <- as.character(merged$dendrite_type)
groupId <- match(type,unique.dt)
cl_dt= c("blue", "red", "gray")
mycolor_dt <- cl_dt[groupId]

unique.types <- as.character(unique(merged$types))
types <- as.character(merged$types)
groupId <- match(types,unique.types)
cl_types <- brewer.pal(length(unique.types),"Spectral")
mycolor_types <-cl_types[groupId]

unique.cre <- as.character(unique(merged$cre_line))
cre <- as.character(merged$cre_line)
groupId <- match(cre,unique.cre)
cl_cre <-brewer.pal(length(unique.cre),"Paired")
mycolor_cre <-cl_cre[groupId]

unique.layer <- as.character(unique(merged$layer_corrected))
layer <- as.character(merged$layer_corrected)
groupId <- match(layer,unique.layer)
cl_layer <- brewer.pal(length(unique.layer),"Blues")
mycolor_layer <-cl_layer[groupId]

clab=cbind(mycolor_types,mycolor_layer,mycolor_cre,mycolor_dt)
colnames(clab)=c("expert type", "layer","cre_line","dendrite type")

data = data.frame(y=zscore_mat, group = groupId)
aov.out = aov(zscore_mat ~ groupId, data)

pvalue = vector(,num_features)
for( i in 1:num_features){
  pvalue[i] = summary(aov.out)[[i]][["Pr(>F)"]][1]
}

rankP = rank(pvalue)
orderP = order(pvalue)


#yy=scale(zscores, center=TRUE, scale=TRUE)
newd=zscore_mat[,orderP]
#






#png(paste(data_dir, "clustering_results/R_heatmap2_dendrite_type.png",sep="/"),  width=800, height=750)
heatmap.3(t(newd),Rowv=FALSE,Colv=TRUE,symm = FALSE,distfun = function(x) dist(x,method = 'euclidean'),hclustfun=function(d) hclust(d,method='ward'), ColSideColors=clab, ColSideColorsSize=4, keysize=1.5,col=redgreen(75), density.info="none", trace="none", cexCol=0.9)

legend("bottomleft",
       legend=c(unique.dt," ", unique.cre, " ", unique.layer, " ",unique.types),
       fill=c(cl_dt,"white",cl_cre,"white",cl_layer,"white",cl_types),
       border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)

#dev.offo
#graphics.off()




#################  cell_shape

#
#

# legend("bottomleft",      # location of the legend on the heatmap plot
#        legend = unique.shape, # category labels
#        col = rainbow(11),  # color key
#        lty= 1,             # line style
#        lwd = 10            # line width
# )

