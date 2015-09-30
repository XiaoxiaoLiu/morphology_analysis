library(gplots)


data_dir<-"/data/mat/xiaoxiaol/data/lims2/0923_pw_aligned"
merged<-read.csv(paste(data_dir, "preprocessed/features_with_db_tags.csv", sep="/"))


dims=dim(merged)



ZSCORE_THRESHOLD<-5


################ remove outliers, normalize data
num_features = 34
y=data.matrix(merged[,8:(8+34-1)]) #-1:remove avg_r


m=apply(y,2,mean)
std=apply(y,2,sd)
mm=t(replicate(dims[1],m))
stdstd=t(replicate(dims[1],std))
zscores= (y-mm)/stdstd
zscores[zscores<-ZSCORE_THRESHOLD] <- -ZSCORE_THRESHOLD
zscores[zscores>ZSCORE_THRESHOLD] <- ZSCORE_THRESHOLD


unique.type <- as.character(unique(merged$dendrite_type))
type <- as.character(merged$dendrite_type)

groupId <- match(type,unique.type)
mycolors= c("gray", "blue", "red")
mycolor1 <- mycolors[groupId]


data = data.frame(y=zscores, group = groupId)
aov.out = aov(zscores ~ groupId, data)

pvalue = vector(,num_features)
for( i in 1:num_features){
  pvalue[i] = summary(aov.out)[[i]][["Pr(>F)"]][1]
}

rankP = rank(pvalue)
orderP = order(pvalue)


#yy=scale(zscores, center=TRUE, scale=TRUE)
newd=zscores[,orderP]

png(paste(data_dir, "clustering_results/R_heatmap2_dendrite_type.png",sep="/"),  width=800, height=750)
heatmap.2(t(newd),Rowv=TRUE,Colv=TRUE,symm = FALSE,hclustfun=function(d) hclust(d,method='ward'), ColSideColors=mycolor1, keysize=1.5,col=redgreen(75), density.info="none", trace="none", cexCol=0.9, margins=c(10,15))

legend("topright",      # location of the legend on the heatmap plot
       legend = c("aspiny", "spiny", "sparse spiny"), # category labels
       col = mycolors,  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)
#dev.offo
graphics.off()




#################  cell_shape
# unique.shape <- as.character(unique(merged$cell_shape))
# shape <- as.character(merged$cell_shape)
# groupId <- match(shape,unique.shape)
# mycolor1 <- rainbow(11)[groupId]
#
#
# data = data.frame(y = zscores, group = groupId)
# aov.out = aov(zscores ~ groupId, data)
#
# pvalue = vector(,21);
# for( i in 1:21){
#   pvalue[i] = summary(aov.out)[[i]][["Pr(>F)"]][1]
# }
#
# rankP = rank(pvalue)
# orderP = order(pvalue)
#
#
# yy=scale(zscores, center=TRUE, scale=TRUE)
# newd=yy[,orderP]
#
#
# heatmap.2(t(newd), Rowv=FALSE, Colv=TRUE, col=redgreen(75), scale="column", key=T, ColSideColors=mycolor1, keysize=1.5,density.info="none", trace="none",cexCol=0.9,margins=c(12,15))
# legend("bottomleft",      # location of the legend on the heatmap plot
#        legend = unique.shape, # category labels
#        col = rainbow(11),  # color key
#        lty= 1,             # line style
#        lwd = 10            # line width
# )
############# spiny nonspinydim(y)
