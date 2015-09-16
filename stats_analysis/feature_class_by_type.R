
library(gplots)

#merged<-read.csv("~/work/data/lims2/neuron_recon_2/merged.csv")
merged<-read.csv("~/work/data/lims2/0903_filtered_ephys_qc/preprocessed/features_with_db_tags.csv")

summary(merged)
names(merged)

dims=dim(merged)


################ remove outliers, normalize data
num_features = 34
y=data.matrix(merged[,8:(8+34-1)]) #-1:remove avg_r


m=apply(y,2,mean)
std=apply(y,2,sd)
mm=t(replicate(dims[1],m))
stdstd=t(replicate(dims[1],std))
zscores= (y-mm)/stdstd
zscores[zscores<-3] <- -3
zscores[zscores>3] <- 3



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

#pdf("~/work/data/lims2/neuron_recon_2/heatmap2_dendrite_type.pdf", height=10, width=10)
heatmap.2(t(newd),
          hclustfun=function(d) hclust(d, method='ward'),
          Colv=TRUE,
          col=redgreen(75), 
          scale="column", 
          #key=T,
          ColSideColors=mycolor1, 
          #keysize=1.5,
          #density.info="none", 
          #trace="none",cexCol=0.9,margins=c(12,15)
          )
#
legend("bottomleft",      # location of the legend on the heatmap plot
       legend = c("aspiny", "spiny", "sparse spiny"), # category labels
       col = mycolors,  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)
#dev.off