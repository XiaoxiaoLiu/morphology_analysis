library(gplots)
##### load data ###########
scoreMatrix<-read.csv("~/work/data/lims2/neuron_recon_2/allFeatures_scoreMatrix.csv")

scoreMat=data.matrix(scoreMatrix)



merged<-read.csv("~/work/data/lims2/neuron_recon_2/merged_allFeatures.csv")


unique.dendrite_type <- as.character(unique(merged$dendrite_type_tag))
dendrite_type <- as.character(merged$dendrite_type_tag)
groupId <- match(dendrite_type,unique.dendrite_type)
mycolor_dendrite_type <- topo.colors(3)[groupId]


unique.cre_line <- as.character(unique(merged$cre_line))
cre_line <- as.character(merged$cre_line)
groupId <- match(cre_line,unique.cre_line)
mycolor_cre_line <- rainbow(13)[groupId]


unique.cell_shape <- as.character(unique(merged$cell_shape))
cell_shape <- as.character(merged$cell_shape)
groupId <- match(cell_shape,unique.cell_shape)
mycolor_cell_shape <- rainbow(11)[groupId]



m=scoreMat[,2:ncol(scoreMat)]
rownames(m) <- paste(merged$cell_shape)
colnames(m) <- paste(merged$orca_path)
m[m>150]<- 150

hm =heatmap(m,
            symm=TRUE,
            col=redgreen(95),
            margins=c(15,15), 
            hclustfun=function(d) hclust(d, method='ward'),
            RowSideColors=mycolor_dendrite_type,
            ColSideColors=mycolor_cell_shape,
            scale="none",cexRow=0.3, cexCol=0.3
            )
title('Hierachical Clustering based on all 34 Morph Features')

legend("bottomright",      # location of the legend on the heatmap plot
       legend = unique.dendrite_type, # category labels
       col = topo.colors(3),  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)

legend("topright",      # location of the legend on the heatmap plot
       legend = unique.cell_shape, # category labels
       col = rainbow(11),  # color key
       lty= 1,             # line style
       lwd = 10            # line width
)

#legend("topright",      # location of the legend on the heatmap plot
#       legend = unique.cre_line, # category labels
#       col =rainbow(13),  # color key
#       lty= 1,             # line style
#       lwd = 10            # line width
#)

#dev.off()               # close the PNG device
