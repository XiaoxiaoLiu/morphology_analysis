library(gplots)
##### load data ###########
scoreMatrix<-read.csv("~/work/data/lims2/neuron_recon_2/gmi_scoreMatrix.csv")

scoreMat=data.matrix(scoreMatrix)



merged<-read.csv("~/work/data/lims2/neuron_recon_2/merged_gmiFeatures.csv")


m=scoreMat[,2:ncol(scoreMat)]
rownames(m) <- paste(merged$dendrite_type_tag)
colnames(m) <- paste(merged$cre_line)
m[m>150]<- 150





# creates a 5 x 5 inch image
#png("~/work/data/lims2/neuron_recon_2/heatmap_distance_matrix.png",    # create PNG for the heat map        
#    width = 5*300,        # 5 x 300 pixels
#    height = 5*300,
#    res = 300,            # 300 pixels per inch
#    pointsize = 8)        # smaller font size

hm =heatmap(m,
            symm=TRUE,
            col=redgreen(95),
            margins=c(12,15), 
            method='ward'
            )

title('GMI Feature Distance Matrix')
#dev.off()               # close the PNG device
