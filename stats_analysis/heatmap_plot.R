features<-read.csv("~/work//data/lims2/neuron_recon_2/zscore_glFeatures.csv")
f<-features[,1:21]
d<-data.matrix(f)

library(gplots)
help(heatmap.2)
heatmap.2(t(d), col=redgreen(75), scale="column", key=T, keysize=1.5,density.info="none", trace="none",cexCol=0.9)
