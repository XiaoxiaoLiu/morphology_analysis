read.csv("~/work//data/lims2/manual/zscore_glFeatures_50.csv")
features<-read.csv("~/work//data/lims2/manual/zscore_glFeatures_50.csv")
f<-features[,2:22]
data.matrix(f)
d<-data.matrix(f)
install.packages("gplots")
library(gplots)
help(heatmap.2)
heatmap.2(d, col=redgreen(75), scale="row", key=T, keysize=1.5,
density.info="none", trace="none",cexCol=0.9)
savehistory("~/work/v3d/vaa3d_tools/blastneuron/scripts/heatmap_r.Rhistory")

