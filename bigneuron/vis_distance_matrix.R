data_dir<-"/data/mat/xiaoxiaol/data/big_neuron/consensus_all/janelia_set1/clustering_result"


scoreMatrix <-read.csv("/data/mat/xiaoxiaol/data/big_neuron/consensus_all/janelia_set1/clustering_result/ward_all_ol_clipped/zscore.csv")
scoreMat=data.matrix(scoreMatrix)
d=scoreMat[,1:ncol(scoreMat)]
dim(d)


#resample only 500 for testing
subset=seq(1, dim(d)[1], dim(d)[1]/500)
d=d[subset,]
dim(d)

library(mclust)
# Run the function to see how many clusters
# it finds to be optimal, set it to search for
# at least 1 model and up 20.
d_clust <- Mclust(d, G=8:11)
m.best <- dim(d_clust$z)[2]
cat("model-based optimal number of clusters:", m.best, "\n")
# 4 clusters
plot(d_clust)





library(apcluster)
d.apclus <- apcluster(negDistMat(r=2), d)
cat("affinity propogation optimal number of clusters:", length(d.apclus@clusters), "\n")
# 4
heatmap(d.apclus)
plot(d.apclus, d)