

######################################
SRCDIR <- "/Users/xiaoxiaoliu/work/src/morphology_analysis/IVSCC"
source(paste0(SRCDIR, "/cv_functions.R"))

data_dir<-"/Users/xiaoxiaoliu/work/data/lims2/1027_pw_aligned/"

AssignedID_FN <- paste0(data_dir, "clustering_results/ward_all_ol_clipped/cluster_id.csv")
outdir<-paste0(data_dir,"clustering_results/cross_validation/random_shuffle")


feature<-read.csv(paste0(data_dir, "meta_merged_allFeatures.csv"), sep=",",header=TRUE,stringsAsFactors=FALSE)


# remove string columns
rnames <- feature$specimen_name
print(rnames)
#shuffle rnames
set.seed(001) # just to make it reproducible
rnames<-sample(rnames)



feature<-subset(feature, select = -c(specimen_name, specimen_id,cre_line,layer,dendrite_type,swc_file) )
print(colnames(feature))

FeatureMatrix<- data.matrix(feature[1:dim(feature)[2]])
rownames(FeatureMatrix)<-rnames


flag.plot<-TRUE
RUN_CVRF(AssignedID_FN, FeatureMatrix, outdir, flag.plot, "morph",100)


