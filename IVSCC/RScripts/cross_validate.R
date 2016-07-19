

######################################

# copy Cre_line_typeC.09102015 to data folder


SRCDIR <- "/Users/xiaoxiaoliu/work/src/morphology_analysis/IVSCC"
source(paste0(SRCDIR, "/cv_functions.R"))

data_dir<-"/Users/xiaoxiaoliu/work/data/lims2/pw_aligned_1223/"

AssignedID_FN <- paste0(data_dir, "0108/clustering_results/no_GMI/ward_aspiny_all_ol_clipped/cluster_id.csv")

outdir<-paste0(data_dir,"0108/cross_validation/no_GMI/ward_aspiny_all_ol_clipped")


feature<-read.csv(paste0(data_dir, "0108/aspiny_selected_features.csv"), sep=",",header=TRUE,stringsAsFactors=FALSE)
# remove string columns
rnames <- feature$specimen_name
feature<-subset(feature, select = -c(specimen_name, specimen_id,cre_line,region_info,dendrite_type,filename,swc_file_name) )

print(colnames(feature))

FeatureMatrix<- data.matrix(feature[1:dim(feature)[2]])
rownames(FeatureMatrix)<-rnames


flag.plot<-TRUE
RUN_CVRF(AssignedID_FN, FeatureMatrix, outdir, flag.plot, "morph",100)


