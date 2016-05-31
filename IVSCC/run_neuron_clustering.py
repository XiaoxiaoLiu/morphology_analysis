import os
from os import sys, path
import platform
import numpy as np
import pandas as pd
from datetime import datetime



machine_name = platform.node()
if "ibs"  in machine_name:
    WORK_PATH = "/local1/xiaoxiaol/work"
if "pstar" in machine_name:
    WORK_PATH = "/home/xiaoxiaol"
if "swk"  in machine_name:
    WORK_PATH = "/Users/xiaoxiaoliu/work"


p=  WORK_PATH + '/src/morphology_analysis'
sys.path.append(p)



import morph_clustering.morph_cluster as fc
import IVSCC_features as features


def filter_featureset(all_feature_file,output_csv_file):

    feature_names  = features.get_feature_names('spiny')

    df_clean = pd.read_csv(all_feature_file)
    col_names = np.append(['swc_file_name','dendrite_type','specimen_id','specimen_name','region_info','cre_line'],feature_names)
    df_clean = df_clean[col_names]

    df_clean.dropna(inplace=True) # drop rows contain NAs
    #df_clean=df_clean[df_clean['max_path_distance'] >0 ]
    #df_clean=df_clean[df_clean['max_branch_order'] > 1 ]
    #df_clean=df_clean[df_clean['total_length'] > 10 ]

    df_clean.to_csv(output_csv_file)

    return feature_names



def cluster_analysis(clean_feature_file,feature_names,output_dir, method='ward'):
    print datetime.now().strftime('starting:%Y-%m-%d %H:%M:%S')

    df_f = pd.read_csv(clean_feature_file)


   # all_feature_names=gl_feature_names
    print "There are %d neurons in this dataset" % df_f.shape[0]


    REMOVE_OUTLIERS = 1 #clipping the dataset
    if REMOVE_OUTLIERS > 0:
        postfix = "_ol_removed"
    else:
        postfix = "_ol_clipped_5_glonly"


    if method == "ap" or method == "all":
        fc.run_affinity_propagation(df_f, feature_names, output_dir, postfix, "/data/mat/xiaoxiaol/data/lims2/ivscc_0519/SWC",REMOVE_OUTLIERS)


    num_clusters = 12
    if method == "ward" or method == "all":
        fc.run_ward_cluster(df_features=df_f, feature_names=feature_names, num_clusters=num_clusters,
                            output_dir= output_dir,
                            output_postfix=postfix,experiment_type='ivscc', low=3, high = 20, plot_heatmap=1,
                            RemoveOutliers=REMOVE_OUTLIERS, swc_path="/data/mat/xiaoxiaol/data/lims2/ivscc_0519/SWC")

    print datetime.now().strftime('end:%Y-%m-%d %H:%M:%S')
    return






def main():
     dataset ='IVSCC'


     if dataset=='IVSCC':
            data_DIR = "/data/mat/xiaoxiaol/data/lims2/ivscc_0519"
            output_dir = data_DIR+'/clustering_result'
            if not os.path.exists(output_dir):
                 os.system('mkdir '+output_dir)
            feature_tag_csv = data_DIR + '/spiny_features_curated.csv'
            output_clean_features_csv = data_DIR + '/spiny_features_filtered.csv'
            feature_names = filter_featureset (feature_tag_csv,output_clean_features_csv)

     #generate linkage
     cluster_analysis(output_clean_features_csv,feature_names,output_dir, method='ward')

     #merge cluster id


if __name__ == "__main__":
    main()
