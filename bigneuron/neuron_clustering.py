import os
from os import sys, path
import platform
import numpy as np
import pandas as pd
from datetime import datetime

if (platform.system() == "Linux"):
    WORK_PATH = "/local1/xiaoxiaol/work"
else:
    WORK_PATH = "/Users/xiaoxiaoliu/work"


p=  WORK_PATH + '/src/morphology_analysis'
sys.path.append(p)

import morph_clustering.morph_cluster as fc
import utilities.morph_nfb_2_csv as nfb


def filter_featureset(all_feature_file,output_csv_file):


    gl_feature_names = nfb.get_GL_feature_names('no_radii')
    gmi_feature_names = nfb.get_GMI_feature_names()
    feature_names = np.append(gl_feature_names, gmi_feature_names)
    print feature_names
    df_clean = pd.read_csv(all_feature_file)
    col_names = np.append(['swc_file_name'],feature_names)
    df_clean = df_clean[col_names]

    df_clean.dropna(inplace=True) # drop rows contain NAs
    df_clean=df_clean[df_clean['max_path_distance'] >0 ]
    df_clean=df_clean[df_clean['max_branch_order'] > 1 ]
    df_clean=df_clean[df_clean['total_length'] > 10 ]


    df_clean.to_csv(output_csv_file)


    return feature_names



def cluster_analysis(clean_feature_file,feature_names,output_dir, method='ward'):
    print datetime.now().strftime('starting:%Y-%m-%d %H:%M:%S')

    df_f = pd.read_csv(clean_feature_file)


   # all_feature_names=gl_feature_names
    print "There are %d neurons in this dataset" % df_f.shape[0]


    REMOVE_OUTLIERS = 0 #clipping the dataset
    if REMOVE_OUTLIERS > 0:
        postfix = "_ol_removed"
    else:
        postfix = "_ol_clipped"


    if method == "ap" or method == "all":
        fc.run_affinity_propagation(df_f, feature_names, output_dir, postfix)


    num_clusters=29
    if method == "ward" or method == "all":
        fc.run_ward_cluster(df_features=df_f, feature_names=feature_names, num_clusters=num_clusters,output_dir= output_dir,
                                          output_postfix=postfix,experiment_type='bigneuron')

    print datetime.now().strftime('end:%Y-%m-%d %H:%M:%S')
    return



def main():
     dataset ='taiwan'


     if dataset=='Janelia':
            data_DIR = "/data/mat/xiaoxiaol/data/big_neuron/consensus_all/janelia_set1"
            output_dir = data_DIR+'/clustering_result'
            if not os.path.exists(output_dir):
                 os.system('mkdir '+output_dir)
            feature_tag_csv = data_DIR + '/j1_smooth_features_with_tags.csv'
            output_clean_features_csv = data_DIR + '/j1_smooth_features_with_tags_cleaned.csv'

     if dataset=='taiwan':
            data_DIR = "/data/mat/xiaoxiaol/data/big_neuron/consensus_all/taiwan"
            output_dir = data_DIR+'/clustering_result'
            if not os.path.exists(output_dir):
                 os.system('mkdir '+output_dir)
            feature_tag_csv = data_DIR + '/taiwan_smooth_features_with_tags.csv'
            output_clean_features_csv = data_DIR + '/taiwan_smooth_features_with_tags_cleaned.csv'

     feature_names=filter_featureset(feature_tag_csv,output_clean_features_csv)
     cluster_analysis(output_clean_features_csv,feature_names,output_dir, method='ward')



if __name__ == "__main__":
    main()
