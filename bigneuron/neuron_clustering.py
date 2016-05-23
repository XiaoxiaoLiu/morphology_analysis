import os
from os import sys, path
import platform
import numpy as np
import pandas as pd
from datetime import datetime



machine_name = platform.node()
if machine_name.contains( "ibs"):
    WORK_PATH = "/local1/xiaoxiaol/work"
if machine_name.contains( "pstar"):
    WORK_PATH = "/home/xiaoxiaol"
if machine_name.contains( "swk"):
    WORK_PATH = "/Users/xiaoxiaoliu/work"


p=  WORK_PATH + '/src/morphology_analysis'
sys.path.append(p)



import morph_clustering.morph_cluster as fc
import utilities.morph_nfb_2_csv as nfb


def filter_featureset(all_feature_file,output_csv_file):


    gl_feature_names = nfb.get_GL_feature_names('no_radii')
    gmi_feature_names = nfb.get_GMI_feature_names()

    feature_names = np.append(gl_feature_names, gmi_feature_names)

    #feature_names = gl_feature_names

    feature_names = feature_names[feature_names !='num_stems']
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


    REMOVE_OUTLIERS = 1 #clipping the dataset
    if REMOVE_OUTLIERS > 0:
        postfix = "_ol_removed"
    else:
        postfix = "_ol_clipped_5_glonly"


    if method == "ap" or method == "all":
        fc.run_affinity_propagation(df_f, feature_names, output_dir, postfix)


    num_clusters = 1000
    if method == "ward" or method == "all":
        fc.run_ward_cluster(df_features=df_f, feature_names=feature_names, num_clusters=num_clusters,output_dir= output_dir,
                                          output_postfix=postfix,experiment_type='bigneuron', low=500, high = 1500, plot_heatmap=0, RemoveOutliers=REMOVE_OUTLIERS)

    print datetime.now().strftime('end:%Y-%m-%d %H:%M:%S')
    return



def add_nt_meta_taiwan(output_clean_features_csv,output_features_csv):
     df_f=pd.read_csv(output_clean_features_csv)

     old_names = ['5HT1bMARCM', 'ChaMARCM', 'DVGLUTMARCM', 'DVGlutMARCM', 'DVGlutMARCMFemale',
                 'DvGlutMARCM', 'FruMARCM', 'G0239', 'GH146MARCM', 'GadMARCM', 'HMARCM',
                 'THMARCM', 'TMARCM', 'TPHMARCH', 'TPHMARCM', 'VAMFlipOut', 'dTdc', 'dTdc2MARCM','npfMARCM']
     clean_names = ['5HT1b', 'Cha', 'DVGluT', 'DVGlut', 'DVGlut',
                 'DVGlut', 'Fru', 'G0239', 'GH146', 'Gad', np.nan,
                 'TH', np.nan, 'TPH', 'TPH', 'VAMFlipOut', 'dTdc', 'dTdc2','npf']
     nt_type_lut = dict(zip(old_names, clean_names))
     df_f['nt_type'] = df_f['swc_file_name'].apply(lambda x :pd.Series(nt_type_lut[x.split('/')[-1].split('.')[1].split('-')[0]]))

     print np.unique (df_f['nt_type'])

     df_f.dropna(inplace=True)
     df_f.to_csv(output_features_csv)


def main():
     dataset ='taiwan'

     if dataset=='Janelia':
            data_DIR = "/data/mat/xiaoxiaol/data/big_neuron/consensus_all/janelia_set1"
            output_dir = data_DIR+'/clustering_result'
            if not os.path.exists(output_dir):
                 os.system('mkdir '+output_dir)
            feature_tag_csv = data_DIR + '/j1_smooth_features_with_tags.csv'
            output_clean_features_csv = data_DIR + '/j1_smooth_features_with_tags_cleaned.csv'

            feature_names=filter_featureset(feature_tag_csv,output_clean_features_csv)

            output_features_csv = output_clean_features_csv

     if dataset=='taiwan':
            data_DIR = "/data/mat/xiaoxiaol/data/big_neuron/consensus_all/taiwan"
            output_dir = data_DIR+'/clustering_result'
            if not os.path.exists(output_dir):
                 os.system('mkdir '+output_dir)
            feature_tag_csv = data_DIR + '/taiwan_smooth_features_with_tags.csv'
            output_clean_features_csv = data_DIR + '/taiwan_smooth_features_with_tags_cleaned.csv'

            feature_names = filter_featureset (feature_tag_csv,output_clean_features_csv)
            output_features_csv = data_DIR + '/taiwan_smooth_gl_features_with_tags_cleaned_nttype.csv'
            add_nt_meta_taiwan(output_clean_features_csv,output_features_csv)

     #generate linkage
     cluster_analysis(output_features_csv,feature_names,output_dir, method='ward')

     #load linkage
     #fc.plot_linkage(output_feature_csv,linkage)



if __name__ == "__main__":
    main()
