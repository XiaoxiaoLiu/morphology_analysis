import os
from os import sys, path
import platform
import numpy as np
import pandas as pd

if (platform.system() == "Linux"):
    WORK_PATH = "/local1/xiaoxiaol/work"
else:
    WORK_PATH = "/Users/xiaoxiaoliu/work"


p=  WORK_PATH + '/src/morphology_analysis'
sys.path.append(p)

import clustering.feature_clustering as fc
import utilities.morph_nfb_2_csv as nfb




def cluster_analysis(all_feature_file,output_dir,feature_type='all', method='ward'):

    merged = pd.read_csv(all_feature_file)

    gl_feature_names = nfb.get_GL_feature_names()
    print gl_feature_names

    gmi_feature_names = nfb.get_GMI_feature_names()

    all_feature_names = np.append(gl_feature_names, gmi_feature_names)

    print "There are %d neurons in this dataset" % merged.shape[0]

    feature_names = all_feature_names
    if feature_type == "all":
        feature_names = all_feature_names
    if feature_type == "gmi":
        feature_names = gmi_feature_names


    REMOVE_OUTLIERS = 0

    postfix = "_" + feature_type


    if REMOVE_OUTLIERS > 0:
        postfix += "_ol_removed"
    else:
        postfix += "_ol_clipped"

    print "??"
    redundancy_removed_features_names = fc.remove_correlated_features(merged, feature_names, 0.95)
    print(" The %d features that are not closely correlated are %s" % (
        len(redundancy_removed_features_names), redundancy_removed_features_names))



    swc_screenshot_folder = None
    if method == "ap" or method == "all":
        num_clusters, dunn_index1 = fc.affinity_propagation(merged, redundancy_removed_features_names,
                                                     output_dir + '/ap' + postfix,
                                                    swc_screenshot_folder, REMOVE_OUTLIERS)

    num_clusters=13
    if method == "ward" or method == "all":
        dunn_index2 = fc.ward_cluster(merged, redundancy_removed_features_names, num_clusters,
                               output_dir + '/ward' + postfix, swc_screenshot_folder,
                               REMOVE_OUTLIERS)




def main():

    data_DIR = "/data/mat/xiaoxiaol/data/big_neuron/consensus_all/janelia_set1"
    output_dir = data_DIR+'/clustering_result'
    if not os.path.exists(output_dir):
         os.system('mkdir '+output_dir)
    feature_nbf="/data/mat/xiaoxiaol/data/big_neuron/consensus_all/janelia_set1/j1_smooth_features.nfb"
    output_feature_csv = data_DIR + '/j1_smooth_features_with_tags.csv'
    nfb.generateALLFeatureCSV(feature_nbf, output_feature_csv)

    logfile = output_dir+'/clustering.log'

    cluster_analysis(output_feature_csv,output_dir,feature_type='all', method='ward')



if __name__ == "__main__":
    main()
