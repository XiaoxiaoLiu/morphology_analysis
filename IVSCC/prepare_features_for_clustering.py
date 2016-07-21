__author__ = 'xiaoxiaol'

import os
from os import sys, path
import platform
import numpy as np


if (platform.system() == "Linux"):
    WORK_PATH = "/local1/xiaoxiaol/work"
else:
    WORK_PATH = "/Users/xiaoxiaoliu/work"


p=  WORK_PATH + '/src/morphology_analysis'
sys.path.append(p)



import pandas as pd
import IVSCC_morph_features as features




### parse the query csv and generate a dataframe with required db tags
def cleanup_query_csv(db_tags_csv_file):
    df_db_tags = pd.read_csv(db_tags_csv_file)
    swc_file_names = []
    cre_lines = []
    layers = []
    for i in range(df_db_tags.shape[0]):
        swc_fn = df_db_tags['filename'][i].split('/')[-1]
        #_m.swc ==> _m_pia.swc
        swc_fn=swc_fn[:-4]+"_pia.swc"
        swc_file_names.append(swc_fn)

        if not pd.isnull(df_db_tags['dendrite_type'][i]):
            df_db_tags.set_value(i, 'dendrite_type', df_db_tags.dendrite_type[i].split(' - ')[-1])

        creline = 'NA'
        if not pd.isnull(df_db_tags['specimen_name'][i]):
            creline = df_db_tags['specimen_name'][i].split(';')[0]
        if creline == 'Sst-IRES-Cre-19768.06.02.01' : # database error
            creline = 'Sst-IRES-Cre'
        cre_lines.append(creline)

        layer = 'NA'
        if not pd.isnull(df_db_tags['region_info'][i]):
            layer = df_db_tags['region_info'][i].split(' ')[-1]
        layers.append(layer)

    df_db_tags['swc_file_name'] = pd.Series(swc_file_names)  ### add swc_file_name tag for merging
    df_db_tags['cre_line'] = pd.Series(cre_lines)
    df_db_tags['layer'] = pd.Series(layers)
    return df_db_tags

    # clean up db tagsdf_db_tags= cleanup_query_csv(db_tags_csv_file)


def main():

    data_DIR = '/data/mat/xiaoxiaol/data/lims2/ivscc_0607'
    feature_file = data_DIR +'/features_morph_2016.06.07.csv'

    df_features_with_tags = cleanup_query_csv(feature_file)

    feature_meta_file= data_DIR+'/ivscc_0607_new_aligned_with_meta.csv'
    df_features_with_tags.to_csv(feature_meta_file, index=False)



    out_dir=data_DIR
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)


    ########################################################
    all_feature_file = feature_meta_file
    #########################################################
    df_features = pd.read_csv(all_feature_file)


    apical_feature_names = features.get_feature_names('apical')



    meta_feature_names = np.array(['specimen_name','specimen_id','dendrite_type','cre_line','region_info','layer','filename','swc_file_name'])


    all_dendritic_feature_names = features.get_feature_names('all_dendrite')

    df_features[all_dendritic_feature_names]= df_features[all_dendritic_feature_names].astype(float)
    print "There are %d neurons in this dataset" % df_features.shape[0]
    print "Dendrite types: ", np.unique(df_features['dendrite_type'])


    # df_features_all = df_features[np.append(meta_feature_names,all_dendritic_feature_names)]
    # df_features_all.to_csv(data_DIR+'/0108/all_dendrite_features.csv')

    df_groups = df_features.groupby(['dendrite_type'])
    df_spiny = df_groups.get_group('spiny')
    print "There are %d neurons are spiny\n\n" % df_spiny.shape[0]

    # chose different sets of features
    df_w_spiny = df_spiny[np.append(meta_feature_names,apical_feature_names)]
    df_w_spiny.to_csv(data_DIR +'/spiny_apical_features.csv', index=False)
    #
    # df_w_spiny = df_spiny[np.append(meta_feature_names,apical_no_z_feature_names)]
    # print len(apical_no_z_feature_names)
    # df_w_spiny.to_csv(data_DIR +'/spiny_apical_no_z_features.csv', index=False)
    #
    #
    # df_w_spiny = df_spiny[np.append(meta_feature_names,spiny_feature_names)]
    # df_w_spiny.to_csv(data_DIR +'/spiny_apical_basal_feature.csv', index=False)
    #
    # df_w_spiny = df_spiny[np.append(meta_feature_names,spiny_no_z_feature_names)]
    # df_w_spiny.to_csv(data_DIR +'/spiny_apical_basal_no_z_feature.csv', index=False)



    df_aspiny =  pd.concat([df_groups.get_group('aspiny'),df_groups.get_group('sparsely spiny')],axis=0)
    print "There are %d neurons are aspiny " % df_aspiny.shape[0]

    aspiny_feature_names =  features.get_feature_names('aspiny')
    df_w_aspiny = df_aspiny[np.append(meta_feature_names,aspiny_feature_names)]
    df_w_aspiny.to_csv(data_DIR +'/aspiny_features.csv', index=False)







if __name__ == "__main__":
        main()
