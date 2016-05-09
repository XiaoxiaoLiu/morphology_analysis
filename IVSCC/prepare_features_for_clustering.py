__author__ = 'xiaoxiaol'

import os
from os import sys, path
import platform


if (platform.system() == "Linux"):
    WORK_PATH = "/local1/xiaoxiaol/work"
else:
    WORK_PATH = "/Users/xiaoxiaoliu/work"


p=  WORK_PATH + '/src/morphology_analysis'
sys.path.append(p)


import utilities.morph_nfb_2_csv as nfb

import blast_neuron.blast_neuron_comp as bn

import glob
import pandas as pd





### parse the query csv and generate a dataframe with required db tags
def cleanup_query_csv(db_tags_csv_file):
    df_db_tags = pd.read_csv(db_tags_csv_file)
    swc_file_names1 = []
    cre_lines = []
    layers = []
    for i in range(df_db_tags.shape[0]):
        swc_fn = df_db_tags['filename'][i].split('/')[-1]
        swc_file_names1.append(swc_fn)

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
            layer = df_db_tags['region_info'][i].split(', ')[-1]
        layers.append(layer)

    df_db_tags['swc_file_name'] = pd.Series(swc_file_names1)  ### add swc_file_name tag for merging
    df_db_tags['cre_line'] = pd.Series(cre_lines)
    df_db_tags['layer'] = pd.Series(layers)
    return df_db_tags

    # clean up db tagsdf_db_tags= cleanup_query_csv(db_tags_csv_file)


def main():

    data_DIR ="/data/mat/xiaoxiaol/data/lims2/ivscc_0411"
    feature_file = data_DIR +'/ivscc_0411.csv'
    df_features_with_tags = cleanup_query_csv(feature_file)
    df_features_with_tags.to_csv(data_DIR+'/ivscc_0411_features_with_meta.csv', index=False)
    dfg= df_features_with_tags.groupby('cre_line')



if __name__ == "__main__":
        main()