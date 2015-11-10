import os
from os import sys, path
import platform


if (platform.system() == "Linux"):
    WORK_PATH = "/local1/xiaoxiaol/work"
else:
    WORK_PATH = "/Users/xiaoxiaoliu/work"


p=  WORK_PATH + '/src/cell-type-analysis'
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
        swc_fn = df_db_tags['orca_path'][i].split('/')[-1]
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


def resample(original_dir,preprocessed_dir):

    if not os.path.exists(preprocessed_dir):
        os.mkdir(preprocessed_dir)

    for input_swc_path in glob.glob(original_dir + "/*.swc"):
        print input_swc_path
        swc_fn = input_swc_path.split('/')[-1]

        preprocessed_swc_fn = preprocessed_dir+'/' + swc_fn
        bn.resample(input_swc_path, preprocessed_swc_fn)  ## due to the pw alignment, no  alignment are necessary
        ##bn.pre_processing(input_swc_path, preprocessed_swc_fn)

    return


def cal_bn_features(preprocessed_dir):
    preprocessed_ANO = preprocessed_dir + "/preprocessed.ano"
    bn.genLinkerFile(preprocessed_dir, preprocessed_ANO)

    ##batch computing  generate features
    feature_file = preprocessed_dir+'/features.nfb'
    bn.batch_compute(preprocessed_ANO, feature_file)

    ###  convert feature file into csv file
    nfb.generateALLFeatureCSV(feature_file, preprocessed_dir + '/features_with_tags.csv')
    return




def main():
    ###############################################################################
    #data_DIR = '/data/mat/xiaoxiaol/data/lims2/0903_filtered_ephys_qc'
    data_DIR = WORK_PATH +'/data/lims2/0923_pw_aligned'
    # /original stores the downloaded swc files
    original_dir = data_DIR + "/pw_aligned"
    preprocessed_dir = data_DIR + "/preprocessed"
    ##
    db_tags_csv_file = data_DIR + '/0923_filtered_ephys_qc.csv'  # 153 neurons

    ###############################################################################

    # run blast neuron  features
    cal_bn_features(original_dir,preprocessed_dir)

    # clean up db tags
    df_db_tags= cleanup_query_csv(db_tags_csv_file)

    # merge all info
    df_features = pd.read_csv(preprocessed_dir + '/features_with_tags.csv')

    swc_file_names2 = []
    for i in range(df_features.shape[0]):
        swc_fn = df_features.swc_file[i].split('/')[-1]
        swc_file_names2.append(swc_fn)
    df_features['swc_file_name'] = pd.Series(swc_file_names2)

    merged = pd.merge(df_db_tags, df_features, how='inner', on=['swc_file_name'])
    merged.drop(['orca_path', 'swc_file_name', 'region_info'], axis=1, inplace=True)


    # add three more tags
    # modify height width depth



    merged.to_csv(preprocessed_dir + '/features_with_db_tags.csv')



    # add three more tags
    # modify height width depth

    return


if __name__ == "__main__":
      main()