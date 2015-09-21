import os
from os import sys, path

p = path.dirname(path.dirname(path.abspath(__file__)))
p=  '/local1/xiaoxiaol/work/src/cell-type-analysis'
sys.path.append(p)

sys.path.append(p + '/utilities')
import morph_nfb_2_csv as nfb

sys.path.append(p + '/blast_neuron')
import blast_neuron_comp as bn

import glob
import pandas as pd


data_DIR = '/data/mat/xiaoxiaol/data/lims2/0903_filtered_ephys_qc'

# /original stores the downloaded swc files
original_dir = data_DIR + "/original"
preprocessed_dir = data_DIR + "/preprocessed"
if not os.path.exists(preprocessed_dir):
    os.mkdir(preprocessed_dir)

for input_swc_path in glob.glob(original_dir + "/*.swc"):
    print input_swc_path
    swc_fn = input_swc_path.split('/')[-1]

    preprocessed_swc_fn = data_DIR + '/preprocessed/' + swc_fn
    bn.pre_processing(input_swc_path, preprocessed_swc_fn)

preprocessed_ANO = preprocessed_dir + "/preprocessed.ano"
bn.genLinkerFile(preprocessed_dir, preprocessed_ANO)

##batch computing  generate features
feature_file = data_DIR + '/preprocessed/features.nfb'
bn.batch_compute(preprocessed_ANO, feature_file)



###  convert feature file into csv file
FEATURE_FILE = preprocessed_dir + '/features.nfb'
nfb.generateALLFeatureCSV(FEATURE_FILE, preprocessed_dir + '/features_with_tags.csv')


### merge with more database tags:
db_tags_csv_file = data_DIR + '/0903_filtered_ephys_qc.csv'
df_db_tags = pd.read_csv(db_tags_csv_file)

swc_file_names1 = []
cre_lines = []
layers = []
for i in range(df_db_tags.shape[0]):
    swc_fn = df_db_tags['orca_path'][i].split('/')[-1]
    swc_file_names1.append(swc_fn)

    if not pd.isnull(df_db_tags['dendrite_type'][i]):
        df_db_tags.set_value(i, 'dendrite_type', df_db_tags.dendrite_type[i].split(' - ')[-1])

    creline = None
    if not pd.isnull(df_db_tags['specimen_name'][i]):
        creline = df_db_tags['specimen_name'][i].split(';')[0]
    cre_lines.append(creline)

    layer = None
    if not pd.isnull(df_db_tags['region_info'][i]):
        layer = df_db_tags['region_info'][i].split(', ')[-1]
    layers.append(layer)

df_db_tags['swc_file_name'] = pd.Series(swc_file_names1)
df_db_tags['cre_line'] = pd.Series(cre_lines)
df_db_tags['layer'] = pd.Series(layers)

# merge all info
df_features = pd.read_csv(preprocessed_dir + '/features_with_tags.csv')

swc_file_names2 = []
for i in range(df_features.shape[0]):
    swc_fn = df_features.swc_file[i].split('/')[-1]
    swc_file_names2.append(swc_fn)
df_features['swc_file_name'] = pd.Series(swc_file_names2)

merged = pd.merge(df_db_tags, df_features, how='inner', on=['swc_file_name'])
merged.drop(['orca_path', 'swc_file_name', 'region_info'], axis=1, inplace=True)
merged.to_csv(preprocessed_dir + '/features_with_db_tags.csv')

