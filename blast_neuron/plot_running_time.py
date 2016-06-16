__author__ = 'xiaoxiaol'



import sys
import os
import platform

if (platform.system() == "Linux"):
    WORK_PATH = "/local1/xiaoxiaol/work"
else:
    WORK_PATH = "/Users/xiaoxiaoliu/work"

p = WORK_PATH + '/src/morphology_analysis'
sys.path.append(p)

import  bigneuron.recon_prescreening as rp
import  bigneuron.plot_distances as plt_dist
import pandas as pd
import numpy as np


data_DIR ="/data/mat/xiaoxiaol/data/reconstructions_2015_1207"
original_dir = data_DIR +"/auto_recons"
resampled_dir = data_DIR+ "/resampled"
sorted_dir = data_DIR +"/sorted"
lookup_image_id_table_file = data_DIR +"/image_name_lookup_table.csv"



########################  increamenatal data ##########################
updated_data_DIR = "/data/mat/xiaoxiaol/data/reconstructions_20151214"
new_data_Dir = "/data/mat/xiaoxiaol/data/reconstructions_2015_1214"

# copy old data to new folder ( to avoid errors while reorganizing)
#os.system( 'cp -r  '+data_DIR + "  "+ new_data_Dir+'/auto_recons')


#replace and add swc files and log files
#shell scripts update_data.sh

#running_time
#vim edit merge running_time.csv into one

data_DIR = new_data_Dir
original_dir = data_DIR +"/auto_recons"
resampled_dir = data_DIR+ "/resampled"
sorted_dir = data_DIR +"/sorted"




neuron_distance_csv = data_DIR+'/nd.csv'

df_nd = pd.read_csv(neuron_distance_csv)
algorithms = np.unique(df_nd.algorithm)
print algorithms

dfg = df_nd.groupby('algorithm')
sample_size_per_algorithm = np.zeros(algorithms.size)
for i in range( algorithms.size):
    print algorithms[i]
    sample_size_per_algorithm[i] = (dfg.get_group(algorithms[i]).shape[0])

order = sample_size_per_algorithm.argsort()
algorithms_ordered = algorithms[order[::-1]]


time_csv="/data/mat/xiaoxiaol/data/reconstructions_2015_1214/auto_recons/running_time.csv"
output_time_csv= "/data/mat/xiaoxiaol/data/reconstructions_2015_1214/running_time_algorithm.csv"
algorithm_plugin_match_csv ="/data/mat/xiaoxiaol/data/reconstructions_2015_1214/ported_neuron_tracing_spreadsheet.csv"
rp.summerize_running_time(time_csv, algorithm_plugin_match_csv,output_time_csv)


df_gold = pd.read_csv(GOLD_CSV)
df_silver = pd.read_csv(SILVER_CSV)
#print df_silver.columns
#print df_gold.columns

df_share = pd.merge(df_silver,df_gold,on="image_file_name")


plt_dist.plot_running_time(output_time_csv, data_DIR,algorithms_ordered)
