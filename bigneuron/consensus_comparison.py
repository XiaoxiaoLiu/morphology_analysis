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

#
# data_DIR ="/data/mat/xiaoxiaol/data/big_neuron/silver/20160113_merged_gold_gt"
# original_dir = data_DIR +"/auto_recons"
# lookup_image_id_table_file = data_DIR +"/../image_name_lookup_table.csv"
# time_csv = data_DIR + "/auto_recons/running_time_merged.csv"
# neuron_distance_csv = data_DIR +'/neuron_distances_with_gold.csv'



# compare valid consensus results with each individual algorithms

# input_log_path_csv=
# output_csv= data_DIR + "consensus.dist.csv"
# lookuptable=
# collect_consensus_distance(input_log_path_csv,output_csv, lookuptable )





#read nd distance csv


# merge by image file name
#
#
#
# # plot
# # ## sort by sample size
# df_nd = pd.read_csv(neuron_distance_csv)
# algorithms = np.unique(df_nd.algorithm)
#
# dfg = df_nd.groupby('algorithm')
# sample_size_per_algorithm = np.zeros(algorithms.size)
# for i in range( algorithms.size):
#     sample_size_per_algorithm[i] = (dfg.get_group(algorithms[i]).shape[0])
#
# order = sample_size_per_algorithm.argsort()
# algorithms_ordered = algorithms[order[::-1]]
#
#
# plt_dist.plot_similarities(neuron_distance_csv, data_DIR,algorithms_ordered,metric='neuron_distance',CASE_BY_CASE_PLOT = 0,value_label='Similarity (0~1) based on Average Neuron Distance (D1)')
