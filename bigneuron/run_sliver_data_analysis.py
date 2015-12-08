__author__ = 'xiaoxiaol'

import sys
import platform

if (platform.system() == "Linux"):
    WORK_PATH = "/local1/xiaoxiaol/work"
else:
    WORK_PATH = "/Users/xiaoxiaoliu/work"

p = WORK_PATH + '/src/morphology_analysis'
sys.path.append(p)

import recon_prescreening as rp
import plot_distances as plt_dist
import pandas as pd
import numpy as np



###### generate gold standard 166 table
gold_dir = "/data/mat/xiaoxiaol/data/gold166/checked_final_swcs"
GOLD_CSV = "/data/mat/xiaoxiaol/data/gold166/gold.csv"
gold_feature_csv= "/data/mat/xiaoxiaol/data/gold166/features_with_tags.csv"

#print "resample and sort gold166s"
#rp.resample_and_sort(gold_dir+"/original",gold_dir+"/resampled",gold_dir+"/sorted")
#rp.gen_gold_feature_csv(gold_dir+"/sorted", GOLD_CSV,gold_feature_csv)



###########################################################
data_DIR ="/data/mat/xiaoxiaol/data/reconstructions_2015_1207"
original_dir = data_DIR +"/auto_recons"
resampled_dir = data_DIR+ "/resampled"
sorted_dir = data_DIR +"/sorted"
lookup_image_id_table_file = data_DIR +"/image_name_lookup_table.csv"


######  resample and sort
rp.resample_and_sort(original_dir,resampled_dir,sorted_dir)



######  sliver data table
SILVER_CSV = data_DIR+'/recon_table.csv'

rp.recon_table_gen(sorted_dir,lookup_image_id_table_file,SILVER_CSV)


##### merge to get the common set between gold and silver
merged_csv_file = data_DIR+'/shared.csv'
rp.merge_gold_silver(GOLD_CSV,SILVER_CSV,merged_csv_file)


###### compute neuron distances
neuron_distance_csv=data_DIR+'/nd.csv'
rp.cal_neuron_dist(merged_csv_file,neuron_distance_csv )

## sort by sample size
df_nd = pd.read_csv(neuron_distance_csv)
algorithms = np.unique(df_nd.algorithm)

dfg = df_nd.groupby('algorithm')
sample_size_per_algorithm=np.zeros(algorithms.size)
for i in range( algorithms.size):
    sample_size_per_algorithm[i] = (dfg.get_group(algorithms[i]).shape[0])

order = sample_size_per_algorithm.argsort()
algorithms_ordered = algorithms[order[::-1]]
plt_dist.plot_sample_size(neuron_distance_csv,data_DIR,algorithms_ordered)



plt_dist.plot_neuron_distance(neuron_distance_csv, data_DIR,algorithms_ordered,CASE_BY_CASE_PLOT = 0)

###### compute blastneuron features
tmp_feature_csv = sorted_dir +'/tmp_features_with_tags.csv'

#rp.cal_bn_features(sorted_dir,tmp_feature_csv)
output_feature_csv = data_DIR +'/features_with_tags.csv'
rp.map_image_name(tmp_feature_csv,lookup_image_id_table_file, output_feature_csv)

bn_csv = data_DIR+"/blastneuron_ssd.csv"

rp.cal_blastneuron_distance(output_feature_csv,gold_feature_csv,merged_csv_file,bn_csv)

# plot
plt_dist.plot_blasneuron_distance(bn_csv,data_DIR,algorithms_ordered,CASE_BY_CASE_PLOT=0)



