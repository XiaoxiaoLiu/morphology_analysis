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


##################    gold standard dataset processing #########################################
###### only need to do it once
gold_dir = "/data/mat/xiaoxiaol/data/gold166/checked_final_swcs"
GOLD_CSV = "/data/mat/xiaoxiaol/data/gold166/gold.csv"  #manually generated
gold_feature_csv= "/data/mat/xiaoxiaol/data/gold166/features_with_tags.csv"

#print "resample and sort gold166s"
##rp.resample_and_sort(gold_dir+"/original",gold_dir+"/resampled",gold_dir+"/sorted")
#rp.gen_gold_feature_csv(gold_dir+"/original",output_gold_csv_file= GOLD_CSV,output_gold_feature_csv=gold_feature_csv)

algorithm_plugin_match_csv = "/data/mat/xiaoxiaol/data/big_neuron/silver/ported_neuron_tracing_spreadsheet.csv"

###########################  preprocessing and organize data ##################################
data_DIR ="/data/mat/xiaoxiaol/data/big_neuron/silver/20160113_merged_gold_gt"
original_dir = data_DIR +"/auto_recons"
#resampled_dir = data_DIR+ "/resampled"
#sorted_dir = data_DIR +"/sorted"
lookup_image_id_table_file = data_DIR +"/../image_name_lookup_table.csv"
time_csv = data_DIR + "/auto_recons/running_time_merged.csv"


neuron_distance_csv = data_DIR +'/neuron_distances_with_gold.csv'


ND = 1
BD = 0
COMPUTE = 1
MEDIAN = 0
################################################################################


if COMPUTE:
        ######  resample  #?and sort
        #rp.resample_and_sort(original_dir,resampled_dir,sorted_dir,GEN_QSUB=0,overwrite_sorted=0)


        ######  sliver data table
        SILVER_CSV = data_DIR+'/recon_table.csv'
        rp.recon_table_gen(original_dir,lookup_image_id_table_file,SILVER_CSV)

        #####  merge to get the common set between gold and silver
        merged_csv_file = data_DIR+'/recon_shared_with_gold_set.csv'
        rp.merge_gold_silver(GOLD_CSV,SILVER_CSV,merged_csv_file)


        #####  report which gold dataset did not have any recons?
        df_merge = pd.read_csv(merged_csv_file)
        df_gold = pd.read_csv(GOLD_CSV)
        m = pd.unique(df_merge.image_file_name)
        g = pd.unique(df_gold.image_file_name)

        print "\n\nGold dataset contains " +str(g.size) +" image dataset"
        print "There are " + str(df_merge.shape[0])+" reconstructions are generated from " + str(pd.unique(df_merge.algorithm).size) +" algorithms."
        for i in g:
           if i not in m:
              print "No reconstructions for image: " + i
        ###########################   distance calculation  ########################################


        #### compute neuron distances

        if ND :
            ##rp.cal_neuron_dist(merged_csv_file,neuron_distance_csv,overwrite_existing = 0, old_output_csv= data_DIR+'/nd_old.csv') # build on top of previous results
            rp.cal_neuron_dist(merged_csv_file,neuron_distance_csv,overwrite_existing = 0,GEN_QSUB = 1) # build on top of previous results


        ## post processing remove invalid
        df_nd = pd.read_csv(neuron_distance_csv)
        df_nd1 = df_nd[df_nd['neuron_distance'] != -1]  # empty node swc will have nd reported as "-1", remove those invalid recons
        df_nd2= df_nd1.dropna(axis=0)
        df_nd2.to_csv(neuron_distance_csv, index=False)




        ###### compute blastneuron features
        if BD > 0:
            tmp_feature_csv = original_dir +'/tmp_features_with_tags.csv'

            rp.cal_bn_features(original_dir,tmp_feature_csv)
            output_feature_csv = data_DIR +'/features_with_tags.csv'
            rp.map_image_name(tmp_feature_csv,lookup_image_id_table_file, output_feature_csv)

            df_m = pd.read_csv(output_feature_csv)
            df_m1 = df_m[df_m['num_nodes'] != 0]
            df_m1.to_csv(output_feature_csv, index=False)


            bn_csv = data_DIR+"/blastneuron_lm_ssd.csv"
            rp.cal_blastneuron_distance(output_feature_csv,gold_feature_csv,merged_csv_file,output_csv=bn_csv, LMEASURE_ONLY = 1)



###########################  plotting ###########################################################

# ## sort by sample size
df_nd = pd.read_csv(neuron_distance_csv)
algorithms = np.unique(df_nd.algorithm)

dfg = df_nd.groupby('algorithm')
sample_size_per_algorithm = np.zeros(algorithms.size)
for i in range( algorithms.size):
    sample_size_per_algorithm[i] = (dfg.get_group(algorithms[i]).shape[0])

order = sample_size_per_algorithm.argsort()
algorithms_ordered = algorithms[order[::-1]]



# plot
if BD >0:
   plt_dist.plot_blasneuron_distance(bn_csv,data_DIR,algorithms_ordered,CASE_BY_CASE_PLOT=1)
   plt_dist.plot_similarities(bn_csv, data_DIR,algorithms_ordered,metric='SSD',CASE_BY_CASE_PLOT = 0,value_label='Similarity (0~1) on Global Morph Feature Score')



plt_dist.plot_similarities(neuron_distance_csv, data_DIR,algorithms_ordered,metric='neuron_distance',CASE_BY_CASE_PLOT = 0,value_label='Similarity (0~1) based on Average Neuron Distance (D1)')
plt_dist.plot_similarities(neuron_distance_csv, data_DIR,algorithms_ordered,metric='neuron_difference',CASE_BY_CASE_PLOT = 0, value_label='Similarity (0~1) based on Neuron Difference Score (D2*D3)')
plt_dist.plot_similarities(neuron_distance_csv, data_DIR,algorithms_ordered,metric='neuron_distance_diff',CASE_BY_CASE_PLOT = 0, value_label='Similarity (0~1) based on Average Neuron Distance on Different Structures (D2)')
plt_dist.plot_similarities(neuron_distance_csv, data_DIR,algorithms_ordered,metric='neuron_distance_perc',CASE_BY_CASE_PLOT = 0, value_label='Similarity (0~1) based on Neuron Different Structure Percentage (D3)')



plt_dist.plot_neuron_distance(neuron_distance_csv, data_DIR,algorithms_ordered,CASE_BY_CASE_PLOT = 0)


########################################  for all generated reconstructions ############################
#plot runnign time

output_time_csv = data_DIR + "/running_time_algorithm.csv"

rp.summerize_running_time(time_csv, algorithm_plugin_match_csv,lookup_image_id_table_file,output_time_csv)


df_gold = pd.read_csv(GOLD_CSV)
df_silver_time = pd.read_csv(output_time_csv)
df_share_time = pd.merge(df_silver_time,df_gold,on="image_file_name")
df_share_time.to_csv(data_DIR+ "/running_time_algorithm_gold.csv", index=False)
plt_dist.plot_running_time(data_DIR+ "/running_time_algorithm_gold.csv", data_DIR,algorithms_ordered)
plt_dist.plot_running_time_validation(data_DIR+ "/running_time_algorithm_gold.csv",neuron_distance_csv, data_DIR,algorithms_ordered)
