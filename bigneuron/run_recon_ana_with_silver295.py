__author__ = 'xiaoxiaol'
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


##################    silver standard dataset processing #########################################
###### only need to do it once
silver_dir = "/data/mat/xiaoxiaol/data/big_neuron/silver/silver_295_swcs"
SILVER_CSV = "/data/mat/xiaoxiaol/data/big_neuron/silver/silver_295_swcs/silver.csv"  #manually generated
silver_feature_csv= "/data/mat/xiaoxiaol/data/big_neuron/silver/silver_295_swcs/features_with_tags.csv"

#print "resample and sort silver gt 295"
##rp.resample_and_sort(silver_dir+"/original",silver_dir+"/resampled",silver_dir+"/sorted")
#rp.gen_gold_feature_csv(silver_dir+"/original",output_gold_csv_file= SILVER_CSV,output_gold_feature_csv=silver_feature_csv)

algorithm_plugin_match_csv = "/data/mat/xiaoxiaol/data/big_neuron/silver/ported_neuron_tracing_spreadsheet.csv"


###########################  preprocessing and organize data ##################################
data_DIR ="/data/mat/xiaoxiaol/data/big_neuron/silver/20160113_merged_silver_gt"
original_dir = data_DIR +"/auto_recons"
resampled_dir = data_DIR+ "/resampled"
sorted_dir = data_DIR +"/sorted"
lookup_image_id_table_file = data_DIR +"/../image_name_lookup_table.csv"
time_csv = data_DIR + "/auto_recons/running_time_merged.csv"
#########################################
ND = 0
BD = 0
COMPUTE = 0
neuron_distance_csv = data_DIR +'/neuron_distances_with_silver_gt.csv'
######  resample  #?and sort
#rp.resample_and_sort(original_dir,resampled_dir,sorted_dir,GEN_QSUB=0,overwrite_sorted=0)


if COMPUTE>0:

    ######  sliver data table
    rc_SILVER_CSV = data_DIR+'/recon_table.csv'
    rp.recon_table_gen(original_dir,lookup_image_id_table_file,rc_SILVER_CSV)


    #####  merge to get the common set between silver gt and silver rc
    merged_csv_file = data_DIR+'/recon_shared_with_silver_gt.csv'
    rp.merge_gold_silver(SILVER_CSV,rc_SILVER_CSV,merged_csv_file)


    #####  report which gold dataset did not have any recons?
    df_merge = pd.read_csv(merged_csv_file)
    df_silver_gt = pd.read_csv(SILVER_CSV)
    m = pd.unique(df_merge.image_file_name)
    g = pd.unique(df_silver_gt.image_file_name)

    print "\n\nSilver gt dataset contains " +str(g.size) +" image dataset"
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
    if BD:
        tmp_feature_csv = original_dir +'/tmp_features_with_tags.csv'

        rp.cal_bn_features(original_dir,tmp_feature_csv)
        output_feature_csv = data_DIR +'/features_with_tags.csv'
        rp.map_image_name(tmp_feature_csv,lookup_image_id_table_file, output_feature_csv)

        df_m = pd.read_csv(output_feature_csv)
        df_m1 = df_m[df_m['num_nodes'] != 0]
        df_m1.to_csv(output_feature_csv, index=False)


        bn_csv = data_DIR+"/blastneuron_lm_ssd.csv"
        rp.cal_blastneuron_distance(output_feature_csv,silver_feature_csv,merged_csv_file,output_csv=bn_csv, LMEASURE_ONLY = 1)



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
if BD:
   plt_dist.plot_blasneuron_distance(bn_csv,data_DIR,algorithms_ordered,CASE_BY_CASE_PLOT=1)
   plt_dist.plot_similarities(bn_csv, data_DIR,algorithms_ordered,metric='SSD',CASE_BY_CASE_PLOT = 0,value_label='Similarity (0~1) on Global Morph Feature Score')



plt_dist.plot_similarities(neuron_distance_csv, data_DIR,algorithms_ordered,metric='neuron_distance',CASE_BY_CASE_PLOT = 0,value_label='Similarity (0~1) on Average Neuron Distance(s1)')
plt_dist.plot_similarities(neuron_distance_csv, data_DIR,algorithms_ordered,metric='neuron_difference',CASE_BY_CASE_PLOT = 0, value_label='Similarity (0~1) on Neuron Difference Score (s2*s3)')
plt_dist.plot_similarities(neuron_distance_csv, data_DIR,algorithms_ordered,metric='neuron_distance_diff',CASE_BY_CASE_PLOT = 0, value_label='Similarity (0~1) on Average Neuron Distance on Different Structures (s2)')
plt_dist.plot_similarities(neuron_distance_csv, data_DIR,algorithms_ordered,metric='neuron_distance_perc',CASE_BY_CASE_PLOT = 0, value_label='Similarity (0~1) on Neuron Different Structure Percentage (s3)')



plt_dist.plot_neuron_distance(neuron_distance_csv, data_DIR,algorithms_ordered,CASE_BY_CASE_PLOT = 0)


########################################  for all generated reconstructions ############################
#plot runnign time
df_silver_gt = pd.read_csv(SILVER_CSV)

output_time_csv = data_DIR + "/running_time_algorithm.csv"

rp.summerize_running_time(time_csv, algorithm_plugin_match_csv,lookup_image_id_table_file,output_time_csv)


df_rc_time = pd.read_csv(output_time_csv)
df_share_time = pd.merge(df_rc_time,df_silver_gt,on="image_file_name")
df_share_time.to_csv(data_DIR+ "/running_time_algorithm_silver_gt.csv", index=False)
plt_dist.plot_running_time(data_DIR+ "/running_time_algorithm_silver_gt.csv", data_DIR,algorithms_ordered)
plt_dist.plot_running_time_validation(data_DIR+ "/running_time_algorithm_silver_gt.csv",neuron_distance_csv, data_DIR,algorithms_ordered)
