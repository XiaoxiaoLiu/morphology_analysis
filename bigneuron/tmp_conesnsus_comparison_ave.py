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
import glob


data_DIR ="/data/mat/xiaoxiaol/data/big_neuron/silver/gold_163_soma_sort_0210"
lookup_image_id_table_file = data_DIR +"/../image_name_lookup_table.csv"

neuron_distance_csv = data_DIR +'/../20160113_merged_gold_gt/neuron_distances_with_gold.csv'


#compare valid consensus results with each individual algorithms
log_file_list= data_DIR + "/consensus_weighted_dist_log_list.txt"
#os.system('ls '+data_DIR+'/*/processed/consensus_p2_dark_pruned_2.eswc.weighted.dist.log >'+log_file_list)
# read all weighted_neuron_distance logs for consensus results
output_csv = data_DIR + "/consensus_weighted_dist.csv"
#rp.collect_consensus_distance(log_file_list, output_csv, lookup_image_id_table_file)






#read nd distance csv
df_nd = pd.read_csv(neuron_distance_csv)
before = df_nd.shape[0]
df_nd = df_nd[df_nd['neuron_distance'] != 0]  # neuron distance bug, one node swc will have nd=0
after = df_nd.shape[0]
if after < before:
    print "\n\nwarning: removing 0 nd entries:", after-before

df_nd_g = df_nd.groupby('image_file_name')


# not all images have generated consensus results
df_consensus_wd = pd.read_csv(output_csv)
all_images = np.unique(df_nd['image_file_name'])
images_have_consensus_results = np.unique(df_consensus_wd['image_file_name'])
print len(images_have_consensus_results)," images that have consensus results"

#print out images that do not have consensus results
topdirs = glob.glob(os.path.join(data_DIR, '*'))
    # print topdirs
print "the following images have no consensus results:"
for subdir in topdirs:
        #print subdir
        if os.path.isdir(subdir) and  ( not os.path.exists(subdir+"/consensus_p2_dark_pruned_2.eswc")):
            print subdir


# compose a spreadsheet for comparison
df_merge = pd.DataFrame(columns=['image_file_name', 'algorithm', 'weighted_ave_neuron_distance'])
i=0
for image in images_have_consensus_results:
     #print image
     df_nd_image = df_nd_g.get_group(image)

     num_rows = df_nd_image.shape[0]
     #print num_rows
     for j in range(num_rows):
         df_merge.loc[i] = [image, df_nd_image.iloc[j]['algorithm'], df_nd_image.iloc[j]['neuron_distance'] ]
         i= i+1
     df_con_matching = df_consensus_wd[df_consensus_wd['image_file_name']==image]
     if df_con_matching.shape[0] >0 :
         consensus_wd_t = df_con_matching.iloc[0]['weighted_neuron_distance_ave']
         df_merge.loc[i] = [image,'consensus',consensus_wd_t]
         i= i+1



merged_csv= data_DIR+'/consensus_compare_wnd.csv'
df_merge.to_csv(merged_csv, index=False)

# plot
# ## sort by sample size
algorithms = np.unique(df_merge.algorithm)

dfg = df_merge.groupby('algorithm')
sample_size_per_algorithm = np.zeros(algorithms.size)
for i in range( algorithms.size):
    sample_size_per_algorithm[i] = (dfg.get_group(algorithms[i]).shape[0])

order = sample_size_per_algorithm.argsort()
algorithms_ordered = algorithms[order[::-1]]





plt_dist.plot_compare_consensus_distance(merged_csv, data_DIR,algorithms_ordered,metric='weighted_ave_neuron_distance',CASE_BY_CASE_PLOT = 0,
                                         value_label='Weighted Average Neuron Distance (bidirectional)')
plt_dist.plot_similarities(merged_csv, data_DIR,algorithms_ordered,metric='weighted_ave_neuron_distance',CASE_BY_CASE_PLOT = 0,
                                        value_label='Similarities on  Weighted Average Neuron Distance (bidirectional)')
