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


subfolder="0401_gold163_all_soma_sort"
#subfolder="gold_163_all_soma_sort_0328"
data_DIR="/data/mat/xiaoxiaol/data/big_neuron/silver/"+subfolder

lookup_image_id_table_file = data_DIR +"/../image_name_lookup_table.csv"

#neuron_distance_csv = data_DIR +'/../20160113_merged_gold_gt/neuron_distances_with_gold.csv'
df_extract_median_consensus_distance = pd.read_csv(data_DIR+'/analysis_results/extracted_median_consensus.csv')

merged_csv= data_DIR+'/consensus_median_algorithm_all_nd_ave.csv'
EXTRACT_NUMBERS=1

if EXTRACT_NUMBERS == 1:
    dfg_t=df_extract_median_consensus_distance.groupby('algorithm')
    dfg_median =dfg_t.get_group('median')
    image_names=[]
    df_lookup_table = pd.read_csv(lookup_image_id_table_file)
    for i in range(len(dfg_median)):
        image_id = dfg_median.iloc[i]['image_id']
        image = df_lookup_table.image_file_name[image_id-1]
        image_names.append(image)
    dfg_median['image_file_name'] = image_names
    dfg_median.to_csv(data_DIR+'/analysis_results/median_with_image_file_name.csv')


    #compare valid consensus results with each individual algorithms
    log_file_list= data_DIR + "/consensus_weighted_dist_log_list.txt"
    median_log_file_list = data_DIR + "/median_dist_log_list.txt"
    os.system('ls '+data_DIR+'/*/consensus.eswc.weighted.dist.log >'+log_file_list)
    # read all weighted_neuron_distance logs for consensus results
    output_csv = data_DIR + "/consensus_weighted_dist.csv"
    rp.collect_consensus_distance(log_file_list, output_csv, lookup_image_id_table_file)




    log_file_list= data_DIR + "/soma_sorted_dist_log_list.txt"
    os.system('ls '+data_DIR+'/*/auto_recons/*.swc.dist.log >'+log_file_list)
    neuron_distance_csv = data_DIR + "/soma_sorted_neuron_dist.csv"
    #rp.collect_neuron_distance(log_file_list, neuron_distance_csv, lookup_image_id_table_file)


    #read nd distance csv
    df_nd = pd.read_csv(neuron_distance_csv)
    before = df_nd.shape[0]
    df_nd = df_nd[df_nd['neuron_distance_12'] != 0]  # neuron distance bug, one node swc will have nd=0
    df_nd = df_nd[df_nd['neuron_distance_12'] != -1]  # neuron distance bug, one node swc will have nd=0
    after = df_nd.shape[0]
    if after < before:
        print "\n\nwarning: removing 0 nd entries:", after-before

    df_nd_g = df_nd.groupby('image_file_name')


    # not all images have generated consensus results
    df_consensus_wd = pd.read_csv(output_csv)
    images_have_consensus_results = np.unique(df_consensus_wd['image_file_name'])
    images = np.unique(df_nd['image_file_name'])
    print len(images_have_consensus_results)," images that have consensus results"

    #print out images that do not have consensus results
    topdirs = glob.glob(os.path.join(data_DIR, '*'))
        # print topdirs
    print "the following images have no consensus results:"
    for subdir in topdirs:
            #print subdir
            if os.path.isdir(subdir) and  ( not os.path.exists(subdir+"/consensus.eswc")):
                print subdir


    # compose a spreadsheet for comparison
    df_merge = pd.DataFrame(columns=['image_file_name', 'algorithm', 'weighted_ave_neuron_distance'])
    i=0
    for image in images:#images_have_consensus_results:
         #print image
         df_nd_image = df_nd_g.get_group(image)
         num_rows = df_nd_image.shape[0]
         #print num_rows
         for j in range(num_rows):
             df_merge.loc[i] = [image, df_nd_image.iloc[j]['algorithm'], df_nd_image.iloc[j]['neuron_distance_12'] ]
             #df_merge.loc[i] = [image, df_nd_image.iloc[j]['algorithm'], df_nd_image.iloc[j]['neuron_distance_ave'] ]
             i= i+1
         df_con_matching = df_consensus_wd[df_consensus_wd['image_file_name']==image]


         median_swc_file = dfg_median[dfg_median['image_file_name'] ==image]
         if len(median_swc_file) ==0:
             continue
         median_alg = rp.matchFileToAlgorithmName( median_swc_file.iloc[0]['swc_file_name'])
         df_median = df_nd_image[df_nd_image['algorithm']== median_alg]
         median_nd = df_median.iloc[0]['neuron_distance_12']

         if df_con_matching.shape[0] >0 :
             #consensus_wd_t = df_con_matching.iloc[0]['weighted_neuron_distance_ave']
             consensus_wd_t = df_con_matching.iloc[0]['weighted_neuron_distance_12']
             #print consensus_wd_t
             df_merge.loc[i] = [image,'consensus',consensus_wd_t]
             i= i+1
             df_merge.loc[i] = [image,'median',median_nd]
             i= i+1

    df_merge.to_csv(merged_csv, index=False)





# plot
# ## sort by sample size
df_merge =pd.read_csv(merged_csv)
algorithms = np.unique(df_merge.algorithm)
dfg = df_merge.groupby('algorithm')
sample_size_per_algorithm = np.zeros(algorithms.size)
for i in range( algorithms.size):
    sample_size_per_algorithm[i] = (dfg.get_group(algorithms[i]).shape[0])

order = sample_size_per_algorithm.argsort()
algorithms_ordered = algorithms[order[::-1]]





#error bar plot
plt_dist.plot_compare_consensus_distance(merged_csv, data_DIR,algorithms_ordered,metric='weighted_ave_neuron_distance',CASE_BY_CASE_PLOT = 0,
                                         value_label='Weighted Average Neuron Distance 12 (to the Gold Standard)')
#plt_dist.plot_similarities(merged_csv, data_DIR,algorithms_ordered,metric='weighted_ave_neuron_distance',CASE_BY_CASE_PLOT = 0,
#                                        value_label='Similarities on  Weighted Average Neuron Distance (to the Gold Standard)')

df_merge_m_c=df_merge[ (df_merge['algorithm'] == "median") |( df_merge['algorithm'] == "consensus")]
plt_dist.plot_compare_median_consensus(output_dir=data_DIR,df_order= df_merge_m_c, metric='weighted_ave_neuron_distance', type = 'ts',DISPLAY = 1)