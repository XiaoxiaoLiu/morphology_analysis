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


data_DIR ="/data/mat/xiaoxiaol/data/20151130_rhea_reconstructions_for_allen300_silver_set"
original_dir = data_DIR +"/auto_recons"
preprocessed_dir = data_DIR+ "/resampled"
sorted_dir = data_DIR +"/sorted"



###### generate gold standard 166 table
gold_dir = "/data/mat/xiaoxiaol/data/gold166/checked_final_swcs/preprocessed"
GOLD_CSV = "/data/mat/xiaoxiaol/data/gold166/gold.csv"
#rp.gen_gold_csv(gold_dir, GOLD_CSV)
gold_feature_csv= gold_dir +"/features_with_tags.csv"
#generate gold standard blastneuron features
# data_DIR= WORK_PATH+"/data/gold79/sorted"
 #    #########################################################
 #    FEATURE_FILE = data_DIR + '/features.nfb'
 #    nfb.generateALLFeatureCSV_gold166(FEATURE_FILE, gold_feature_csv')



######  resample and sort
#rp.resample_and_sort(original_dir,preprocessed_dir,sorted_dir)



###### genearte sliver data table
SILVER_CSV = data_DIR+'/recon_table.csv'
lookup_image_id_table_file = data_DIR +"/image_name_lookup_table.csv"
#rp.recon_table_gen(sorted_dir,lookup_image_id_table_file,SILVER_CSV)


##### merge to get the common set between gold and silver
merged_csv_file = data_DIR+'/shared.csv'
#rp.merge_gold_silver(GOLD_CSV,SILVER_CSV,merged_csv_file)


# run with commented computation!!!!!!!!!!!
neuron_distance_csv=data_DIR+'/nd.csv'
rp.cal_neuron_dist(merged_csv_file,neuron_distance_csv )



algorithms = np.unique(df_nd.algorithm)
#plot
plt_dist.plot_neuron_distance(neuron_distance_csv, data_DIR,algorithms,CASE_BY_CASE_PLOT = 0)
# collect neuron distance


###### compute blastneuron features
output_feature_csv = sorted_dir +'/features_with_tags.csv'
#rp.cal_bn_features(sorted_dir,output_feature_csv)
bn_csv = data_DIR+"/blastneuron_ssd.csv"
rp.cal_blastneuron_distance(output_feature_csv,gold_feature_csv,merged_csv_file,bn_csv)


# plot
plt_dist.plot_blasneuron_distance(bn_csv,data_DIR,algorithms,CASE_BY_CASE_PLOT=0)



