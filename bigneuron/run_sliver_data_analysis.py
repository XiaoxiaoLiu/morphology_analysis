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



data_DIR ="/data/mat/xiaoxiaol/data/20151130_rhea_reconstructions_for_allen300_silver_set"
original_dir = data_DIR +"/auto_recons"
preprocessed_dir = data_DIR+ "/resampled"
sorted_dir = data_DIR +"/sorted"



###### generate gold standard 166 table
gold_dir = "/data/mat/xiaoxiaol/data/gold166/checked_final_swcs/preprocessed"
GOLD_CSV = "/data/mat/xiaoxiaol/data/gold166/gold.csv"
#rp.gen_gold_csv(gold_dir, GOLD_CSV)



######  resample and sort
#rp.resample_and_sort(original_dir,preprocessed_dir,sorted_dir)



###### genearte sliver data table
SILVER_CSV = data_DIR+'/recon_table.csv'
lookup_image_id_table_file = data_DIR +"/image_name_lookup_table.csv"
#rp.recon_table_gen(sorted_dir,lookup_image_id_table_file,SILVER_CSV)


#####merge to get the common set between gold and silver
merged_csv_file = data_DIR+'/shared.csv'
#rp.merge_gold_silver(GOLD_CSV,SILVER_CSV,merged_csv_file)

rp.cal_neuron_dist(merged_csv_file)


