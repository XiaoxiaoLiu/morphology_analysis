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

import  bigneuron.median_compare_to_consensus as mc2c

import pandas as pd



subfolder="0401_gold163_all_soma_sort"
#subfolder="gold_163_all_soma_sort_0328"
data_DIR="/data/mat/xiaoxiaol/data/big_neuron/silver/"+subfolder


df_nd = pd.read_csv(data_DIR+'/list.txt')
imageIDs = df_nd['image_id'].apply(str)


mc2c.pipe(input_data_dir=data_DIR, output_dir=data_DIR+"/analysis_results", imageIDs=imageIDs,distance_file_postfix='median_distances.csv', COLLECT_FROM_DISTANCE_MATRIX=1,EXTRACT_MEDIAN_CONSENSUS=1, DISPLAY=0)
print "\n\n\n"



