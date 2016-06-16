import sys
import os
import platform


import pandas as pd
import numpy as np



if (platform.system() == "Linux"):
    WORK_PATH = "/local1/xiaoxiaol/work"
else:
    WORK_PATH = "/Users/xiaoxiaoliu/work"

p = WORK_PATH + '/src/morphology_analysis'
sys.path.append(p)

import  bigneuron.median_compare_to_consensus as mc2c


#DATA='/home/xiaoxiaol/work/data'
#DATA='/data/mat/xiaoxiaol/data/big_neuron'
DATA='/mnt/BigNeuron/data'
test = 0


JANELIA_Set1 = 0
JANELIA_Set2 = 1

smooth=1


if JANELIA_Set1:
    data_DIR= DATA+"/Janelia/set1_extract_single"


    fn_list = '~/work/data/jen1_image_file_name_list.csv'
    df_i = pd.read_csv(fn_list)
    imageIDs = df_i['image_file_name']
    if test:
         imageIDs=imageIDs[0:100]

    subfolder="consensus_0330"
    if smooth>0:
           subfolder="consensus_0330_anisosmooth"
    output_dir = data_DIR+'/'+subfolder+"/analysis_results"
    mc2c.pipe(data_DIR+'/'+subfolder, output_dir, imageIDs,'median_distances.csv',COLLECT_FROM_DISTANCE_MATRIX=1,EXTRACT_MEDIAN_CONSENSUS=1)



if JANELIA_Set2:
    data_DIR=DATA+"/Janelia/set2_accepted_single"

    fn_list = '~/work/data/jen2_image_file_name_list.csv'
    df_i = pd.read_csv(fn_list)
    imageIDs = df_i['image_file_name']

    if test:
         imageIDs=imageIDs[0:100]

    subfolder="consensus_0330"
    if smooth>0:
           subfolder="consensus_0330_anisosmooth"
    output_dir = data_DIR+'/'+subfolder+"/analysis_results"
    mc2c.pipe(data_DIR+'/'+subfolder, output_dir, imageIDs,'median_distances.csv',COLLECT_FROM_DISTANCE_MATRIX=1,EXTRACT_MEDIAN_CONSENSUS=1)
    print "\n\n\n"




