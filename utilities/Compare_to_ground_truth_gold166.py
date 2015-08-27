# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 16:25:35 2015

@author: xiaoxiaol
"""
import pandas as pd
import os
import numpy as np


#V3D="/Users/xiaoxiaoliu/work/v3d/v3d_external/bin/vaa3d64.app/Contents/MacOS/vaa3d64"
V3D="/local1/xiaoxiaol/work/v3d/v3d_external/bin/vaa3d"

def neuron_dist(inputswc_path1, inputswc_path2, logfile):
#Distance between neuron 1 /home/xiaoxiaol/work/data/test_frog2-2.swc and neuron 2 /home/xiaoxiaol/work/data/test_frog2-2.swc is: 
#entire-structure-average = 8.20009e-07
#differen-structure-average = 0
#percent of different-structure = 0

# log file format
# file1 file2   8.20009e-07  0 0

    cmd = V3D +  " -x neuron_distance -f neuron_distance -i "+ inputswc_path1 +" " + inputswc_path2 +" -o "+logfile
    os.system(cmd)
    print cmd
    
    # read log file
    f = open(logfile, 'r')    
    line=f.readline()
    
    diff = line.split(' ')[2]  # read the average difference number
    return  diff
    
    
    
    
    

def SSD(feature_array1, feature_array2):
    diff_v = np.array(feature_array1) - np.array(feature_array2 ) 
    ssd = np.sum(np.abs(diff_v)**2)

    return ssd




results_csv = "/data/mat/xiaoxiaol/data/gold166/preprocessed/features_with_tags.csv"
gold_csv="/data/mat/xiaoxiaol/data/gold166/checked_final_swcs/preprocessed/features_with_tags.csv"


df_results = pd.read_csv(results_csv)
df_gold = pd.read_csv(gold_csv)


# find common sets
images_gold=np.unique(df_gold.image)
images_results = np.unique(df_results.image)

# common set of images
final_images = np.intersect1d(images_gold, images_results)
df_filtered = df_results[df_results.image.isin(final_images)]

#calcualte std for each image, calculate normalized results
df_normalized = pd.DataFrame()


for i in range(final_images.size):
      imageName = final_images[i]
      df_image = df_filtered[df_filtered.image == imageName]
      gld_features = df_gold.loc[df_gold.image == imageName
      



for rowIdx in range(df_results_normalized.shape[0]):
     imageName = df_results_normalized.image[rowIdx]
     gold_features = df_gold.loc[df_gold.image == imageName]

     #normalize restuls with mean and std
     results_features = df_results_normalized.values[rowIdx]
     
     
     difference = SSD(results_features[3:],gold_features.values[0,2:])
     

     
     