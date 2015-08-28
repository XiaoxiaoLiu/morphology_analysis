# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 17:18:00 2015

@author: xiaoxiaol
"""


import pandas as pd
import numpy as np
from os import path, sys

from os import sys, path
p = path.dirname(path.dirname(path.abspath(__file__)))
sys.path.append(p)
print sys.path

sys.path.append(p+'/blast_neuron')

import blast_neuron_comp as bn




data_DIR =  "/data/mat/xiaoxiaol/data/gold166/preprocessed"
results_csv = "/data/mat/xiaoxiaol/data/gold166/preprocessed/features_with_tags.csv"
gold_csv="/data/mat/xiaoxiaol/data/gold166/checked_final_swcs/preprocessed/features_with_tags.csv"


df_results = pd.read_csv(results_csv)
df_gold = pd.read_csv(gold_csv)



# find common sets
images_gold = np.unique(df_gold.image)
images_results = np.unique(df_results.image)

# common set of images
final_images = np.intersect1d(images_gold, images_results)


df_neuron_distance = pd.DataFrame( columns=('swc_file', 'gold_swc_file', 'neuron_distance'))
idx = 0
for i in range(final_images.size):
      imageName = final_images[i]
      df_image = df_results[df_results.image == imageName]

      df_gold_image = df_gold[df_gold.image == imageName].iloc[0]

      for j in range(df_image.shape[0]):
          df_swc = df_image.iloc[j]

          nd = bn.neuron_dist(df_swc.swc_file, df_gold_image.swc_file, df_swc.swc_file+".log")
          df_neuron_distance.loc[idx] = [df_swc.swc_file, df_gold_image.swc_file, nd['ave'] ]
          idx = idx +1

df_neuron_distance.to_csv(data_DIR +"/neuron_distance.csv")

