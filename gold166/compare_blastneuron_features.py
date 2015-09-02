# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 16:25:35 2015

@author: xiaoxiaol
"""
import pandas as pd
import os
import numpy as np



def SSD(feature_array1, feature_array2):
    diff_v = np.array(feature_array1) - np.array(feature_array2 ) 
    ssd = np.sum(np.abs(diff_v)**2)

    return ssd


data_DIR =  "/data/mat/xiaoxiaol/data/gold166/preprocessed"
results_csv = "/data/mat/xiaoxiaol/data/gold166/gold_results_combined/sorted/features_with_tags.csv"
gold_csv="/data/mat/xiaoxiaol/data/gold166/checked_final_swcs/preprocessed/features_with_tags.csv"


df_results = pd.read_csv(results_csv)
df_gold = pd.read_csv(gold_csv)


########## filter out spanning tree, which generates huge numbers
df_results = df_results[df_results.algorithm != 'spanningtree']



# find common sets
images_gold = np.unique(df_gold.image)
images_results = np.unique(df_results.image)

# common set of images
final_images = np.intersect1d(images_gold, images_results)




#calcualte std for each image, calculate normalized results
df_normalized = pd.DataFrame()
df_gold_normalized = pd.DataFrame()

for i in range(final_images.size):
      imageName = final_images[i]
      df_image = df_results[df_results.image == imageName]
      if df_image.shape[0] > 5:  # too few samples
          df_gold_image = df_gold[df_gold.image == imageName]
       
          df_image_normalized = pd.DataFrame()
          df_gold_image_normalized = pd.DataFrame()

          cols = df_image.columns[3:]


          for col in cols:
              if (df_image[col].std() == 0.0 ):
                  print "std = 0!", col, " ",imageName
              df_image_normalized[col] = (df_image[col]- df_image[col].median() )/(df_image[col].std()+0.0000001)
              df_gold_image_normalized[col] =( df_gold_image[col] - df_image[col].median() )/(df_image[col].std()+0.0000001)

          df_gold_image_normalized['image'] = imageName
          df_image_normalized[df_image.columns[0:3]] = df_image[df_image.columns[0:3]]

          df_gold_normalized = df_gold_normalized.append(df_gold_image_normalized,ignore_index=True)
          df_normalized = df_normalized.append(df_image_normalized,ignore_index=True)



df_gold_normalized = df_gold_normalized[[ u'image',u'num_nodes', u'soma_surface', u'num_stems', u'num_bifurcations', u'num_branches',
                                         u'num_of_tips', u'overall_width', u'overall_height', u'overall_depth', u'average_diameter',
                                         u'total_length', u'total_surface', u'total_volume', u'max_euclidean_distance',
                                         u'max_path_distance', u'max_branch_order', u'average_contraction', u'average fragmentation',
                                         u'parent_daughter_ratio', u'bifurcation_angle_local', u'bifurcation_angle_remote', u'moment1',
                                         u'moment2', u'moment3', u'moment4', u'moment5', u'moment6', u'moment7', u'moment8', u'moment9',
                                         u'moment10', u'moment11', u'moment12', u'moment13', u'avgR']]

df_normalized = df_normalized[[ u'image', u'algorithm', u'swc_file',u'num_nodes', u'soma_surface', u'num_stems', u'num_bifurcations', u'num_branches',
                                         u'num_of_tips', u'overall_width', u'overall_height', u'overall_depth', u'average_diameter',
                                         u'total_length', u'total_surface', u'total_volume', u'max_euclidean_distance',
                                         u'max_path_distance', u'max_branch_order', u'average_contraction', u'average fragmentation',
                                         u'parent_daughter_ratio', u'bifurcation_angle_local', u'bifurcation_angle_remote', u'moment1',
                                         u'moment2', u'moment3', u'moment4', u'moment5', u'moment6', u'moment7', u'moment8', u'moment9',
                                         u'moment10', u'moment11', u'moment12', u'moment13', u'avgR']]

ssd_metrics=[]
for rowIdx in range(df_normalized.shape[0]):
     imageName = df_normalized.image[rowIdx]
     gold_df = df_gold_normalized[df_gold_normalized.image == imageName]

     #normalize restuls with mean and std
     result= df_normalized.values[rowIdx]
     
     ssd_metrics.append(  SSD(result[3:],gold_df.values[0,1:]) )
     

df_normalized['SSD'] = ssd_metrics
df_normalized.to_csv(data_DIR + '/normalized_w_dist.csv')
