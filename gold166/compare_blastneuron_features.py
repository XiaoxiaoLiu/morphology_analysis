# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 16:25:35 2015

@author: xiaoxiaol
"""
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb


def SSD(feature_array1, feature_array2):
    diff_v = np.array(feature_array1) - np.array(feature_array2)
    ssd = np.sum(np.abs(diff_v) ** 2)

    return ssd


data_DIR = "/data/mat/xiaoxiaol/data/gold166/gold166_results_combined/sorted"
results_csv = data_DIR + "/features_with_tags.csv"
gold_csv = "/data/mat/xiaoxiaol/data/gold166/checked_final_swcs/preprocessed/features_with_tags.csv"

df_results = pd.read_csv(results_csv)
df_gold = pd.read_csv(gold_csv)


########## filter out spanning tree, which generates huge numbers
# df_results = df_results[df_results.algorithm != 'spanningtree']


# find common sets
images_gold = np.unique(df_gold.image)
images_results = np.unique(df_results.image)

# common set of images
final_images = np.intersect1d(images_gold, images_results)

feature_cols = [u'num_nodes', u'soma_surface', u'num_stems', u'num_bifurcations', u'num_branches',\
                u'num_of_tips', u'overall_width', u'overall_height', u'overall_depth', u'average_diameter',\
                u'total_length', u'total_surface', u'total_volume', u'max_euclidean_distance',\
                u'max_path_distance', u'max_branch_order', u'average_contraction', u'average fragmentation',\
                u'parent_daughter_ratio', u'bifurcation_angle_local', u'bifurcation_angle_remote', u'moment1',\
                u'moment2', u'moment3', u'moment4', u'moment5', u'moment6', u'moment7', u'moment8', u'moment9',\
                u'moment10', u'moment11', u'moment12', u'moment13', u'avgR']

selected_cols = [u'num_nodes', u'num_bifurcations', u'num_branches', u'num_of_tips', u'overall_width',u'overall_height', u'overall_depth',u'total_length', u'average fragmentation',
u'bifurcation_angle_local', u'bifurcation_angle_remote']#, u'moment1',u'moment2', u'moment3', u'moment4', u'moment5', u'moment6', u'moment7', u'moment8', u'moment9',u'moment10', u'moment11', u'moment12', u'moment13']

my_cols1 = [u'image', u'algorithm', u'swc_file']
my_cols1.extend(selected_cols)
df_results_s = df_results[my_cols1]
del df_results
my_cols2 = [u'image']
my_cols2.extend(selected_cols)
df_gold_s = df_gold[my_cols2]
del df_gold

# calcualte std for each image, calculate normalized results
df_normalized = pd.DataFrame(columns=my_cols1)
df_gold_normalized = pd.DataFrame(columns=my_cols2)

for i in range(final_images.size):
    imageName = final_images[i]
    df_image = df_results_s[df_results_s.image == imageName]

    if df_image.shape[0] > 5:  # too few samples
        df_image[selected_cols] = df_image[selected_cols].astype(float)  # some moment values are interpreted as strings
        df_gold_image = df_gold_s[df_gold_s.image == imageName]

        df_normalized_per_image = pd.DataFrame(columns=my_cols1)
        df_gold_normalized_per_image = pd.DataFrame(columns=my_cols2)

        df_normalized_per_image[[u'image', u'algorithm', u'swc_file']] = df_image[[u'image', u'algorithm', u'swc_file']]
        df_gold_normalized_per_image[u'image'] = df_gold_image.image
        print imageName

        for col in selected_cols:
            if (df_image[col].std() == 0.0 ):
                print "std = 0!", col, " ", imageName
            df_normalized_per_image[col] = (df_image[col] - df_image[col].median() ) / (df_image[col].std() + 0.0000001)
            df_gold_normalized_per_image[col] = ( df_gold_image[col] - df_image[col].median() ) / (df_image[col].std() + 0.0000001)

        # append the results for this image
        df_normalized = df_normalized.append(df_normalized_per_image, ignore_index=True)
        df_gold_normalized= df_gold_normalized.append(df_gold_normalized_per_image, ignore_index=True)


#calculated sit
ssd_metrics = []
for rowIdx in range(df_normalized.shape[0]):
    imageName = df_normalized.image[rowIdx]
    gold_df = df_gold_normalized.loc[df_gold_normalized.image == imageName]

    #normalize restuls with mean and std
    result = df_normalized.iloc[rowIdx]
    sumsquare = SSD(result[selected_cols], gold_df[selected_cols])
    ssd_metrics.append(sumsquare)

df_normalized['SSD'] = ssd_metrics
## reordering columns
df_normalized.to_csv(data_DIR + '/normalized_bn_dist.csv')



#####################  plots################
df_nd = df_normalized

outputDir = data_DIR + '/bn_dist'
if not os.path.exists(outputDir):
    os.mkdir(outputDir)

CASE_BY_CASE_PLOT = 1
sb.set_context("talk", font_scale=1.5)

if CASE_BY_CASE_PLOT:
    images = np.unique(df_nd.image)
    for image in images:
        df_image_cur = df_nd[df_nd.image == image]
        if df_image_cur.shape[0] > 0:
            plt.figure()
            plt.bar(range(df_image_cur.swc_file.size), df_image_cur['SSD'])
            plt.xticks(range(df_image_cur.swc_file.size), df_image_cur['algorithm'].values[:], rotation="60")
            plt.ylabel('BlastNeuron Feature SSD')
            plt.subplots_adjust(bottom=0.3)
            plt.savefig(outputDir + '/' + image.split('/')[-1] + '_bn_dist.png', format='png')
            plt.close()

ALL_ALGORITHM_PLOT = 1
if ALL_ALGORITHM_PLOT:
    plt.figure()
    sb.barplot(x='algorithm', y='SSD', data=df_nd)

    plt.xticks(rotation="80")
    plt.subplots_adjust(bottom=0.5)
    plt.savefig(outputDir + '/all_algorithm_distance.png', format='png')
    plt.close()


print "done plotting"