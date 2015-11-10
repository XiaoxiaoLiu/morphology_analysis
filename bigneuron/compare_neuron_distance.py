# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 17:18:00 2015

@author: xiaoxiaol
"""

import pandas as pd
import numpy as np
import os
from os import path, sys

from os import sys, path

p = path.dirname(path.dirname(path.abspath(__file__)))
sys.path.append(p)
print sys.path

sys.path.append(p + '/blast_neuron')

import blast_neuron_comp as bn
import matplotlib.pyplot as plt
import seaborn as sb


data_DIR = "/data/mat/xiaoxiaol/data/gold166/gold166_results_combined/sorted"
results_csv = data_DIR + "/features_with_tags.csv"
gold_csv = "/data/mat/xiaoxiaol/data/gold166/checked_final_swcs/preprocessed/features_with_tags.csv"

df_results = pd.read_csv(results_csv)
df_gold = pd.read_csv(gold_csv)



# find common sets
images_gold = np.unique(df_gold.image)
images_results = np.unique(df_results.image)

# common set of images
final_images = np.intersect1d(images_gold, images_results)
print "There are ", final_images.size, " images reconstructed in total"

#################################################

RUN_DISTANCE_CALC = 0
if RUN_DISTANCE_CALC:
    for i in range(final_images.size):
        imageName = final_images[i]
        df_image = df_results[df_results.image == imageName]
        df_gold_image = df_gold[df_gold.image == imageName].iloc[0]

        for j in range(df_image.shape[0]):
            df_swc = df_image.iloc[j]
            # bn.run_neuron_dist(df_swc.swc_file, df_gold_image.swc_file, df_swc.swc_file+".log", 1, "nd")
            bn.run_neuron_dist(df_gold_image.swc_file, df_swc.swc_file, df_swc.swc_file + ".r.log", 1, "nd")


###############################################
COLLECT_NEURON_DISTANCE = 1
if COLLECT_NEURON_DISTANCE:
    df_neuron_distance = pd.DataFrame(columns=('swc_file', 'gold_swc_file', 'algorithm', 'neuron_distance'))
    idx = 0
    for i in range(final_images.size):
        imageName = final_images[i]
        df_image = df_results[df_results.image == imageName]
        df_gold_image = df_gold[df_gold.image == imageName].iloc[0]

        for j in range(df_image.shape[0]):
            df_swc = df_image.iloc[j]

            logfile = df_swc.swc_file + ".r.log"
            if path.isfile(logfile):
                nd = bn.read_neuron_dist_log(logfile)
                algorithm = df_swc.algorithm
                df_neuron_distance.loc[idx] = [df_swc.swc_file, df_gold_image.swc_file, algorithm, nd['ave']]
                idx = idx + 1
    df_neuron_distance.to_csv(data_DIR + "/neuron_distance.r.csv")


######  plot##################################################
df_nd = pd.read_csv(data_DIR + "/neuron_distance.r.csv")
outputDir = data_DIR + "/neuron_dist_plots_r"
if not path.exists(outputDir):
    os.mkdir(outputDir)

CASE_BY_CASE_PLOT = 1
if CASE_BY_CASE_PLOT:
    images = np.unique(df_nd.gold_swc_file)
    for image in images:
        df_image_cur = df_nd[df_nd.gold_swc_file == image]
        if df_image_cur.shape[0] > 10:
            plt.figure()
            plt.bar(range(df_image_cur.swc_file.size), df_image_cur.neuron_distance)
            plt.xticks(range(df_image_cur.swc_file.size), df_image_cur['algorithm'].values[:], rotation="60")
            plt.ylabel('Neuron Distance')
            plt.subplots_adjust(bottom=0.3)
            plt.savefig(outputDir + '/' + image.split('/')[-1] + '_nd.png', format='png')
            plt.close()

ALL_ALGORITHM_PLOT = 1
if ALL_ALGORITHM_PLOT:
    plt.figure()
    sb.barplot(x='algorithm', y='neuron_distance', data=df_nd)
    sb.set_context("talk", font_scale=3.0)
    plt.xticks(rotation="80")
    plt.subplots_adjust(bottom=0.5)
    plt.savefig(outputDir + '/all_algorithm_distance.png', format='png')
    plt.close()

