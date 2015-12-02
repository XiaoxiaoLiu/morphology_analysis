__author__ = 'xiaoxiaol'

import matplotlib.pyplot as plt
import seaborn as sb
import os
import sys.path as path
import numpy as np
import pandas as pd





def plot_neuron_distance(neuron_distance_csv, outputDir,CASE_BY_CASE_PLOT = 0,ALL_ALGORITHM_PLOT =1):
    #neuron_distance_csv = data_DIR + "/neuron_distance.r.csv"
    #outputDir = data_DIR + "/neuron_dist_plots_r"
    df_nd = pd.read_csv(neuron_distance_csv)

    if not path.exists(outputDir):
        os.mkdir(outputDir)


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

    if ALL_ALGORITHM_PLOT:
        plt.figure()
        sb.barplot(x='algorithm', y='neuron_distance', data=df_nd)
        sb.set_context("talk", font_scale=3.0)
        plt.xticks(rotation="80")
        plt.subplots_adjust(bottom=0.5)
        plt.savefig(outputDir + '/all_algorithm_distance.png', format='png')
        plt.close()

    return




def plot_blasneuron_distance(bn_csv,outputDir,CASE_BY_CASE_PLOT=0,ALL_ALGORITHM_PLOT=1):
    #outputDir = data_DIR + '/bn_dist'
    df_nd = pd.read_csv(bn_csv)

    if not os.path.exists(outputDir):
        os.mkdir(outputDir)

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


    if ALL_ALGORITHM_PLOT:
        plt.figure()
        sb.barplot(x='algorithm', y='SSD', data=df_nd)

        plt.xticks(rotation="80")
        plt.subplots_adjust(bottom=0.5)
        plt.savefig(outputDir + '/all_algorithm_distance.png', format='png')
        plt.close()


    print "done plotting"