__author__ = 'xiaoxiaoliu'


import pandas as pd
import numpy as np
import os
from os import sys, path
import seaborn as sb


import matplotlib.pyplot as plt
WORK_PATH = "/Users/xiaoxiaoliu/work"
p =  WORK_PATH + '/src/morphology_analysis'
sys.path.append(p)
sb.set_context("poster")


data_DIR = WORK_PATH+"/data/20151030_rhea_reconstructions_for_allen300_silver_set"
original_dir = data_DIR + "/auto_recons"
preprocessed_dir = data_DIR +"/79/resampled"
sorted_dir = data_DIR +"/79/sorted"
gold_image_dir = WORK_PATH+"/data/gold79/origin_data"
votemaps_dir = data_DIR+"/votemaps"


######  plot##################################################
NEURON_DISTANCE_PLOT = 1
df_nd = pd.read_csv(data_DIR + "/neuron_distance.r.csv")
outputDir = data_DIR +"/neuron_distance_plots"
algorithms = np.unique(df_nd.algorithm)

if not path.exists(outputDir):
    os.mkdir(outputDir)
if NEURON_DISTANCE_PLOT:
    CASE_BY_CASE_PLOT =0
    if CASE_BY_CASE_PLOT:
        images = np.unique(df_nd.gold_swc_file)
        for image in images:
            df_image_cur = df_nd[df_nd.gold_swc_file == image]
            if df_image_cur.shape[0] > 10:
                plt.figure()
                plt.bar(range(df_image_cur.swc_file.size), df_image_cur.neuron_distance)
                plt.xticks(range(df_image_cur.swc_file.size), df_image_cur['algorithm'].values[:],
                           rotation="90")
                plt.ylabel('Neuron Distance')
                plt.subplots_adjust(bottom=0.3)
                plt.savefig(outputDir + '/' + image.split('/')[-1] + '_nd.png', format='png')
                #plt.show()
                plt.close()
    ### all algorithm plot
    plt.figure()
    sb.barplot(x='algorithm', y='neuron_distance', data=df_nd,order=algorithms)
    sb.set_context("talk", font_scale=3.0)
    plt.xticks(rotation="90")
    plt.subplots_adjust(bottom=0.5)
    plt.savefig(outputDir + '/all_algorithm_distance.png', format='png')
    plt.show()
    plt.close()






sb.set_context('poster')
NPEARSON_PLOT=0
if NPEARSON_PLOT:
    df_pr = pd.read_csv(votemaps_dir+"/pearsonr.csv")
    #plt.figure(figsize=(10,30))
    plt.subplots_adjust(bottom=0.4,top=0.95)
    sb.tsplot(data=df_pr.perasonr)

    image_name_short=[]
    for i in range(df_pr.image.size):
        image_name_short.append(df_pr.image[i][0:15])
        # st=df_pr.gold_image_file[i].split('/')[-3]
        # st=st.split('checked')[-1]
        # image_name_short.append( st[2:])

    plt.xticks(range(df_pr.image.size),image_name_short,rotation=90)
    plt.show()
    plt.savefig(votemaps_dir + '/pearsonr.png', format='png')
    plt.close()




BLASTNEURON_PLOT=1
if BLASTNEURON_PLOT:
    #####################  plots################
    df_nd = pd.read_csv(data_DIR+"/normalized_bn_dist.csv")

    outputDir = data_DIR + '/bn_dist'
    if not os.path.exists(outputDir):
        os.mkdir(outputDir)

    CASE_BY_CASE_PLOT = 0
    images = np.unique(df_nd.image)
    print "there are "+str(images.size)+" images"
    if CASE_BY_CASE_PLOT:

        for image in images:
            df_image_cur = df_nd[df_nd.image == image]
            if df_image_cur.shape[0] > 0:
                plt.figure()
                plt.bar(range(df_image_cur.swc_file.size), df_image_cur['SSD'])
                plt.xticks(range(df_image_cur.swc_file.size), df_image_cur['algorithm'].values[:], rotation="60")
                plt.ylabel('BlastNeuron Feature SSD')
                plt.subplots_adjust(bottom=0.7)
                plt.savefig(outputDir + '/' + image.split('/')[-1] + '_bn_dist.png', format='png')
                plt.close()

    ALL_ALGORITHM_PLOT = 1

    if ALL_ALGORITHM_PLOT:
        plt.figure()
        sb.barplot(x='algorithm', y='SSD', data=df_nd,order=algorithms)

        #sb.stripplot(y='algorithm', x='SSD', data=df_nd,jitter=True, edgecolor="gray")
        plt.xticks(rotation="90")
        plt.subplots_adjust(bottom=0.4)

        #plt.show()
        plt.savefig(outputDir + '/all_algorithm_distance.png', format='png')
        plt.close()

        plt.figure()

        sample_size=[]
        print "there are "+str(algorithms.size)+" algorithms"
        for alg in algorithms:
                df_alg = df_nd[df_nd.algorithm == alg]
                sample_size.append(df_alg.image.size)
        sb.barplot(range(algorithms.size),np.array(sample_size))
        plt.xticks(range(algorithms.size), algorithms,
                   rotation="90")
        plt.subplots_adjust(bottom=0.6,top=0.9)
        plt.ylabel('Number of reconstructions')
        plt.savefig(outputDir + '/all_algorithm_sample_size.png', format='png')
        plt.show()
        plt.close()