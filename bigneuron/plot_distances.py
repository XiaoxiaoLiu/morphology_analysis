__author__ = 'xiaoxiaol'

import matplotlib.pyplot as plt
import seaborn as sb
import os
import os.path as path
import numpy as np
import pandas as pd



def plot_neuron_distance(neuron_distance_csv, outputDir,algorithms,CASE_BY_CASE_PLOT = 0):
    #neuron_distance_csv = data_DIR + "/neuron_distance.r.csv"
    #outputDir = data_DIR + "/neuron_dist_plots_r"
    df_nd = pd.read_csv(neuron_distance_csv)
    all_images = np.unique(df_nd.image_file_name)
    if not path.exists(outputDir):
        os.mkdir(outputDir)

    shared_images = []
    if CASE_BY_CASE_PLOT:
        for image in all_images:
            df_image_cur = df_nd[df_nd.image_file_name == image]
            if df_image_cur.shape[0] > 0:
                plt.figure()
                plt.bar(range(df_image_cur.swc_file.size), df_image_cur.neuron_distance)
                plt.xticks(range(df_image_cur.swc_file.size), df_image_cur['algorithm'].values[:],
                           rotation="90")
                plt.ylabel('Neuron Distance')
                plt.subplots_adjust(bottom=0.3)
                plt.savefig(outputDir + '/sorted/' + image.split('/')[-1] + '_nd.png', format='png')
                #plt.show()
                plt.close()
            else:
                print image
            if df_image_cur.shape[0] > 20:
                shared_images.append(image)
        print "there are "+ str(len(shared_images)) +" images that have reconstructions from all " + str(20) +" algorithms."


    ### all algorithm plot
    dfg = df_nd.groupby('algorithm')
    sample_size_per_algorithm=[]
    for alg in algorithms:
        sample_size_per_algorithm.append(dfg.get_group(alg).shape[0])

    plt.figure()
    a=sb.barplot(x='algorithm', y='neuron_distance', data=df_nd,order=algorithms)


    a.set_xticklabels(['%s ($n$=%d )'%(algorithms[i], sample_size_per_algorithm[i]) for i in range(algorithms.size) ])
    #sb.set_context("talk", font_scale=3.0)
    plt.xticks(rotation="90")
    plt.subplots_adjust(bottom=0.5)
    plt.savefig(outputDir + '/all_algorithm_nd_distance.png', format='png')
    plt.show()
    plt.close()






def plot_blasneuron_distance(bn_csv,outputDir,algorithms,CASE_BY_CASE_PLOT=0):

    #outputDir = data_DIR + '/bn_dist'

    df_nd = pd.read_csv(bn_csv)

    if not os.path.exists(outputDir):
        os.mkdir(outputDir)

    images = np.unique(df_nd.image_file_name)
    print "there are "+str(images.size)+" images"


    if CASE_BY_CASE_PLOT:
       for image in images:
            df_image_cur = df_nd[df_nd.image_file_name == image]
            if df_image_cur.shape[0] > 0:
                plt.figure()
                plt.bar(range(df_image_cur.swc_file.size), df_image_cur['SSD'])
                plt.xticks(range(df_image_cur.swc_file.size), df_image_cur['algorithm'].values[:], rotation="60")
                plt.ylabel('BlastNeuron Feature SSD')
                plt.subplots_adjust(bottom=0.7)
                plt.savefig(outputDir + '/' + image.split('/')[-1] + '_bn_dist.png', format='png')
                plt.close()




    myalgorithms = np.unique(df_nd.algorithm)
    if myalgorithms.size !=  algorithms.size:
        print "error: algorithms size is wrong"



    dfg = df_nd.groupby('algorithm')
    sample_size_per_algorithm=[]
    for alg in algorithms:
        sample_size_per_algorithm.append(dfg.get_group(alg).shape[0])

    plt.figure()
    a=sb.barplot(x='algorithm', y='SSD', data=df_nd,order=algorithms)
    #sb.set_context("talk", font_scale=3.0)
    a.set_xticklabels(['%s ($n$=%d )'%(algorithms[i], sample_size_per_algorithm[i]) for i in range(algorithms.size) ])

    #sb.stripplot(y='algorithm', x='SSD', data=df_nd,jitter=True, edgecolor="gray")
    plt.xticks(rotation="90")
    plt.subplots_adjust(bottom=0.5)

    plt.savefig(outputDir + '/all_algorithm_bn_distance.png', format='png')
    plt.show()
    plt.close()




def plot_sample_size(bn_csv,outputDir,algorithms):
    df_nd = pd.read_csv(bn_csv)
    plt.figure()

    sample_size=[]
    print "there are "+str(algorithms.size)+" algorithms"
    for alg in algorithms:
            df_alg = df_nd[df_nd.algorithm == alg]
            sample_size.append(df_alg.image_file_name.size)
    sb.barplot(range(algorithms.size),np.array(sample_size))
    #sb.set_context("paper", font_scale=1.0)
    plt.xticks(range(algorithms.size), algorithms,
               rotation="90")
    plt.subplots_adjust(bottom=0.6,top=0.9)
    plt.ylabel('Number of reconstructions')
    plt.savefig(outputDir + '/all_algorithm_sample_size.png', format='png')
    plt.show()
    plt.close()





def plot_common_set():

    return

