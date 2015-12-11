__author__ = 'xiaoxiaol'

import matplotlib.pyplot as plt
import seaborn as sb
import os
import os.path as path
import numpy as np
import pandas as pd




def calculate_similarities(neuron_distance_csv,metric='neuron_distance', output_similarity_csv =None):
    df_nd = pd.read_csv(neuron_distance_csv)
    all_images = np.unique(df_nd.image_file_name)
    all_algorithms = np.unique(df_nd.algorithm)

    dfg = df_nd.groupby('image_file_name')

    df_out = pd.DataFrame()
    #sample_size_per_algorithm=[]
    for image in all_images:
        df_image = dfg.get_group(image)
        #sample_size_per_image.append(df_image.shape[0])
        # df_image['similarity'] = np.exp(-(df_image[metric] - df_image[metric].min()))
        df_image['similarity'] = np.exp(-(df_image[metric] - df_nd[metric].min())/(df_nd[metric].max()-df_nd[metric].min()+0.000000001))


        # construct a complete table, and fill the valid reults
        df_image_filled_template = pd.DataFrame(columns = df_image.columns)
        df_image_filled_template.algorithm = all_algorithms
        df_image_filled_template.image_file_name = image
        df_image_filled_template['similarity'] = 0.0

        for i in range(df_image.shape[0]):
             alg=   df_image.iloc[i].algorithm
             id = df_image_filled_template[df_image_filled_template.algorithm ==alg].index[0]
             df_image_filled_template.ix[id] = df_image.iloc[i]


        df_out = df_out.append(df_image_filled_template)

    if not output_similarity_csv == None:
        df_out.to_csv(output_similarity_csv)
    return  df_out




def plot_similarities(neuron_distance_csv, outputDir,algorithms,metric='neuron_distance',CASE_BY_CASE_PLOT = 0, value_label=None):
    #df_nd = pd.read_csv(neuron_distance_csv)

    df_nd = calculate_similarities(neuron_distance_csv,metric)
    all_images = np.unique(df_nd.image_file_name)
    if not path.exists(outputDir):
        os.mkdir(outputDir)

    if CASE_BY_CASE_PLOT:
        dfg = df_nd.groupby('image_file_name')

        #sample_size_per_algorithm=[]
        for image in all_images:
            df_image_cur = dfg.get_group(image)
            if df_image_cur.shape[0] > 0:
                plt.figure()
                plt.bar(range(df_image_cur.swc_file.size), df_image_cur['similarity'])
                plt.xticks(range(df_image_cur.swc_file.size), df_image_cur['algorithm'].values[:],
                           rotation="90")
                plt.ylabel(' Similarity (0~1) by ' +metric)
                plt.subplots_adjust(bottom=0.3)
                plt.savefig(outputDir + '/sorted/' + image.split('/')[-1] + '_'+metric+'_similarity.png', format='png')
                #plt.show()
                plt.close()
            else:
                print image+" has no valid reconstructions"


    dfg = df_nd.groupby('algorithm')
    rate_per_algorithm=[]
    for alg in algorithms:
        df_a = dfg.get_group(alg)
        sucessrate= float(np.count_nonzero(df_a['similarity']))/df_a.shape[0] * 100
        number = np.count_nonzero(df_a['similarity'])
        rate_per_algorithm.append(number)


    plt.figure()
    a=sb.barplot(y='algorithm', x='similarity', data=df_nd,order=algorithms)
    a.set_yticklabels(['%s ($n$=%d )'%(algorithms[i], rate_per_algorithm[i]) for i in range(algorithms.size) ])
    #sb.set_context("talk", font_scale=3.0)
    #plt.xticks(rotation="90")
    if value_label == None:
        value_label = ' Similarity (0~1) by '+ metric
    plt.ylabel('algorithm (n = # recons)')
    plt.xlabel(value_label)
    plt.subplots_adjust(left=0.4)
    plt.savefig(outputDir + '/all_algorithm_'+metric+'_similarity.png', format='png')
    plt.show()
    plt.close()

    return




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


    #plot the average node distances
    plt.figure()
    a=sb.barplot(y='algorithm', x='neuron_distance', data=df_nd,order=algorithms)
    a.set_yticklabels(['%s ($n$=%d )'%(algorithms[i], sample_size_per_algorithm[i]) for i in range(algorithms.size) ])
    #sb.set_context("talk", font_scale=3.0)
    #plt.xticks(rotation="90")
    plt.xlabel('Average Neuron Distance')
    plt.subplots_adjust(left=0.4)
    plt.savefig(outputDir + '/all_algorithm_average_neuron_distance.png', format='png')
    plt.show()
    plt.close()

    #plot the differences
    plt.figure()
    #df_nd['neuron_difference'] = df_nd['neuron_distance_diff'] *df_nd['neuron_distance_perc']
    a=sb.barplot(y='algorithm', x='neuron_difference', data=df_nd,order=algorithms,orient='h')
    a.set_yticklabels(['%s ($n$=%d )'%(algorithms[i], sample_size_per_algorithm[i]) for i in range(algorithms.size) ])
    #sb.set_context("talk", font_scale=3.0)
    #plt.xticks(rotation="90")
    plt.xlabel('Neuron Difference Score')
    plt.subplots_adjust(left=0.4)
    plt.savefig(outputDir + '/all_algorithm_neuron_difference.png', format='png')
    plt.show()
    plt.close()
    return



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
                plt.savefig(outputDir + '/sorted/' + image.split('/')[-1] + '_bn_dist.png', format='png')
                plt.close()

    myalgorithms = np.unique(df_nd.algorithm)
    if myalgorithms.size !=  algorithms.size:
        print "error: algorithms size is wrong"
        print myalgorithms.size
        print algorithms.size

    print myalgorithms

    dfg = df_nd.groupby('algorithm')
    sample_size_per_algorithm=[]
    for alg in algorithms:
        sample_size_per_algorithm.append(dfg.get_group(alg).shape[0])

    plt.figure()
    a=sb.barplot(y='algorithm', x='SSD', data=df_nd,order=algorithms,orient="h")
    #sb.set_context("talk", font_scale=3.0)
    a.set_yticklabels(['%s ($n$=%d )'%(algorithms[i], sample_size_per_algorithm[i]) for i in range(algorithms.size) ])

    #sb.stripplot(y='algorithm', x='SSD', data=df_nd,jitter=True, edgecolor="gray")
    #plt.xticks(rotation="90")
    plt.xlabel('BlastNeuron Feature Score Distance')
    plt.subplots_adjust(left=0.4)

    plt.savefig(outputDir + '/all_algorithm_blastneuron_feature_ssd.png', format='png')
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
    sb.barplot(y=range(algorithms.size),x=np.array(sample_size),orient="h")
    #sb.set_context("paper", font_scale=1.0)
    plt.yticks(range(algorithms.size), algorithms)
    plt.subplots_adjust(left=0.4,right=0.9)
    plt.xlabel('Number of reconstructions')
    plt.savefig(outputDir + '/all_algorithm_valid_reconstruction_number.png', format='png')
    plt.show()
    plt.close()



def plot_common_set():

    return

