__author__ = 'xiaoxiaol'

import matplotlib.pyplot as plt
import seaborn as sb
import os
import os.path as path
import numpy as np
import pandas as pd



algorithm_plugin_match_csv ="/data/mat/xiaoxiaol/data/reconstructions_2015_1214/ported_neuron_tracing_spreadsheet.csv"
df_check_table = pd.read_csv(algorithm_plugin_match_csv)
keys = df_check_table['algorithm']
values = df_check_table['better_algorithm_name']
algorithm_name_mapping = dict(zip(keys, values))



def calculate_similarities(neuron_distance_csv,metric='neuron_distance', output_similarity_csv =None):
    df_nd = pd.read_csv(neuron_distance_csv)
    all_images = np.unique(df_nd.image_file_name)
    all_algorithms = np.unique(df_nd.algorithm)

    print "\n\nCalculate similarity based on " +metric
    print neuron_distance_csv + " has :"
    print str(all_algorithms.size) + " algorithms"
    print str(all_images.size) +" images"
    print all_algorithms

    dfg = df_nd.groupby('image_file_name')

    df_out = pd.DataFrame()
    #sample_size_per_algorithm=[]
    for image in all_images:
        df_image = dfg.get_group(image)
        #sample_size_per_image.append(df_image.shape[0])
        # df_image['similarity'] = np.exp(-(df_image[metric] - df_image[metric].min()))
        # similarity == nan:  metric reports nan
        # similarity = 0 : missing entry ( missing recons)
        df_image['similarity'] = np.exp(-(df_image[metric] - df_nd[metric].min())/(df_nd[metric].max()-df_nd[metric].min()+0.000000001))


        # construct a complete table, and fill the valid results
        df_image_filled_template = pd.DataFrame(columns = df_image.columns)
        df_image_filled_template.algorithm = all_algorithms
        df_image_filled_template.image_file_name = image
        df_image_filled_template['similarity'] = 0.0

        for i in range(df_image.shape[0]):
             alg=   df_image.iloc[i].algorithm
             id = df_image_filled_template[df_image_filled_template.algorithm ==alg].index[0]
             df_image_filled_template.ix[id] = df_image.iloc[i]


        df_out = df_out.append(df_image_filled_template,ignore_index=True)

    if not output_similarity_csv == None:
        df_out.to_csv(output_similarity_csv, index=False)
        print "output "+ output_similarity_csv
    return  df_out




def plot_similarities(neuron_distance_csv, outputDir,algorithms,metric='neuron_distance',CASE_BY_CASE_PLOT = 0, value_label=None):
    #df_nd = pd.read_csv(neuron_distance_csv)

    df_nd = calculate_similarities(neuron_distance_csv,metric,output_similarity_csv=neuron_distance_csv+".similarity.csv")
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
                algorithm_names = [algorithm_name_mapping[x] for x in  df_image_cur['algorithm']]
                plt.xticks(range(df_image_cur.swc_file.size), algorithm_names,
                           rotation="90")
                plt.ylabel(' Similarity (0~1) by ' +metric)
                plt.subplots_adjust(bottom=0.3)
                plt.savefig(outputDir + '/sorted/figs/' + image.split('/')[-1] + '_'+metric+'_similarity.png', format='png')
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
    algorithm_names = [algorithm_name_mapping[x] for x in algorithms]
    a.set_yticklabels(['%s ($n$=%d )'%(algorithm_names[i], rate_per_algorithm[i]) for i in range(algorithms.size) ])
    #sb.set_context("talk", font_scale=3.0)
    #plt.xticks(rotation="90")
    if value_label == None:
        value_label = ' Similarity (0~1) by '+ metric
    plt.ylabel('algorithm (n = # recons)')
    plt.xlabel(value_label)
    plt.subplots_adjust(left=0.4)
    plt.savefig(outputDir + '/all_algorithm_'+metric+'_similarity.png', format='png')
    #plt.show()
    plt.close()

    return



def plot_two_algorithms(neuron_distance_csv, alg1,alg2):
    df_nd = pd.read_csv(neuron_distance_csv)
    dfg = df_nd.groupby('algorithm')
    df_1 = dfg.get_group(alg1)
    df_2 = dfg.get_group(alg2)
    df_merge = pd.merge(df_1,df_2,on="image_file_name")
    plt.figure()
    plt.plot(range(df_merge.image_file_name.size),df_merge['neuron_distance_x'],'r*-')
    plt.plot(range(df_merge.image_file_name.size),df_merge['neuron_distance_y'],"b*-")
    plt.legend()
    plt.show()
    plt.close()



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

                sb.barplot(y=range(df_image_cur['algorithm'].size),x=df_image_cur['neuron_distance'],orient="h")
                #plt.xticks(range(df_image_cur.swc_file.size), df_image_cur['algorithm'].values[:], rotation="90")
                algorithm_names = [algorithm_name_mapping[x] for x in df_image_cur['algorithm']]
                plt.yticks(range(df_image_cur['algorithm'].size), np.array(algorithm_names))

                plt.xlabel('Average Neuron Distance (s1)')
                plt.subplots_adjust(left=0.4,right=0.9)

                plt.savefig(outputDir + '/sorted/figs/' + image.split('/')[-1] + '_nd.png', format='png')

                plt.close()


            else:
                print image
        #     if df_image_cur.shape[0] > 10:
        #         shared_images.append(image)
        # print "there are "+ str(len(shared_images)) +" images that have reconstructions from all " + str(20) +" algorithms."
        #

    ### all algorithm plot
    dfg = df_nd.groupby('algorithm')
    sample_size_per_algorithm=[]
    for alg in algorithms:
        sample_size_per_algorithm.append(dfg.get_group(alg).shape[0])



    #plot the average node distances
    plt.figure()
    a=sb.barplot(y='algorithm', x='neuron_distance', data=df_nd,order=algorithms)
    algorithm_names = [algorithm_name_mapping[x] for x in algorithms]
    a.set_yticklabels(['%s ($n$=%d )'%(algorithm_names[i], sample_size_per_algorithm[i]) for i in range(algorithms.size) ])
    #sb.set_context("talk", font_scale=3.0)
    #plt.xticks(rotation="90")
    plt.xlabel('Average Neuron Distance (s1)')
    plt.subplots_adjust(left=0.4)
    plt.savefig(outputDir + '/all_algorithm_average_neuron_distance_s1.png', format='png')
    plt.show()
    plt.close()

    #plot the differences
    plt.figure()
    #df_nd['neuron_difference'] = df_nd['neuron_distance_diff'] *df_nd['neuron_distance_perc']
    a=sb.barplot(y='algorithm', x='neuron_difference', data=df_nd,order=algorithms,orient='h')
    algorithm_names = [algorithm_name_mapping[x] for x in algorithms]
    a.set_yticklabels(['%s ($n$=%d )'%(algorithm_names[i], sample_size_per_algorithm[i]) for i in range(algorithms.size) ])
    #sb.set_context("talk", font_scale=3.0)
    #plt.xticks(rotation="90")
    plt.xlabel('Neuron Difference Score (s2*s3)')
    plt.subplots_adjust(left=0.4)
    plt.savefig(outputDir + '/all_algorithm_neuron_difference_score_s2s3.png', format='png')
    #plt.show()
    plt.close()
    return



def plot_blasneuron_distance(bn_csv,outputDir,algorithms,CASE_BY_CASE_PLOT=0):

    #outputDir = data_DIR + '/bn_dist'
    df_nd = pd.read_csv(bn_csv)

    if not os.path.exists(outputDir):
        os.mkdir(outputDir)

    images = np.unique(df_nd.image_file_name)
    print "Global Morph feature plotting: there are "+str(images.size)+" images"


    if CASE_BY_CASE_PLOT:
       for image in images:
            df_image_cur = df_nd[df_nd.image_file_name == image]
            if df_image_cur.shape[0] > 0:
                plt.figure()

                sb.barplot(y=range(df_image_cur['algorithm'].size),x=df_image_cur['SSD'],orient="h")
                algorithm_names = [algorithm_name_mapping[x] for x in df_image_cur['algorithm']]
                plt.yticks(range(df_image_cur['algorithm'].size), np.array(algorithm_names))

                plt.xlabel('Global Morph Feature SSD')
                plt.subplots_adjust(left=0.4,right=0.9)
                plt.savefig(outputDir + '/sorted/figs/' + image.split('/')[-1] + '_bn_dist.png', format='png')

                plt.close()


    myalgorithms = np.unique(df_nd.algorithm)
    if myalgorithms.size !=  algorithms.size:
        print "error: algorithms size is wrong"
        print myalgorithms

    dfg = df_nd.groupby('algorithm')
    sample_size_per_algorithm=[]
    for alg in algorithms:
        sample_size_per_algorithm.append(dfg.get_group(alg).shape[0])

    plt.figure()
    a=sb.barplot(y='algorithm', x='SSD', data=df_nd,order=algorithms,orient="h")
    #sb.set_context("talk", font_scale=3.0)
    algorithm_names = [algorithm_name_mapping[x] for x in algorithms]
    a.set_yticklabels(['%s ($n$=%d )'%(algorithm_names[i], sample_size_per_algorithm[i]) for i in range(algorithms.size) ])

    #sb.stripplot(y='algorithm', x='SSD', data=df_nd,jitter=True, edgecolor="gray")
    #plt.xticks(rotation="90")
    plt.xlabel('Global Morph Feature Score Distance')
    plt.subplots_adjust(left=0.4)

    plt.savefig(outputDir + '/all_algorithm_glogbal_morph_feature_ssd.png', format='png')
    #plt.show()
    plt.close()



def plot_sample_size(bn_csv,outputDir,algorithms):
    df_nd = pd.read_csv(bn_csv)

    df_nd['algorithm'] = [algorithm_name_mapping[x] for x in df_nd['algorithm'] ]

    plt.figure()

    sample_size=[]
    print "there are "+str(algorithms.size)+" algorithms"
    for alg in algorithms:
            df_alg = df_nd[df_nd.algorithm == alg]
            sample_size.append(df_alg.image_file_name.size)
    sb.barplot(y=range(algorithms.size),x=np.array(sample_size),orient="h")
    #sb.set_context("paper", font_scale=1.0)

    algorithm_names = [algorithm_name_mapping[x] for x in algorithms]
    plt.yticks(range(algorithms.size), algorithm_names)
    plt.subplots_adjust(left=0.4,right=0.9)
    plt.xlabel('Number of reconstructions')
    plt.savefig(outputDir + '/all_algorithm_valid_reconstruction_number.png', format='png')
    #plt.show()
    plt.close()
    return



def plot_running_time(time_csv, outputDir, algorithms):

    df_time = pd.read_csv(time_csv)
    m_algorithms = np.unique(df_time.algorithm)



    df_time['running_time'] =  df_time['running_time'] /1000.0
    dfg = df_time.groupby('algorithm')

    sample_size_per_algorithm=[]
    for alg in algorithms:
        #print alg
        if alg in m_algorithms:
            sample_size_per_algorithm.append(dfg.get_group(alg).shape[0])
        else:
            sample_size_per_algorithm.append(0)

    plt.figure()
    a=sb.barplot(y='algorithm', x='running_time', data=df_time,order=algorithms,orient="h")

    algorithm_names = [algorithm_name_mapping[x] for x in algorithms]

    a.set_yticklabels(['%s ($n$=%d )'%(algorithm_names[i], sample_size_per_algorithm[i]) for i in range(algorithms.size) ])
    plt.subplots_adjust(left=0.4,right=0.9)
    plt.xlabel('Running Time (seconds)')
    plt.savefig(outputDir + '/runningtime_goldset.png', format='png')
    #plt.show()
    plt.close()


    return

    #plt.show()
    plt.close()
    return



def plot_running_time_validation(time_csv,neuron_distance_csv, outputDir, algorithms):


    df_nd = pd.read_csv(neuron_distance_csv)
    df_nd_group_by_image = df_nd.groupby('image_file_name')


    df_time = pd.read_csv(time_csv)

    df_time['running_time'] =  df_time['running_time'] /1000.0


    dfg = df_time.groupby('image_file_name')
    all_images = np.unique(df_time['image_file_name'])

    df_out = pd.DataFrame()
    #sample_size_per_algorithm=[]
    for image in all_images:
        df_image = dfg.get_group(image)
        #valid results
        df_nd_image = df_nd_group_by_image.get_group(image)

        # construct a complete table, and fill the valid results
        df_time_filled_template = pd.DataFrame(columns = df_time.columns)
        df_time_filled_template.algorithm = algorithms
        df_time_filled_template.image_file_name = image
        df_time_filled_template['running_time'] = 60*60

        for i in range(df_image.shape[0]):

             alg =   df_image.iloc[i]['algorithm']

             if alg in np.unique(df_nd_image['algorithm']):
                 id = df_time_filled_template[df_time_filled_template.algorithm == alg].index[0]
                 df_time_filled_template.ix[id] = df_image.iloc[i]

        df_out = df_out.append(df_time_filled_template,ignore_index=True)

    a=sb.barplot(y='algorithm', x='running_time', data=df_out, order=algorithms, orient="h")

    algorithm_names = [algorithm_name_mapping[x] for x in algorithms]

    a.set_yticklabels(algorithm_names)
    plt.subplots_adjust(left=0.4,right=0.9)
    plt.xlabel('Running Time (seconds): all 163 images, 1 hour wall time is used for N/A entries')
    plt.savefig(outputDir + '/runningtime_goldset_1hourForNA.png', format='png')
    #plt.show()
    plt.close()




    return



def plot_common_set():

    return

