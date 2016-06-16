__author__ = 'xiaoxiaol'

import matplotlib.pyplot as plt
import seaborn as sb
import os
import os.path as path
import numpy as np
import pandas as pd


data_DIR =  "/data/mat/xiaoxiaol/data/big_neuron/silver"

algorithm_plugin_match_csv = data_DIR + "/ported_neuron_tracing_spreadsheet.csv"

df_check_table = pd.read_csv(algorithm_plugin_match_csv)
keys = df_check_table['algorithm']
values = df_check_table['better_algorithm_name']
algorithm_name_mapping = dict(zip(keys, values))

sb.set_context("talk", font_scale=0.7)

def calculate_similarities(neuron_distance_csv,metric='neuron_distance', output_similarity_csv =None):
    df_nd = pd.read_csv(neuron_distance_csv)
    all_images = np.unique(df_nd.image_file_name)
    all_algorithms = np.unique(df_nd.algorithm)

    print "\n\nCalculate similarity based on " +metric
    print neuron_distance_csv + " has :"
    print str(all_algorithms.size) + " algorithms"
    print str(all_images.size) +" images"
    #print all_algorithms

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
        #df_image['similarity'] = np.exp(-df_image[metric]/df_nd[metric].max())


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
    df_nd_ori = pd.read_csv(neuron_distance_csv)

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
                #plt.subplots_adjust(bottom=0.3)
                plt.savefig(outputDir + '/sorted/figs/' + image.split('/')[-1] + '_'+metric+'_similarity.png', format='png')
                #plt.show()
                plt.close()
            else:
                print image+" has no valid reconstructions"


    # dfg = df_nd.groupby('algorithm')
    # rate_per_algorithm=[]
    # for alg in algorithms:
    #     df_a = dfg.get_group(alg)
    #     sucessrate= float(np.count_nonzero(df_a['similarity']))/df_a.shape[0] * 100
    #     number = np.count_nonzero(df_a['similarity'])
    #     rate_per_algorithm.append(number)
    dfg = df_nd_ori.groupby('algorithm')
    sample_size_per_algorithm=[]
    for alg in algorithms:
        sample_size_per_algorithm.append(dfg.get_group(alg).shape[0])

    plt.figure()
    sb.set_context("talk", font_scale=0.7)
    a=sb.barplot(y='algorithm', x='similarity', data=df_nd,order=algorithms)
    algorithm_names = [algorithm_name_mapping[x] for x in algorithms]
    a.set_yticklabels(['%s ($n$=%d )'%(algorithm_names[i], sample_size_per_algorithm[i]) for i in range(algorithms.size) ])
    #sb.set_context("talk", font_scale=3.0)
    #plt.xticks(rotation="90")
    if value_label == None:
        value_label = ' Similarity (0~1) by '+ metric
    plt.ylabel('algorithm (n = # recons)')
    plt.xlabel(value_label)
    plt.subplots_adjust(left=0.4,right=0.9, bottom=0.1, top=0.9)
    plt.savefig(outputDir + '/'+value_label+'.png', format='png')
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
    #plt.show()
    plt.close()





def plot_compare_consensus_distance2(distance_csv, outputDir,algorithms,metric, CASE_BY_CASE_PLOT = 0,value_label=None):
    df_nd = pd.read_csv(distance_csv)
    all_images = np.unique(df_nd.image_file_name)
    if not path.exists(outputDir):
        os.mkdir(outputDir)


    ### all algorithm plot
    dfg = df_nd.groupby('algorithm')
    sample_size_per_algorithm=[]
    for alg in algorithms:
        sample_size_per_algorithm.append(dfg.get_group(alg).shape[0])



    #plot the weighted average node distances
    plt.figure()
    sb.set_context("talk", font_scale=0.7)
   # a=sb.barplot(y='algorithm', x=metric, data=df_nd,order=algorithms)

    algorithm_names = [algorithm_name_mapping[x] for x in algorithms]
    a.set_yticklabels(['%s ($n$=%d )'%(algorithm_names[i], sample_size_per_algorithm[i]) for i in range(algorithms.size) ])
    #sb.set_context("talk", font_scale=3.0)
    #plt.xticks(rotation="90")
    plt.xlabel(value_label)
    plt.subplots_adjust(left=0.4, bottom=0.1, top=0.9)
    plt.savefig(outputDir + '/'+value_label+'.png', format='png')
    #plt.show()
    plt.close()

    return


def TBD_plot_compare_consensus_distance_box_plot(distance_csv, outputDir,algorithms,metric, CASE_BY_CASE_PLOT = 0,value_label=None):
    df_nd = pd.read_csv(distance_csv)
    all_images = np.unique(df_nd.image_file_name)
    if not path.exists(outputDir):
        os.mkdir(outputDir)

    shared_images = []
    # if CASE_BY_CASE_PLOT:
    #     for image in all_images:
    #         df_image_cur = df_nd[df_nd.image_file_name == image]
    #         if df_image_cur.shape[0] > 0:
    #             plt.figure()
    #             ax = sb.boxplot(x=metric, y="algorithm", data=df_image_cur,notch=True,whis=1.5, color="c")
    #
    #             # Add in points to show each observation
    #             sb.stripplot(x=metric, y="algorithm", data=df_image_cur,jitter=False, size=2, color=".3", linewidth=0)
    #             ax.set_xscale("log")
    #             sb.despine(trim=True)
    #
    #
    #             algorithm_names = [algorithm_name_mapping[x] for x in df_image_cur['algorithm']]
    #             plt.yticks(range(df_image_cur['algorithm'].size), np.array(algorithm_names))
    #             plt.xlabel(value_label)
    #             #plt.subplots_adjust(left=0.4, bottom=0.1, top=0.9)
    #             plt.savefig(outputDir + '/sorted/figs/' + image.split('/')[-1] + '_nd_boxplot.png', format='png')
    #
    #             plt.close()
    #
    #
    #         else:
    #             print image
    #     #     if df_image_cur.shape[0] > 10:
    #     #         shared_images.append(image)
    #     # print "there are "+ str(len(shared_images)) +" images that have reconstructions from all " + str(20) +" algorithms."
    #     #

    ### all algorithm plot
    dfg = df_nd.groupby('algorithm')
    sample_size_per_algorithm=[]
    for alg in algorithms:
        sample_size_per_algorithm.append(dfg.get_group(alg).shape[0])


    #plot the weighted average node distances
    plt.figure()
    sb.set_context("talk", font_scale=0.7)

    ax = sb.boxplot(x=metric, y="algorithm", data=df_nd,whis=1.5, notch=True,color="c")

    # Add in points to show each observation
    sb.stripplot(x=metric, y="algorithm", data=df_nd,jitter=False, size=2, color=".3", linewidth=0)

    sb.despine(trim=True)

    algorithm_names = [algorithm_name_mapping[x] for x in algorithms]
    ax.set_yticklabels(['%s ($n$=%d )'%(algorithm_names[i], sample_size_per_algorithm[i]) for i in range(algorithms.size) ])
    #sb.set_context("talk", font_scale=3.0)
    #plt.xticks(rotation="90")
    plt.xlabel(value_label)
    plt.subplots_adjust(left=0.5, bottom=0.1, top=0.9)
    plt.savefig(outputDir + '/'+value_label+'.box_plot.png', format='png')
    plt.show()
    plt.close()

    return



def plot_compare_median_consensus(output_dir, df_order, metric,DISPLAY = 0):
   TYPE = 'diff'


   if DISPLAY:
        plt.figure()


   if TYPE == 'box':
       ax = sb.boxplot(x=metric, y="algorithm", data=df_order,
                 whis=1.5, color="c")

       # Add in points to show each observation
       sb.stripplot(x=metric, y="algorithm", data=df_order,
                jitter=True, size=3, color=".3", linewidth=0)
       ax.set_xscale("log")
       sb.despine(trim=True)
       # plt.xlabel('images sorted by the average neuron distance of the median reconstruction')
       plt.savefig(output_dir + '/compare_median_with_consensus_'+metric+'.png', format='png')

       if DISPLAY:
             plt.show()
       plt.close()



   df_order = df_order.sort([metric], ascending=[0])


   if DISPLAY:
         plt.figure()

   if TYPE =='ts':
        #sb.factorplot(x='order', y=metric, hue="algorithm", data=df_order, ci=None, kind="point",join=True)
        sb.pointplot(x='image_file_name', y=metric,data = df_order,hue="algorithm", join=True)
        #sb.pointplot(x='image_file_name', y=metric,data = df_consensus, join=True)

        plt.xlabel('images sorted by the average neuron distance of the median reconstruction')
        plt.xticks(rotation='vertical')
        plt.subplots_adjust(bottom=0.5, right=0.8, top=0.9)

        plt.savefig(output_dir + '/compare_median_with_consensus_'+metric+'_ts.png', format='png')

   if TYPE =='diff':
        df_median = df_order[df_order['algorithm']== 'median']
        df_consensus = df_order[df_order['algorithm']== 'consensus']
        df_median.sort(['image_file_name'], ascending=[0], inplace=True)
        df_consensus.sort(['image_file_name'], ascending=[0], inplace=True)

        df_diff=pd.DataFrame(columns=[metric, 'image_file_name'])

        y = np.array(df_consensus[metric].values - df_median[metric].values)
        #df_diff['image_file_name'] = df_consensus['image_file_name']
        #df_diff.sort([metric], ascending=[0], inplace=True)
        f_big = len(np.nonzero(y<0)[0])

        print "consensus is closer to each reconstructions than the median reconstructions in %.2f percent of the %d total images"  %( (100.0* f_big)/len(y), len(y))

        #sb.pointplot(x='image_file_name', y=metric, data = df_diff, join=True)
        y.sort()
        plt.plot(range(len(y)),y,"b.-")
        plt.plot(range(len(y)), np.zeros(len(y)),'r-')
        plt.xlabel('images sorted by the average neuron distance of the median reconstruction')
        plt.ylabel("d(conensus, gs) - d(median, gs) ")

        plt.savefig(output_dir + '/compare_median_with_consensus_'+metric+'_diff.png', format='png')

        print "investigate the following cases, where consensus is worse than median:"
        for im in df_median['image_file_name']:
              dif = df_median[df_median['image_file_name'] == im].iloc[0][metric] - df_consensus[df_consensus['image_file_name'] == im].iloc[0][metric]
              if  dif < -5.0:
                   print im +"  ,with a diff = " + str(dif)



   if DISPLAY:
        plt.show()
        plt.close()






def plot_compare_consensus_distance(distance_csv, outputDir,algorithms,metric, CASE_BY_CASE_PLOT = 0,value_label=None):
    df_nd = pd.read_csv(distance_csv)
    all_images = np.unique(df_nd.image_file_name)
    if not path.exists(outputDir):
        os.mkdir(outputDir)

    shared_images = []
    if CASE_BY_CASE_PLOT:
        for image in all_images:
            df_image_cur = df_nd[df_nd.image_file_name == image]
            if df_image_cur.shape[0] > 0:
                plt.figure()

                df_image_cur.sort([metric], ascending=[1], inplace=True)
                sb.barplot(y=range(df_image_cur['algorithm'].size),x=df_image_cur[metric],orient="h")
                #plt.xticks(range(df_image_cur.swc_file.size), df_image_cur['algorithm'].values[:], rotation="90")
                algorithm_names = [algorithm_name_mapping[x] for x in df_image_cur['algorithm']]
                for index, item in enumerate(algorithm_names):
                    if (item == "Consensus"):
                            algorithm_names[index] = "==>Consensus"
                    if (item == "median"):
                            algorithm_names[index] = "==>Median"
                plt.yticks(range(df_image_cur['algorithm'].size), np.array(algorithm_names))

                plt.xlabel(value_label)
                plt.subplots_adjust(left=0.4, bottom=0.1, top=0.9)
                plt.savefig(outputDir + '/case_by_case/' + image.split('/')[-1] + '_nd.png', format='png')
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


    #plot the weighted average node distances
    plt.figure()
    sb.set_context("talk", font_scale=0.7)
    a=sb.barplot(y='algorithm', x=metric, data=df_nd,order=algorithms)
    algorithm_names = [algorithm_name_mapping[x] for x in algorithms]
    for index, item in enumerate(algorithm_names):
                    if (item == "Consensus"):
                            algorithm_names[index] = "==>Consensus"
                    if (item == "median"):
                            algorithm_names[index] = "==>Median"
    a.set_yticklabels(['%s ($n$=%d )'%(algorithm_names[i], sample_size_per_algorithm[i]) for i in range(algorithms.size) ])
    #sb.set_context("talk", font_scale=3.0)
    #plt.xticks(rotation="90")
    plt.xlabel(value_label)
    plt.subplots_adjust(left=0.4, bottom=0.1, top=0.9)
    plt.savefig(outputDir + '/'+value_label+'.png', format='png')
    #plt.show()
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

                sb.barplot(y=range(df_image_cur['algorithm'].size),x=df_image_cur['neuron_distance'],orient="h")
                #plt.xticks(range(df_image_cur.swc_file.size), df_image_cur['algorithm'].values[:], rotation="90")
                algorithm_names = [algorithm_name_mapping[x] for x in df_image_cur['algorithm']]
                plt.yticks(range(df_image_cur['algorithm'].size), np.array(algorithm_names))

                plt.xlabel('Average Neuron Distance (D1)')
                plt.subplots_adjust(left=0.4, bottom=0.1, top=0.9)
                plt.savefig(outputDir + '/sorted/figs/' + image.split('/')[-1] + '_nd.png', format='png')

                plt.close()


            else:
                print image
        #     if df_image_cur.shape[0] > 10:
        #         shared_images.append(image)
        # print "there are "+ str(len(shared_images)) +" images that have reconstructions from all " + str(20) +" algorithms."
        #
    if algorithms== None:
        algorithms = np.unique(df_nd.algorithm)
    ### all algorithm plot
    dfg = df_nd.groupby('algorithm')
    sample_size_per_algorithm=[]
    for alg in algorithms:
        sample_size_per_algorithm.append(dfg.get_group(alg).shape[0])



    #plot the average node distances
    plt.figure()
    sb.set_context("talk", font_scale=0.7)
    a=sb.barplot(y='algorithm', x='neuron_distance_ave', data=df_nd,order=algorithms)
    algorithm_names = [algorithm_name_mapping[x] for x in algorithms]
    a.set_yticklabels(['%s ($n$=%d )'%(algorithm_names[i], sample_size_per_algorithm[i]) for i in range(algorithms.size) ])
    #sb.set_context("talk", font_scale=3.0)
    #plt.xticks(rotation="90")
    plt.xlabel('Average Neuron Distance')
    plt.subplots_adjust(left=0.4, bottom=0.1, top=0.9)
    plt.savefig(outputDir + '/average_neuron_distance.png', format='png')
    #plt.show()
    plt.close()

    # #plot the differences
    # plt.figure()
    # sb.set_context("talk", font_scale=0.7)
    # #df_nd['neuron_difference'] = df_nd['neuron_distance_diff'] *df_nd['neuron_distance_perc']
    # a=sb.barplot(y='algorithm', x='neuron_difference', data=df_nd,order=algorithms,orient='h')
    # algorithm_names = [algorithm_name_mapping[x] for x in algorithms]
    # a.set_yticklabels(['%s ($n$=%d )'%(algorithm_names[i], sample_size_per_algorithm[i]) for i in range(algorithms.size) ])
    # #sb.set_context("talk", font_scale=3.0)
    # #plt.xticks(rotation="90")
    # plt.xlabel('Neuron Difference Score (D2*D3)')
    # plt.subplots_adjust(left=0.4,right=0.9, bottom=0.1, top=0.9)
    # plt.savefig(outputDir + '/neuron_difference_score_D2D3.png', format='png')
    # #plt.show()
    # plt.close()
    return



def plot_blasneuron_distance(bn_csv,outputDir,algorithms=None,CASE_BY_CASE_PLOT=0):

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
                plt.subplots_adjust(left=0.4,right=0.7)
                plt.savefig(outputDir + '/sorted/figs/' + image.split('/')[-1] + '_bn_dist.png', format='png')

                plt.close()

    if algorithms== None:
        algorithms = np.unique(df_nd.algorithm)


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
    plt.subplots_adjust(left=0.4,right=0.9, bottom=0.1, top=0.9)

    plt.savefig(outputDir + '/global_morph_feature_ssd.png', format='png')
    #plt.show()
    plt.close()



# def plot_sample_size(bn_csv,outputDir,algorithms):
#     df_nd = pd.read_csv(bn_csv)
#
#     df_nd['algorithm'] = [algorithm_name_mapping[x] for x in df_nd['algorithm'] ]
#
#     plt.figure()
#
#     dfg = df_nd.groupby('algorithm')
#
#     sample_size_per_algorithm=[]
#     for alg in algorithms:
#         if alg in np.unique(df_nd['algorithm']):
#             sample_size_per_algorithm.append(dfg.get_group(alg).shape[0])
#         else:
#
#             sample_size_per_algorithm.append(0)
#
#
#     sb.barplot(y=range(algorithms.size),x=np.array(sample_size_per_algorithm),orient="h")
#     #sb.set_context("talk", font_scale=1.0)
#
#     algorithm_names = [algorithm_name_mapping[x] for x in algorithms]
#     plt.yticks(range(algorithms.size), algorithm_names)
#     plt.subplots_adjust(left=0.4,right=0.7)
#     plt.xlabel('Number of reconstructions')
#     plt.savefig(outputDir + '/valid_reconstruction_number.png', format='png')
#     #plt.show()
#     plt.close()
#     return
#


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
    sb.set_context("talk", font_scale=0.7)
    a=sb.barplot(y='algorithm', x='running_time', data=df_time,order=algorithms,orient="h")

    algorithm_names = [algorithm_name_mapping[x] for x in algorithms]

    a.set_yticklabels(['%s ($n$=%d )'%(algorithm_names[i], sample_size_per_algorithm[i]) for i in range(algorithms.size) ])
    plt.subplots_adjust(left=0.4,right=0.9, bottom=0.1, top=0.9)

    plt.xlabel('Running Time (seconds)')
    plt.savefig(outputDir + '/runningtime_n=logfiles.png', format='png')
    #plt.show()
    plt.close()

    return



def plot_running_time_validation(time_csv,neuron_distance_csv, outputDir, algorithms):

    df_nd = pd.read_csv(neuron_distance_csv)
    df_nd_group_by_image = df_nd.groupby('image_file_name')

    df_time = pd.read_csv(time_csv)

    df_time['running_time'] = df_time['running_time'] /1000.0


    dfg = df_time.groupby('image_file_name')
    all_images = np.unique(df_time['image_file_name'])

    df_out = pd.DataFrame()
    sample_size_per_algorithm=[]
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

             alg =  df_image.iloc[i]['algorithm']

             if alg in np.unique(df_nd_image['algorithm']):
                 id = df_time_filled_template[df_time_filled_template.algorithm == alg].index[0]
                 df_time_filled_template.ix[id] = df_image.iloc[i]

        df_out = df_out.append(df_time_filled_template,ignore_index=True)

    a=sb.barplot(y='algorithm', x='running_time', data=df_out, order=algorithms, orient="h")

    dfgt = df_nd.groupby('algorithm')
    sample_size_per_algorithm=[]
    m_algorithms = np.unique(df_time.algorithm)
    for alg in algorithms:
        #print alg
        if alg in m_algorithms:
            sample_size_per_algorithm.append(dfgt.get_group(alg).shape[0])
        else:
            sample_size_per_algorithm.append(0)

    algorithm_names = [algorithm_name_mapping[x] for x in algorithms]
    a.set_yticklabels(['%s ($n$=%d )'%(algorithm_names[i], sample_size_per_algorithm[i]) for i in range(algorithms.size) ])

    #a.set_yticklabels(algorithm_names)
    plt.subplots_adjust(left=0.4,right=0.9, bottom=0.1, top=0.9)
    plt.xlabel('Running Time (seconds): 1 hour wall time is used for missing/invalid reconstructions')
    plt.savefig(outputDir + '/runningtime_1hourForNA.png', format='png')
    #plt.show()
    plt.close()
    return



def plot_common_set():

    return

