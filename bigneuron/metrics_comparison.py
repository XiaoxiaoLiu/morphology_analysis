__author__ = 'xiaoxiaol'


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb
import os
import os.path as path
import numpy as np
import recon_prescreening as rp







def calculate_similarities(neuron_distance_csv,metric='neuron_distance', output_similarity_csv =None):
    df_nd = pd.read_csv(neuron_distance_csv)
    all_images = np.unique(df_nd.image_id)
    all_algorithms = np.unique(df_nd.algorithm)

    print "\n\nCalculate similarity based on " +metric
    print neuron_distance_csv + " has :"
    print str(all_algorithms.size) + " algorithms"
    print str(all_images.size) +" images"
    #print all_algorithms

    dfg = df_nd.groupby('image_id')

    df_out = pd.DataFrame()
    #sample_size_per_algorithm=[]
    for image in all_images:
        df_image = dfg.get_group(image)
        #sample_size_per_image.append(df_image.shape[0])

        # similarity == nan:  metric reports nan
        # similarity = 0 : missing entry ( missing recons)
        #df_image['similarity'] = np.exp(-(df_image[metric] - df_nd[metric].min()+0.000000001)/(df_nd[metric].max()-df_nd[metric].min()+0.000000001))
        df_image['similarity'] = (df_image[metric] - df_image[metric].min()+0.000000001)/(df_nd[metric].max()-df_nd[metric].min()+0.000000001)


        # construct a complete table, and fill the valid results
        df_image_filled_template = pd.DataFrame(columns = df_image.columns)
        df_image_filled_template.algorithm = all_algorithms
        df_image_filled_template.image_id = image
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







def plot_similarities(neuron_distance_csv, outputDir,algorithms=None,metric='neuron_distance',CASE_BY_CASE_PLOT = 0, value_label=None):
    df_nd_ori = pd.read_csv(neuron_distance_csv)
    algorithm_name_mapping = rp.get_algorithm_name_dict()

    df_nd = calculate_similarities(neuron_distance_csv,metric,output_similarity_csv=neuron_distance_csv+".similarity.csv")
    all_images = np.unique(df_nd.image_id)
    if not path.exists(outputDir):
        os.mkdir(outputDir)



    if algorithms is None:
        algorithms= order_algorithms_by_size(df_nd_ori)

    if CASE_BY_CASE_PLOT:
        dfg = df_nd.groupby('image_id')

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


def plot_nblast_distance(neuron_distance_csv, outputDir,algorithms=None):

    df_nd = pd.read_csv(neuron_distance_csv)

    if not path.exists(outputDir):
        os.mkdir(outputDir)


    if algorithms is None:
        algorithms= order_algorithms_by_size(df_nd)
    ### all algorithm plot
    dfg = df_nd.groupby('algorithm')
    sample_size_per_algorithm=[]
    for alg in algorithms:
        sample_size_per_algorithm.append(dfg.get_group(alg).shape[0])

    #plot the average node distances
    plt.figure()
    sb.set_context("talk", font_scale=0.7)
    a=sb.barplot(y='algorithm', x='nblast_bi_score', data=df_nd,order=algorithms)
    algorithm_name_mapping = rp.get_algorithm_name_dict()
    algorithm_names = [algorithm_name_mapping[x] for x in algorithms]
    a.set_yticklabels(['%s ($n$=%d )'%(algorithm_names[i], sample_size_per_algorithm[i]) for i in range(algorithms.size) ])
    #sb.set_context("talk", font_scale=3.0)
    #plt.xticks(rotation="90")
    plt.xlabel('Nblast Scores')
    plt.subplots_adjust(left=0.6,right=0.95,bottom=0.1, top=0.9)
    plt.savefig(outputDir + '/NBlast_Score.png', format='png')
    #plt.show()
    plt.close()

    return


def plot_diadem_score(neuron_distance_csv, outputDir,algorithms=None):

    df_nd = pd.read_csv(neuron_distance_csv)

    if not path.exists(outputDir):
        os.mkdir(outputDir)


    if algorithms is None:
       algorithms= order_algorithms_by_size(df_nd)
    ### all algorithm plot
    dfg = df_nd.groupby('algorithm')
    sample_size_per_algorithm=[]
    for alg in algorithms:
        sample_size_per_algorithm.append(dfg.get_group(alg).shape[0])



    #plot the average node distances
    plt.figure()
    sb.set_context("talk", font_scale=0.7)
    a=sb.barplot(y='algorithm', x='diadem_score', data=df_nd,order=algorithms)
    algorithm_name_mapping = rp.get_algorithm_name_dict()
    algorithm_names = [algorithm_name_mapping[x] for x in algorithms]
    a.set_yticklabels(['%s ($n$=%d )'%(algorithm_names[i], sample_size_per_algorithm[i]) for i in range(algorithms.size) ])
    #sb.set_context("talk", font_scale=3.0)
    #plt.xticks(rotation="90")
    plt.xlabel('Diadem Scores')
    plt.subplots_adjust(left=0.5,right=0.95,bottom=0.1, top=0.9)
    plt.savefig(outputDir + '/Diadem_score.png', format='png')
    #plt.show()
    plt.close()


    return







def plot_all_score(neuron_distance_csv, outputDir,algorithms=None):


    df_nd = pd.read_csv(neuron_distance_csv)


    if not path.exists(outputDir):
        os.mkdir(outputDir)


    if algorithms is None:
        algorithms= order_algorithms_by_size(df_nd)
    ### all algorithm plot
    dfg = df_nd.groupby('algorithm')
    sample_size_per_algorithm=[]
    for alg in algorithms:
        sample_size_per_algorithm.append(dfg.get_group(alg).shape[0])



    #plot the average node distances
    plt.figure()
    sb.set_style("white")

    #g = sb.lmplot(x="image_id", y="diadem_score", hue="algorithm", data=df_nd,fit_reg=False)


    a=sb.barplot(y='algorithm', x=['diadem_score','nblast_bi_score'], data=df_nd,order=algorithms)
    algorithm_name_mapping = rp.get_algorithm_name_dict()
    algorithm_names = [algorithm_name_mapping[x] for x in algorithms]
    a.set_yticklabels(['%s ($n$=%d )'%(algorithm_names[i], sample_size_per_algorithm[i]) for i in range(algorithms.size) ])
    sb.set_context("talk", font_scale=3.0)
    plt.xticks(rotation="90")
    plt.xlabel('All Scores')
    plt.subplots_adjust(left=0.5,right=0.95,bottom=0.1, top=0.9)
    plt.savefig(outputDir + '/All_scores.png', format='png')
    plt.show()
   # plt.close()

    return

def order_algorithms_by_size(df_data):

    algorithms = np.unique(df_data.algorithm)
    dfg = df_data.groupby('algorithm')
    sample_size_per_algorithm=[]
    for alg in algorithms:
         sample_size_per_algorithm.append(dfg.get_group(alg).shape[0])

    sorted_idex = np.argsort(np.array(sample_size_per_algorithm))

    algorithms=algorithms[sorted_idex]
    return algorithms



data_DIR="/data/mat/xiaoxiaol/data/big_neuron/BTU"

#output
result_csv_file = data_DIR+"/all_metrics.csv"
nblast_csv_file =data_DIR+"/nblast/nblast_score_with_meta.csv"
diadem_csv_file =data_DIR+"/diadem/diadem_score_with_meta.csv"

df_diadem = pd.read_csv(data_DIR+"/diadem/result_combined.csv")
df_diadem.dropna(inplace=True)
df_diadem = rp.parse_and_add_algorithm_info(df_diadem)
df_diadem.to_csv(diadem_csv_file,index=False)

df_nblast = pd.read_csv(data_DIR+"/nblast/nblast_score_after_postprocessing_6_29.csv")#
df_nblast.dropna(inplace=True)
df_nblast = rp.parse_and_add_algorithm_info(df_nblast)
df_nblast.to_csv(nblast_csv_file,index=False)

df_all= pd.merge(df_diadem, df_nblast, on=['swc_file_name','algorithm','image_id'])
df_all.dropna(inplace=True)

df_all.to_csv(result_csv_file,index=False)




plot_nblast_distance(nblast_csv_file, data_DIR+"/nblast")
plot_similarities(nblast_csv_file, data_DIR+"/nblast",metric='nblast_bi_score')

plot_diadem_score(diadem_csv_file, data_DIR+"/diadem")
plot_similarities(diadem_csv_file, data_DIR+"/diadem",metric='diadem_score')

#plot_all_score(data_DIR+"/all_metrics.csv", data_DIR)

