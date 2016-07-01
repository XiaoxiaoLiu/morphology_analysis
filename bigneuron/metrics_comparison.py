__author__ = 'xiaoxiaol'


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb
import os
import os.path as path
import numpy as np
import recon_prescreening as rp



def plot_nblast_distance(neuron_distance_csv, outputDir,algorithms=None):

    df_nd = pd.read_csv(neuron_distance_csv)

    if not path.exists(outputDir):
        os.mkdir(outputDir)


    if algorithms is None:
        algorithms = np.unique(df_nd.algorithm)
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
        algorithms = np.unique(df_nd.algorithm)
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
        algorithms = np.unique(df_nd.algorithm)
    ### all algorithm plot
    dfg = df_nd.groupby('algorithm')
    sample_size_per_algorithm=[]
    for alg in algorithms:
        sample_size_per_algorithm.append(dfg.get_group(alg).shape[0])



    #plot the average node distances
    plt.figure()
    sb.set_style("white")

    g = sb.lmplot(x="image_id", y="diadem_score", hue="algorithm", data=df_nd,fit_reg=False)


    #a=sb.barplot(y='algorithm', x=['diadem_score','nblast_bi_score'], data=df_nd,order=algorithms)
   # algorithm_name_mapping = rp.get_algorithm_name_dict()
   # algorithm_names = [algorithm_name_mapping[x] for x in algorithms]
    #a.set_yticklabels(['%s ($n$=%d )'%(algorithm_names[i], sample_size_per_algorithm[i]) for i in range(algorithms.size) ])
    #sb.set_context("talk", font_scale=3.0)
    #plt.xticks(rotation="90")
    #plt.xlabel('All Scores')
    #plt.subplots_adjust(left=0.5,right=0.95,bottom=0.1, top=0.9)
    plt.savefig(outputDir + '/All_scores.png', format='png')
    plt.show()
   # plt.close()


    return




data_DIR="/data/mat/xiaoxiaol/data/big_neuron/BTU"

# df_diadem = pd.read_csv(data_DIR+"/diadem/Result_strict.csv")
# df_nblast = pd.read_csv(data_DIR+"/nblast/nblast_score_after_postprocessing.csv")#
# df_diadem = rp.parse_and_add_algorithm_info(df_diadem)
# df_diadem.to_csv(data_DIR+"/diadem/diadem_score_with_meta.csv",index=False)
# df_nblast = rp.parse_and_add_algorithm_info(df_nblast)
# df_nblast.to_csv(data_DIR+"/nblast/nblast_score_with_meta.csv",index=False)

#df_all= pd.merge(df_diadem, df_nblast, on=['swc_file_name','algorithm','image_id'])
#df_all.to_csv(data_DIR+"/all_metrics.csv",index=False)

# df_diadem.dropna()
# df_nblast.dropna()
# df_diadem = pd.read_csv(data_DIR+"/diadem/diadem_score_with_meta.csv")
# df_nblast = pd.read_csv(data_DIR+"/nblast/nblast_score_with_meta.csv")



df_all= pd.read_csv(data_DIR+"/all_metrics.csv")
algorithms = np.unique(df_all['algorithm'])



#plot_nblast_distance(data_DIR+"/nblast/nblast_score_with_meta.csv", data_DIR+"/nblast",algorithms)
#plot_diadem_score(data_DIR+"/diadem/diadem_score_with_meta.csv", data_DIR+"/diadem",algorithms)


plot_all_score(data_DIR+"/all_metrics.csv", data_DIR,algorithms)

