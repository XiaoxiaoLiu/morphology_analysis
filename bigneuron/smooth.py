__author__ = 'xiaoxiaol'
__author__ = 'xiaoxiaol'
__author__ = 'xiaoxiaol'
import sys
import os
import platform
import matplotlib.pyplot as plt
import seaborn as sb

if (platform.system() == "Linux"):
    WORK_PATH = "/local1/xiaoxiaol/work"
else:
    WORK_PATH = "/Users/xiaoxiaoliu/work"

p = WORK_PATH + '/src/morphology_analysis'
sys.path.append(p)

import  bigneuron.recon_prescreening as rp
import  bigneuron.plot_distances as plt_dist
import pandas as pd
import numpy as np





data_DIR="/home/xiaoxiaol/work/data/consensus_0225_anisosmooth"

#taiwan dataset
#00001 ~16266

imageIDs=[]
TOTAL_NUM=16227
for n in range(1,TOTAL_NUM):
    imageIDs.append(str(n).rjust(5,'0'))


LOOKUP_TABLE_FILE = "/data/mat/xiaoxiaol/data/big_neuron/silver/ported_neuron_tracing_spreadsheet.csv"
df_algorithms = pd.read_csv(LOOKUP_TABLE_FILE)

all_algorithms = pd.unique(df_algorithms.algorithm)



COLLECT=1
if COLLECT:

    df_all = pd.DataFrame(columns=['image_id', 'algorithm','swc_file_name','total_average_distance','total_structure_difference','total_max_distance'])

    for image_id in imageIDs:

        df_image_filled_template = pd.DataFrame(columns = df_all.columns)
        df_image_filled_template['algorithm'] = all_algorithms
        df_image_filled_template['image_id'] = image_id

        df_image_filled_template['swc_file_name'] = np.nan
        df_image_filled_template['total_average_distance'] = np.nan
        df_image_filled_template['total_structure_difference'] = np.nan
        df_image_filled_template['total_max_distance'] = np.nan


        csv_file = data_DIR + '/'+image_id+'_median_distances.csv'
        if os.path.exists(csv_file):
           df_f = pd.read_csv(csv_file)

           for i in range(df_f.shape[0]):
                if type(df_f.iloc[i].swc_file_name ) == str:
                    alg = rp.matchFileToAlgorithmName((df_f.iloc[i].swc_file_name).split('/')[-1])
                else:
                    print "nan input swc_file_name"
                    print df_f.iloc[i].swc_file_name
                    continue
                df_a=df_image_filled_template[df_image_filled_template.algorithm == alg]
                if df_a.shape[0] >0:
                  id = df_a.index[0]

                  df_image_filled_template.loc[id,'swc_file_name'] =df_f.iloc[i]['swc_file_name']
                  df_image_filled_template.loc[id,'total_average_distance'] =df_f.iloc[i]['total_average_distance']
                  df_image_filled_template.loc[id,'total_structure_difference'] = df_f.iloc[i]['total_structure_difference']
                  df_image_filled_template.loc[id,'total_max_distance'] = df_f.iloc[i]['total_max_distance']
                else:
                    print alg
                    print "no match!"
                    print df_f.iloc[i].swc_file_name


        df_all = df_all.append(df_image_filled_template,ignore_index=True)

    df_all.to_csv(data_DIR+'/../consensus_0225_anisosmooth_all_median_distances.csv', index=False)

PLOT=1
if PLOT:
    df_all = pd.read_csv(data_DIR+'/../consensus_0225_anisosmooth_all_median_distances.csv')
    plt.figure()
    sb.set_context("talk", font_scale=0.7)

    dfg = df_all.groupby('algorithm')

    sample_size_per_algorithm = np.zeros(all_algorithms.size)

    jj = 0
    for alg in all_algorithms:
        df_a = dfg.get_group(alg)
        df_a = df_a[df_a['total_average_distance']>=0]
        sample_size_per_algorithm[jj] = df_a.shape[0]
        jj = jj+1

    order = sample_size_per_algorithm.argsort()
    algorithms_ordered = all_algorithms[order[::-1]]
    sample_size_per_algorithm =sample_size_per_algorithm[order[::-1]]


    a = sb.barplot(y='algorithm', x='total_average_distance', data=df_all, order = algorithms_ordered)
    #a = sb.tsplot(data=df_all, time='image_id', value='total_average_distance')




    algorithm_names = [rp.map_better_algorithm_name(x) for x in algorithms_ordered]
    a.set_yticklabels(['%s ($n$=%d )'%(algorithm_names[i], sample_size_per_algorithm[i]) for i in range(algorithms_ordered.size) ])

    plt.subplots_adjust(left=0.4, bottom=0.1, top=0.9)
    plt.savefig(data_DIR + '/../test2.png', format='png')
    plt.show()
    plt.close()

