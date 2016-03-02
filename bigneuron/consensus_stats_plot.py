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





data_DIR="/home/xiaoxiaol/work/data/consensus_0225"


#taiwan dataset
#00001 ~16266

imageIDs=[]
for n in range(1,16227):
    imageIDs.append(str(n).rjust(5,'0'))


LOOKUP_TABLE_FILE = "/data/mat/xiaoxiaol/data/big_neuron/silver/ported_neuron_tracing_spreadsheet.csv"
df_algorithms = pd.read_csv(LOOKUP_TABLE_FILE)

all_algorithms = pd.unique(df_algorithms.algorithm)


def calculate_average(df_in, hasConsensus=True):
    #df_in is the output csv from median_swc plugin
    #it contains unique pair-wise distances
    #consensus is usually the last input the median_swc() inputs, so it won't show up in the "swc_file_name1" column
    #output the average distances array

    #remove invalid results
    df_in = df_in[df_in[' average_distance'] >0]



    df_out = pd.DataFrame(columns = ['swc_file_name','average_distance','average_structure_difference','average_max_distance'])
    dfg1 = df_in.groupby('swc_file_name1')
    dfg2 = df_in.groupby('swc_file_name2')

    swc_names = pd.unique(df_in['swc_file_name1'])
    swc_names_2 = pd.unique(df_in['swc_file_name2'])
    consensus_file_name = df_in['swc_file_name2'].tail(1).values[0]
    print consensus_file_name


    row = 0
    for swc_name in swc_names:
        a = dfg1.get_group(swc_name)
        a = a[a['swc_file_name2']!=consensus_file_name]

        b = pd.DataFrame(columns = ['swc_file_name1','swc_file_name2',' average_distance','structure_difference','max_distance']) #empty
        if swc_name in swc_names_2:
            b = dfg2.get_group(swc_name)


        num_of_swcs = len(a) +len(b)
        df_out.loc[row,'swc_file_name']= swc_name
        df_out.loc[row,'average_distance'] = (a[' average_distance'].sum() + b[' average_distance'].sum())/ num_of_swcs
        df_out.loc[row,'average_structure_difference'] = a['structure_difference'].sum() + b['structure_difference'].sum()/num_of_swcs
        df_out.loc[row,'average_max_distance'] = a['max_distance'].sum() + b['max_distance'].sum()/num_of_swcs
        row = row +1


    df_out.loc[row,'swc_file_name']= consensus_file_name
    consensus_group = dfg2.get_group(consensus_file_name)
    df_out.loc[row,'average_distance'] = consensus_group[' average_distance'].sum() / (num_of_swcs+1)
    df_out.loc[row,'average_structure_difference'] = consensus_group['structure_difference'].sum() / (num_of_swcs+1)
    df_out.loc[row,'average_max_distance'] = consensus_group['max_distance'].sum() / (num_of_swcs+1)

    return df_out




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
           df_ff = calculate_average_all_pair_distance(df_f, hasConsensus = True)

           for i in range(df_f.shape[0]):
                if type(df_f.iloc[i].swc_file_name ) == str:
                    alg = rp.matchFileToAlgorithmName((df_ff.iloc[i].swc_file_name).split('/')[-1])
                else:
                    print "nan input swc_file_name"
                    print df_f.iloc[i].swc_file_name
                    continue
                df_a=df_image_filled_template[df_image_filled_template.algorithm == alg]
                if df_a.shape[0] >0:
                  id = df_a.index[0]

                  df_image_filled_template.loc[id,'swc_file_name'] =df_ff.iloc[i]['swc_file_name']
                  df_image_filled_template.loc[id,'average_distance'] =df_ff.iloc[i]['average_distance']
                  df_image_filled_template.loc[id,'average_structure_difference'] = df_ff.iloc[i]['average_structure_difference']
                  df_image_filled_template.loc[id,'average_max_distance'] = df_ff.iloc[i]['average_max_distance']
                else:
                    print alg
                    print "no match!"
                    print df_f.iloc[i].swc_file_name


        df_all = df_all.append(df_image_filled_template,ignore_index=True)

    df_all.to_csv(data_DIR+'/../consensus_0225_all_median_distances.csv', index=False)

PLOT=1
if PLOT:
    df_all = pd.read_csv(data_DIR+'/../consensus_0225_all_median_distances.csv')
    plt.figure()
    sb.set_context("talk", font_scale=0.7)

    dfg = df_all.groupby('algorithm')

    sample_size_per_algorithm = np.zeros(all_algorithms.size)

    jj = 0
    for alg in all_algorithms:
        df_a = dfg.get_group(alg)
        df_a = df_a[df_a['total_average_distance']>=0]
        print df_a.shape
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
    plt.savefig(data_DIR + '/../test.png', format='png')
    plt.show()
    plt.close()

