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
import pandas as pd
import numpy as np
import time




def calculate_average_all_pair_distance(csv_file, hasConsensus=True):
    #df_in is the output csv from median_swc plugin
    #it contains unique pair-wise distances
    #consensus is usually the last input the median_swc() inputs, so it won't show up in the "swc_file_name1" column
    #output the average distances array

    #remove invalid results
    df_out=pd.DataFrame()
    if not os.path.exists(csv_file):
        return df_out

    df_f = pd.read_csv(csv_file)
    if df_f.empty:
        return df_out
    df_in = df_f[df_f['sum_distance'] >0]



    df_out = pd.DataFrame(columns = ['swc_file_name','average_sum_distance','average_structure_difference','average_max_distance'])

    dfg1 = df_in.groupby('swc_file_name1')
    dfg2 = df_in.groupby('swc_file_name2')

    swc_names = pd.unique(df_in['swc_file_name1'])
    swc_names_2 = pd.unique(df_in['swc_file_name2'])
    consensus_file_name = df_in['swc_file_name2'].tail(1).values[0]
    if 'consensus' not in consensus_file_name:
        #print  "missing consensus"
        return df_out


    row = 0
    for swc_name in swc_names:
        a = dfg1.get_group(swc_name)
        a = a[a['swc_file_name2']!=consensus_file_name]


        b = pd.DataFrame(columns = ['swc_file_name1','swc_file_name2','sum_distance','structure_difference','max_distance']) #empty
        if swc_name in swc_names_2:
            b = dfg2.get_group(swc_name)


        num_of_swcs = len(a) +len(b)
        df_out.loc[row,'swc_file_name']= swc_name.split('/')[-1]

        df_out.loc[row,'average_sum_distance'] = (a['sum_distance'].sum() + b['sum_distance'].sum())/ num_of_swcs
        df_out.loc[row,'average_structure_difference'] = a['structure_difference'].sum() + b['structure_difference'].sum()/num_of_swcs
        df_out.loc[row,'average_max_distance'] = a['max_distance'].sum() + b['max_distance'].sum()/num_of_swcs

        row = row +1


    df_out.loc[row,'swc_file_name']= consensus_file_name.split('/')[-1]
    consensus_group = dfg2.get_group(consensus_file_name)
    df_out.loc[row,'average_sum_distance'] = consensus_group['sum_distance'].sum() / (num_of_swcs+1)
    df_out.loc[row,'average_structure_difference'] = consensus_group['structure_difference'].sum() / (num_of_swcs+1)
    df_out.loc[row,'average_max_distance'] = consensus_group['max_distance'].sum() / (num_of_swcs+1)


    return df_out





def plot_compare_median_consensus(output_dir, df_order, metric, type = 'ts',DISPLAY = 0):
    plt.figure()


    if type =='ts':
        #sb.tsplot(data=df_order, value=metric,time='order',unit="algorithm",condition="algorithm",err_style="unit_traces")

        ax = sb.boxplot(x=metric, y="algorithm", data=df_order,
                 whis=np.inf, color="c")

        # Add in points to show each observation
        sb.stripplot(x=metric, y="algorithm", data=df_order,
                jitter=True, size=3, color=".3", linewidth=0)
        ax.set_xscale("log")
        sb.despine(trim=True)

       # plt.xlabel('images sorted by the average neuron distance of the median reconstruction')
        plt.savefig(output_dir + '/ts_compare_median_with_consensus_'+metric+'.png', format='png')


    if type =='lm':
        sb.lmplot(x="order", y=metric, hue="algorithm", data=df_order)
        plt.xlabel('images sorted by the average neuron distance of the median reconstruction')
        plt.savefig(output_dir + '/lm_compare_median_with_consensus_'+metric+'.lm.png', format='png')

    if DISPLAY:
         plt.show()
    plt.close()







###################################################################################################################################################
def pipe(input_data_dir, output_dir, imageIDs, distance_file_postfix='median_distances.csv',COLLECT_FROM_DISTANCE_MATRIX=1,EXTRACT_MEDIAN_CONSENSUS=1, DISPLAY=0):
    BeginTime = time.time()
    all_average_csv = output_dir+'/all_averaged_distances.csv'
    ######################################

    if COLLECT_FROM_DISTANCE_MATRIX:
         #remove empty files
        os.system('find '+ input_data_dir +' -size 0 -delete')
        os.system('mkdir '+ output_dir)

    else:
        print "Distance csv files had been collected. Skip collecting."


    if COLLECT_FROM_DISTANCE_MATRIX:
        df_all = pd.DataFrame(columns=['image_id', 'algorithm','swc_file_name','average_sum_distance','average_structure_difference','average_max_distance'])
        count = 0
        for image_id in imageIDs:

            df_image_filled_template = pd.DataFrame(columns = df_all.columns)

            df_image_filled_template['image_id'] = image_id

            csv_file = input_data_dir +'/'+image_id+'_'+distance_file_postfix

            if not os.path.exists(csv_file): # for gold 163, consensus results are stored in individual image folders
                 csv_file = input_data_dir +'/'+image_id+'/'+distance_file_postfix
                 if not os.path.exists(csv_file):
                     print  "missing "+csv_file
                     continue

            print "read "+csv_file
            df_ff = calculate_average_all_pair_distance(csv_file, hasConsensus = True)
            if not df_ff.empty:
               for i in range(df_ff.shape[0]):
                    if isinstance(df_ff.iloc[i].swc_file_name, basestring) :
                        alg = rp.matchFileToAlgorithmName((df_ff.iloc[i].swc_file_name).split('/')[-1])
                    else:
                        print "nan input swc_file_name"
                        print df_ff.iloc[i].swc_file_name
                        continue

                    df_image_filled_template.loc[count] = [image_id, alg,df_ff.iloc[i]['swc_file_name'],df_ff.iloc[i]['average_sum_distance'],df_ff.iloc[i]['average_structure_difference'],
                                                        df_ff.iloc[i]['average_max_distance']]
                    count = count +1
            else:
                print "empty df_ff:"+csv_file

            df_all = df_all.append(df_image_filled_template,ignore_index=True)

        df_all.to_csv(all_average_csv, index=False)

        print "Done collecting median distances"
        print "Output:" + all_average_csv



    #####################################################################
    PLOT_algorithm_consensus = 1
    metric = 'average_sum_distance'
    if PLOT_algorithm_consensus:
        df_all = pd.read_csv(all_average_csv)

        all_algorithms = np.unique(df_all.algorithm)
        plt.figure()
        sb.set_context("talk", font_scale=0.7)

        dfg = df_all.groupby('algorithm')

        sample_size_per_algorithm = np.zeros(all_algorithms.size)

        jj = 0
        for alg in all_algorithms:
            df_a = dfg.get_group(alg)
            df_a = df_a[df_a[metric]>=0]
            sample_size_per_algorithm[jj] = df_a.shape[0]
            jj = jj+1

        order = sample_size_per_algorithm.argsort()
        algorithms_ordered = all_algorithms[order[::-1]]
        sample_size_per_algorithm =sample_size_per_algorithm[order[::-1]]


        a = sb.barplot(y='algorithm', x=metric, data=df_all, order = algorithms_ordered)
        #a = sb.tsplot(data=df_all, time='image_id', value='total_average_distance')


        algorithm_names = [rp.map_better_algorithm_name(x) for x in algorithms_ordered]
        a.set_yticklabels(['%s ($n$=%d )'%(algorithm_names[i], sample_size_per_algorithm[i]) for i in range(algorithms_ordered.size) ])

        plt.subplots_adjust(left=0.4, bottom=0.1, top=0.9)
        plt.savefig(output_dir + '/compare_distance_plot.png', format='png')
        if DISPLAY:
             plt.show()
        plt.close()
        print "Done plotting algorithm comparison."
    #####################################################




    metric = 'average_sum_distance'

    output_median_and_consensus_csv_file = output_dir +'/extracted_median_consensus.csv'
    if not EXTRACT_MEDIAN_CONSENSUS :
        print "Median and consensus swc info had been collected. Skip this step."

    if EXTRACT_MEDIAN_CONSENSUS:
        df_all = pd.read_csv(all_average_csv)

        print df_all.shape
        dfg = df_all.groupby('image_id')
        df_median_and_consensus = pd.DataFrame(columns=['image_id', 'algorithm','swc_file_name','average_sum_distance','average_structure_difference','average_max_distance'])
        PLOT_imageIDs = pd.unique(df_all['image_id'])
        count = 0
        for image_id in PLOT_imageIDs:
            #print "image_id: "+ str( image_id)

            df_image = dfg.get_group(image_id)
            #drop nans
            df_image.dropna(axis=0, how="any", inplace =True)
            if len(df_image) <1:
                print "no valid recons found for :" + image_id
                continue

            i = 0
            for fn in df_image['swc_file_name']:

                if 'consensus' in fn:
                    break
                i= i+1
            if i>= len(df_image):
                print 'wrong consensus id found' + str(i)  +"mage id:" +image_id
                continue
            df_median_and_consensus.loc[count] =[image_id, 'consensus',df_image.iloc[i]['swc_file_name'],df_image.iloc[i]['average_sum_distance'],
                                                 df_image.iloc[i]['average_structure_difference'],df_image.iloc[i]['average_max_distance']]
            count = count +1


            df_image.drop(df_image.index[[i]], axis=0, inplace =True)

            df_image.sort(columns=['average_sum_distance', 'average_structure_difference','average_max_distance'], ascending = True,inplace=True)
            df_median_and_consensus.loc[count] =[image_id, 'median',df_image.iloc[0]['swc_file_name'],df_image.iloc[0]['average_sum_distance'],
                                                 df_image.iloc[0]['average_structure_difference'],df_image.iloc[0]['average_max_distance']]

            count = count +1


        df_median_and_consensus.to_csv(output_median_and_consensus_csv_file)
        print "Done extracting median distances"
        print "Output median and consensus distances for each image:"+output_median_and_consensus_csv_file



    PLOT_MEDIAN_CONSENSUS = 1
    if PLOT_MEDIAN_CONSENSUS:
          # reorer by distance
        df_median_and_consensus = pd.read_csv(output_median_and_consensus_csv_file)
        dfg = df_median_and_consensus.groupby('algorithm')
        df_consensus = dfg.get_group('consensus')
        df_median = dfg.get_group('median')

        df_median.reset_index(inplace=True)
        df_consensus.reset_index(inplace=True)


        #sort by average distance
        df_median.sort(columns=['average_sum_distance'], inplace=True)
        df_median['order'] = range(0,len(df_median))

        df_consensus = df_consensus.iloc[df_median.index]
        df_consensus['order'] = range(0,len(df_median))

        df_diff=pd.DataFrame(columns=['image_id','difference'])
        df_diff['difference'] = df_median['average_sum_distance'] -df_consensus['average_sum_distance']
        df_diff['image_id'] = df_median['image_id']

        df_ff =  df_diff['difference']
        print "median - consensus:"

        max_idx = np.nanargmin(np.abs(df_ff-df_ff.max()))
        print "max value: %f, image : %s"  % (df_ff.max(), df_median.iloc[max_idx]['image_id'])

        min_idx = np.nanargmin(np.abs(df_ff-df_ff.min()))
        print "min value: %f, image id: %s"  % (df_ff.min(),df_median.iloc[min_idx]['image_id'])

        median_idx = np.nanargmin(np.abs(df_ff-df_ff.median()))
        print "median value: %f, image : %s"  % (df_ff.median(),df_median.iloc[median_idx]['image_id'])


        df_ff_big = df_diff[df_diff['difference']>0]
        print "consensus is closer to each reconstructions than the median reconstructions in %.2f percent of the %d total images"  %( 100*float(len(df_ff_big))/len(df_diff), len(df_diff))
        df_ff_small = df_diff[df_diff['difference']<0]
        df_ff_small.to_csv(output_dir+'/investigate.csv')
        print "investigate the following cases list in:"+ output_dir+'/investigate.csv'

        #make sure the image_ids are matching
        for i in range(0,len(df_median)):
            if df_consensus.iloc[i]['image_id'] != df_median.iloc[i]['image_id']:
                print "error matching"
                print df_consensus.iloc[i]['image_id']
                print  df_median.iloc[i]['image_id']
                print exit()



        frames=[df_consensus,df_median]
        df_order = pd.concat(frames)
        for type in ['ts']:
            plot_compare_median_consensus(output_dir,df_order, 'average_sum_distance',type, DISPLAY)
            plot_compare_median_consensus(output_dir,df_order, 'average_max_distance',type, DISPLAY)
            plot_compare_median_consensus(output_dir,df_order, 'average_structure_difference',type, DISPLAY)


    endtime = time.time()
    dur=(endtime-BeginTime)/60
    print "All takes  %d min" %(dur)


