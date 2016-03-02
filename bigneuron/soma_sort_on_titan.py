__author__ = 'xiaoxiaoliu'



import pandas as pd
import os
import sys
import platform

if (platform.system() == "Linux"):
    WORK_PATH = "/local1/xiaoxiaol/work"
else:
    WORK_PATH = "/Users/xiaoxiaoliu/work"

p =  WORK_PATH + '/src/morphology_analysis'
sys.path.append(p)
import pandas as pd
import numpy as np
import os
import blast_neuron.blast_neuron_comp as bn


data_DIR = "/lustre/atlas2/nro101/proj-shared/BigNeuron/data/gold_163_all_soma_sort_s1"
output_dir = data_DIR


neuron_distance_csv = "/data/mat/xiaoxiaol/data/big_neuron/silver/20160113_merged_gold_gt/neuron_distances_with_gold.csv"
df_nd = pd.read_csv(neuron_distance_csv)


images = np.unique(df_nd['image_file_name'])


dfg = df_nd.groupby('image_file_name')

count=0
for im in images:

     df_image = dfg.get_group(im)

     df_image=df_image.sort(['neuron_distance'])


     tmp = df_image.iloc[0]['swc_file']


     im_id = tmp.split('/')[-2]  # 2.v3dpbd

     out_dir = output_dir  + '/' + im_id.split('.')[0]


     gold_swc =  df_image.iloc[0]['gold_swc_file']

     output_gold_swc = out_dir+'/00_'+gold_swc.split('/')[-1]


     i=1

     for swc_file in df_image['swc_file']:
          string=str(i)
          if i < 10:
                string = '0'+str(i)
          out_swc = out_dir +'/auto_recons/' + string +'_'+ swc_file.split('/')[-1]

          soma_sorted_swc = out_dir +'/processed/'+ string +'_'+ swc_file.split('/')[-1]

          bn.soma_sorting(output_gold_swc, inputswc_path = out_swc, outputswc_path = soma_sorted_swc, step_size = 1 ,
                               logfile=soma_sorted_swc+'.log', GEN_QSUB = 2, qsub_script_dir= "./txt_jobs", id= count)

          i=i+1
          count= count+1
