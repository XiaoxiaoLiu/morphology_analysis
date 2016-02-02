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


data_DIR = "/data/mat/xiaoxiaol/data/big_neuron/silver"

output_dir = "/data/mat/xiaoxiaol/data/big_neuron/silver/gold_163_all_somasorted"
os.system("mkdir "+output_dir)

neuron_distance_csv = "/data/mat/xiaoxiaol/data/big_neuron/silver/20160113_merged_gold_gt/neuron_distances_with_gold.csv"

#num_of_selected_swc = 14




df_image_location =  pd.read_csv('/data/mat/xiaoxiaol/data/Hanchuan_curated/image_file_location_checkup.csv')
keys = df_image_location['image_file_name']
values = df_image_location['file_path']
image_checkup = dict(zip(keys, values))


df_nd = pd.read_csv(neuron_distance_csv)


#merge with the gold79 subset
#df_79 = pd.read_csv('/Users/xiaoxiaoliu/work/data/gold79/gold.csv')
#images = np.unique(df_79['image_file_name'])
#print images.size

images = np.unique(df_nd['image_file_name'])


dfg = df_nd.groupby('image_file_name')

#df_feature_79=pd.DataFrame()
for im in images:

     df_image = dfg.get_group(im)
#     df_feature_79=df_feature_79.append(df_image,ignore_index=True)

     #print df_image['swc_file']
     #sort by distance
     df_image.sort_values(['neuron_distance'], ascending=[1], inplace=True)
     #print df_image['swc_file']

     tmp = df_image.iloc[0]['swc_file']


     im_id = tmp.split('/')[-2]  # 2.v3dpbd
     # ano_file= im_id+".recons.ano"
     # median_log = im_id+".recons.ano.median.log
     # median_fn = bn.read_median_swc_log(ano_file, median_log)
     # print median_fn


     out_dir = output_dir  + '/' + im_id.split('.')[0]


     if  not os.path.exists(out_dir):
        os.mkdir(out_dir)

     gold_swc =  df_image.iloc[0]['gold_swc_file']


     image_file =  image_checkup[im]
     #print image_file

     output_swc = out_dir+'/00_'+gold_swc.split('/')[-1]

     os.system("cp "+gold_swc + " "+ output_swc)

     output_image = out_dir +'/'+im
     #copy image
     #os.system("cp "+image_file + " "+ output_image)

     i=1
     os.system('mkdir '+ out_dir+'/auto_recons')
     os.system('mkdir '+ out_dir+'/processed')
     for swc_file in df_image['swc_file']:
          string=str(i)
          if i < 10:
                string = '0'+str(i)
          out_swc = out_dir +'/auto_recons/' + string +'_'+ swc_file.split('/')[-1]
          os.system("cp "+ swc_file + " "+ out_swc)
          soma_sorted_swc = out_dir +'/processed/'+ string +'_'+ swc_file.split('/')[-1]
          if  (not os.path.exists(soma_sorted_swc) ) and  os.path.getsize(out_swc) < 1000000:
              bn.soma_sorting(gold_swc, inputswc_path = out_swc, outputswc_path = soma_sorted_swc, step_size = 3 )
          i=i+1
     bn.genLinkerFile( out_dir+'/auto_recons', out_dir+"/auto_recons/"+im_id+'.ano')
     bn.genLinkerFile( out_dir+'/processed', out_dir+"/processed/"+im_id+'.ano')

     bn.consensus(input_ano_path= out_dir+"/processed/"+im_id+'.ano', output_eswc_path=out_dir+"/processed/consensus_p2.eswc", method=2, GEN_QSUB = 0, qsub_script_dir= ".")


#df_feature_79.to_csv(data_DIR+"/gold_trainning_subset/neuron_distances.csv")
#print df_feature_79.algorithm
#print df_feature_79.image_file_name





