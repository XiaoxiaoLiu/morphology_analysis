import pandas as pd
import numpy as np
import os



data_DIR = "/data/mat/xiaoxiaol/data/big_neuron/silver"

output_dir = "/data/mat/xiaoxiaol/data/big_neuron/silver/ready_for_tile_vis"

neuron_distance_csv = "/data/mat/xiaoxiaol/data/big_neuron/silver/20160113_merged_gold_gt/neuron_distances_with_gold_filtered.csv"

num_of_selected_swc = 14




df_image_location =  pd.read_csv('/data/mat/xiaoxiaol/data/Hanchuan_curated/image_file_location_checkup.csv')
keys = df_image_location['image_file_name']
values = df_image_location['file_path']
image_checkup = dict(zip(keys, values))



df_nd = pd.read_csv(neuron_distance_csv)


#merge with the gold79 subset
df_79 = pd.read_csv('/Users/xiaoxiaoliu/work/data/gold79/gold.csv')
images = np.unique(df_79['image_file_name'])



dfg = df_nd.groupby('image_file_name')

for im in images:

     df_image = dfg.get_group(im)


     #print df_image['swc_file']
     #sort by distance
     df_image.sort_values(['neuron_distance'], ascending=[1], inplace=True)
     #print df_image['swc_file']

     tmp = df_image.iloc[0]['swc_file']


     im_id = tmp.split('/')[-2]  # 2.v3dpbd


     out_dir = output_dir  + '/' + im_id
     if  not os.path.exists(out_dir):
        os.mkdir(out_dir)

     gold_swc =  df_image.iloc[0]['gold_swc_file']


     image_file =  image_checkup[im]
     print image_file

     output_swc = out_dir+'/00_'+gold_swc.split('/')[-1]
     os.system("cp "+gold_swc + " "+ output_swc)

     output_image = out_dir +'/'+im
     #copy image
     #os.system("cp "+image_file + " "+ output_image)

     i=1
     for swc_file in df_image['swc_file'][0:num_of_selected_swc]:
          string=str(i)
          if i < 10:
                string = '0'+str(i)
          out_swc = out_dir +'/' + string +'_'+ swc_file.split('/')[-1]
          os.system("cp "+ swc_file + " "+ out_swc)
          i=i+1








