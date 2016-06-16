__author__ = 'xiaoxiaol'
__author__ = 'xiaoxiaol'

import sys
import os
import platform

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
import glob

## merge data for plotting


new_data_DIR ="/data/mat/xiaoxiaol/data/big_neuron/silver/reconstructions_20160114_LCMBoost_v2"

#destination folder
data_DIR ="/data/mat/xiaoxiaol/data/big_neuron/silver/20160113_merged_silver_gt"


lookup_image_id_table_file = "/data/mat/xiaoxiaol/data/big_neuron/silver/image_name_lookup_table.csv"

#postfix = "_updated"
postfix=""
#rename the files so that it won't overwrite the previous results



############### merge files to the destination
swcfiles = glob.glob(os.path.join(new_data_DIR+'/*/', '*.swc'))
print "contains ",len(swcfiles), "  swc files"

for m_file in swcfiles:
    fn = m_file.split('/')[-1]
    new_file = data_DIR +'/auto_recons/'+ m_file.split('/')[-2]+'/'+fn.split('.swc')[0]+postfix +".swc"
    #print new_file
    os.system("cp "+ m_file + " "+ new_file)




################ reformat running time
time_csv = new_data_DIR +"/running_time.csv"

df_rt = pd.DataFrame()
plugins=[]
times=[]
images=[]
with open(time_csv) as f:
    for line in f:
        logfile = line.split(' ')[0]
        time = line.split(' ')[-1]

        times.append(time[:-1])  # chop \n
        tmp =logfile.split('_time')[0]
        if 'v3dpbd' in tmp:
            plugins.append(tmp.split('v3dpbd_')[-1] +postfix)
            images.append(tmp.split('v3dpbd_')[0] +'v3dpbd')
        else:
            plugins.append(tmp.split('v3draw_')[-1]+postfix)
            images.append(tmp.split('v3draw_')[0] +'v3raw')

df_rt['image']=images
df_rt['plugin'] = plugins
df_rt['running_time'] = times

df_rt.to_csv(new_data_DIR+'/running_time_reorg.csv', index=False)

################
#merge time spreasheets


os.system('mv '+ data_DIR+'/auto_recons/running_time_merged.csv ' + data_DIR+'/auto_recons/running_time.csv')

df_old_rt = pd.read_csv(data_DIR +"/auto_recons/running_time.csv")
df_new_rt = pd.concat([df_old_rt,df_rt])
df_new_rt.to_csv(data_DIR+'/auto_recons/running_time_merged.csv', index=False)





#running_time
#vim edit merge running_time.csv into one
