__author__ = 'xiaoxiaol'

# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 17:58:24 2015

@author: xiaoxiaoliu
"""



import pandas as pd
import platform
import os

def read_swc(infile):
    ###n,type,x,y,z,radius,parent
     swc_df = pd.read_csv(infile, sep=" ",skiprows = 3, header =None, names =['n','type','x','y','z','radius','parent'] )
     return swc_df



def write_swc(df_out,outfile):
     ###n,type,x,y,z,radius,parent
     df_out.to_csv(outfile, index=False, header= False, sep=" ")




# 0 = undefined
# 1 = soma
# 2 = axon
# 3 = dendrite
# 4 = apical dendrite
# 5 = fork point
# 6 = end point
# 7 = custom
def remove_swc_by_id(type_id, swc_file, output_swc_file):
      if (output_swc_file == None):
          output_swc_file = swc_file[0:-4] +'_'+ str(type_id) +'.swc'
      swc_df = read_swc(swc_file)
      extracted = swc_df[ swc_df['type'] != type_id ]


      write_swc(extracted,output_swc_file)


if (platform.system() == "Linux"):
    WORK_PATH = "/local1/xiaoxiaol/work"
else:
    WORK_PATH = "/Users/xiaoxiaoliu/work"

data_DIR = WORK_PATH + "/data/lims2/0923_pw_aligned"
df = pd.read_csv(data_DIR +'/preprocessed-with_axon/list.csv')
output_dir = data_DIR +'/axon_removed'
if not os.path.exists(output_dir):
        os.mkdir(output_dir)

for i in range(df.shape[0]):
     swc_file = df['swc_file'][i]
     output_swc_file = output_dir+'/'+swc_file[0:-4]+'.swc'
     remove_swc_by_id(2, data_DIR+'/preprocessed-with_axon/'+swc_file, output_swc_file)
