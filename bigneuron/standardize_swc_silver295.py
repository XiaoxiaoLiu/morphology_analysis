__author__ = 'xiaoxiaol'
__author__ = 'xiaoxiaol'

# run standardize swc to make sure swc files have one single root, and sorted, and has the valid type id ( 1~4)




import matplotlib.pyplot as plt
import seaborn as sb
import os
import os.path as path
import numpy as np
import pandas as pd
import platform
import sys
import glob


if (platform.system() == "Linux"):
    WORK_PATH = "/local1/xiaoxiaol/work"
else:
    WORK_PATH = "/Users/xiaoxiaoliu/work"

p = WORK_PATH + '/src/morphology_analysis'
sys.path.append(p)


import blast_neuron.blast_neuron_comp as bn


### main
data_DIR = "/data/mat/xiaoxiaol/data/big_neuron/silver/20160113_merged_silver_gt"
output_dir = data_DIR

silver_dir = "/data/mat/xiaoxiaol/data/big_neuron/silver/silver_295_swcs"
SILVER_CSV = "/data/mat/xiaoxiaol/data/big_neuron/silver/silver_295_swcs/silver.csv"  #manually generated
df_silver= pd.read_csv(SILVER_CSV)
keys = df_silver['image_file_name']
values = df_silver['silver_swc_file']
map_image_to_silver_swc = dict(zip(keys, values))


lookup_image_id_table_file = data_DIR +"/../image_name_lookup_table.csv"
df_lookup_table = pd.read_csv(lookup_image_id_table_file)

os.system("rm "+data_DIR+"/qsub/*.qsub")
os.system("rm "+data_DIR+"/qsub/*.o*")

subdirs = [x[0] for x in os.walk(data_DIR)]

count = 0
gap_threshold =1
for recon_dir in subdirs[1:]:
        if 'auto_recons' in  recon_dir:
                  folder_name = recon_dir.split('/')[-1]
                  if '.v3d' in folder_name:
                      print folder_name
                  else:
                      continue

                  image_id = int(folder_name.split(".")[0])
                  image_file_name = df_lookup_table.image_file_name[image_id-1]


                  gs_swc_file = map_image_to_silver_swc[image_file_name]


                  out_gs_file = gs_swc_file.split('.swc')[0]+'.strict.swc'
                  if not os.path.exists(out_gs_file):
                      bn.standardize(input_swc_path=gs_swc_file, ref_swc_file=gs_swc_file,output_swc_path=out_gs_file, gap_threshold=1, new_type=3, GEN_QSUB = 1, qsub_script_dir=  data_DIR+"/qsub", id=None)

                  swc_files = glob.glob(recon_dir+'/*.swc')
                  for i in range(len(swc_files)) :
                      input = swc_files[i]
                      fn=(input.split('/')[-1]).split('.swc')[0]+'.strict.swc'
                      out_fn=recon_dir+"/../processed/"+fn
                      if not os.path.exists(out_fn):
                           if os.path.getsize(input) < 5*1024*1024:
                              bn.standardize(input_swc_path=input, ref_swc_file=gs_swc_file,output_swc_path=out_fn, gap_threshold=1, new_type=3, GEN_QSUB = 1, qsub_script_dir= data_DIR+"/qsub", id=None)
                           else:
                              print "swc file size too big, skip:",input

                  #bn.standardize(input_swc_path=recon_dir+'/../consensus3.eswc', ref_swc_file=gs_swc_file,output_swc_path=recon_dir+'/../consensus3.strict.swc', gap_threshold=1, new_type=3, GEN_QSUB = 1, qsub_script_dir=  data_DIR+"/qsub", id=None)


