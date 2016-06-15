__author__ = 'xiaoxiaol'

# run standardize swc to make sure swc files have one single root, and sorted, and has the valide type id ( 1~4)




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

import  bigneuron.recon_prescreening as rp
import  bigneuron.plot_distances as plt_dist
import blast_neuron.blast_neuron_comp as bn


### main
data_DIR = "/data/mat/xiaoxiaol/data/big_neuron/silver/0401_gold163_all_soma_sort"
output_dir = data_DIR
#run_consensus(data_DIR, output_dir)



subdirs = [x[0] for x in os.walk(data_DIR)]
count = 0
for recon_dir in subdirs[1:]:
        folder_name = recon_dir.split('/')[-1]
        if 'auto_recons' in  folder_name:
                  print recon_dir
                  folder_name = recon_dir.split('/')[-1]
                  files = glob.glob(recon_dir+'/../00_*.swc')
                  if len(files)>0:
                       gs_swc_file =files[0]
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



