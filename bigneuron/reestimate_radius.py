__author__ = 'xiaoxiaol'
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
import  blast_neuron.blast_neuron_comp as bn


### main
data_DIR = "/data/mat/xiaoxiaol/data/big_neuron/silver/0401_gold163_all_soma_sort"
output_dir = data_DIR
#run_consensus(data_DIR, output_dir)


os.system("rm "+data_DIR+"/qsub2/*.qsub")
os.system("rm "+data_DIR+"/qsub2/*.o*")


for item in os.listdir(data_DIR):
    folder_name = os.path.join(data_DIR, item)

    if os.path.isdir(folder_name):
        print folder_name
        imagefile = glob.glob(folder_name+'/*.v3dpbd')
        imagefile.extend(glob.glob(folder_name+'/*.v3draw'))
        files =glob.glob(folder_name+'/*.strict.swc')

        if len(files)>0 and len(imagefile)>0:
            gs_swc_file =files[0]
            if not os.path.exists(gs_swc_file+".out.swc"):

                bn.estimate_radius(input_image=imagefile[0], input_swc_path=gs_swc_file,bg_th=40, GEN_QSUB = 0, qsub_script_dir= output_dir+"/qsub2", id=None)


