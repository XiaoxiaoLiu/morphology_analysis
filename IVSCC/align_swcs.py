__author__ = 'xiaoxiaoliu'

import os
from os import sys, path
import platform


if (platform.system() == "Linux"):
    WORK_PATH = "/local1/xiaoxiaol/work"
else:
    WORK_PATH = "/Users/xiaoxiaoliu/work"


p=  WORK_PATH + '/src/cell-type-analysis'
sys.path.append(p)

sys.path.append(p + '/utilities')
import morph_nfb_2_csv as nfb

sys.path.append(p + '/blast_neuron')
import blast_neuron_comp as bn

import glob
import pandas as pd




###############################################################################
#data_DIR = '/data/mat/xiaoxiaol/data/lims2/0903_filtered_ephys_qc'
data_DIR ='/Users/xiaoxiaoliu/work/data/lims2/modified_neurons/r1.0_x1.0_y1.0_z1.0_p0'
# /original stores the downloaded swc files
original_dir = data_DIR
preprocessed_dir = "/Users/xiaoxiaoliu/work/data/lims2/modified_neurons/preprocessed"


###############################################################################


if not os.path.exists(preprocessed_dir):
    os.mkdir(preprocessed_dir)

for input_swc_path in glob.glob(original_dir + "/*.swc"):
    print input_swc_path
    swc_fn = input_swc_path.split('/')[-1]

    preprocessed_swc_fn = preprocessed_dir+'/' + swc_fn
    ##bn.resample(input_swc_path, preprocessed_swc_fn)  ## due to the pw alignment, no  alignment are necessary
    bn.pre_processing(input_swc_path, preprocessed_swc_fn)

preprocessed_ANO = preprocessed_dir + "/preprocessed.ano"
bn.genLinkerFile(preprocessed_dir, preprocessed_ANO)

