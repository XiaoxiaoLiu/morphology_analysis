__author__ = 'xiaoxiaol'


import pandas as pd
import numpy as np
import os


import blast_neuron.blast_neuron_comp as bn
import matplotlib.pyplot as plt
import seaborn as sb




###############################
data_DIR = "/data/mat/xiaoxiaol/data/20151030_rhea_reconstructions_for_allen300_silver_set"

original_dir = data_DIR + "/auto_recons"
preprocessed_dir = data_DIR +"/resampled"
sorted_dir = data_DIR +"/sorted"

for dirpath, dirnames, files in os.walk(original_dir):
    bn.genLinkerFile( dirpath, original_dir+'/ano/'+dirpath.split('/')[-2]+'.'+dirpath.split('/')[-1]+'_linker.ano')

# move the group anos into group/


# filtering extreme large swc files ( deu to tracing errors)



#


#  ##batch computing
# feature_file =  sorted_dir+ "/features.nfb"
# bn.batch_compute (sorted_ANO,feature_file)
