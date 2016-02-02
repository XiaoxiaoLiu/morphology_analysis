__author__ = 'xiaoxiaoliu'


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

import utilities.morph_nfb_2_csv as nfb
import blast_neuron.blast_neuron_comp as bn
import glob
import numpy as np

import os.path as path

##################    gold standard dataset processing #########################################
###### only need to do it once
gold_dir = "/data/mat/xiaoxiaol/data/gold166/checked_final_swcs"
GOLD_CSV = "/data/mat/xiaoxiaol/data/gold166/gold.csv"  #manually generated
gold_feature_csv= "/data/mat/xiaoxiaol/data/gold166/features_with_tags.csv"

#print "resample and sort gold166s"
##rp.resample_and_sort(gold_dir+"/original",gold_dir+"/resampled",gold_dir+"/sorted")
#rp.gen_gold_feature_csv(gold_dir+"/original",output_gold_csv_file= GOLD_CSV,output_gold_feature_csv=gold_feature_csv)

algorithm_plugin_match_csv = "/data/mat/xiaoxiaol/data/big_neuron/silver/ported_neuron_tracing_spreadsheet.csv"

###########################  preprocessing and organize data ##################################
data_DIR ="/data/mat/xiaoxiaol/data/big_neuron/silver/20160113_merged_gold_gt"
ano_dir = data_DIR +"/auto_recons/ano"

anofiles = glob.glob(os.path.join(ano_dir, '*.ano'))
for ano_file in anofiles:
      logfile = ano_file+".median.log"
      #if not os.path.exists(logfile):
      #       bn.median_swc(ano_file, GEN_QSUB = 0, qsub_script_dir= data_DIR + "/median_qsub")
      print logfile

      swc_fn = bn.read_median_swc_log(ano_file, logfile)

      print swc_fn
