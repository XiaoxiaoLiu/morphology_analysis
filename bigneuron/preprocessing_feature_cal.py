import os
from os import path, sys
import glob
import pandas as pd


WORK_PATH = "/Users/xiaoxiaoliu/work"
p =  WORK_PATH + '/src/morphology_analysis'
sys.path.append(p)
import blast_neuron.blast_neuron_comp as bn

data_DIR =WORK_PATH+"/data/20151030_rhea_reconstructions_for_allen300_silver_set"
original_dir = data_DIR + "/auto_recons"
preprocessed_dir = data_DIR +"/resampled"
sorted_dir = data_DIR +"/sorted"
if  not os.path.exists(preprocessed_dir):
      os.mkdir(preprocessed_dir)
if  not os.path.exists(sorted_dir):
      os.mkdir(sorted_dir)

###############################  run on cluster ################
RUN_RESULTS = 1
GEN_QSUB = 0
failure_file = open(data_DIR +"/failurecases_via_size.txt","w")
i = 0
if RUN_RESULTS:
  for input_swc_path in glob.glob(original_dir+"/*/*/*.swc"):
     i = i+1
     print "\n\n "
     print i
     print "resample and sort : "+ input_swc_path
     if(( os.path.getsize(input_swc_path) > 1000) and( os.path.getsize(input_swc_path) < 1024*1024/2)):
        if not "tmp_cache_img"  in input_swc_path:   # skip the tmp files
          swc_fn = "/".join (input_swc_path.split("/")[-3:]) # to keep the subfolder structure

          # resample
          preprocessed_swc_path = preprocessed_dir+ '/'+swc_fn
          bn.resample(input_swc_path, preprocessed_swc_path,3,GEN_QSUB,data_DIR+'/qsub/resample')  # generate QSUB scripts

          # sort
          sorted_swc_path = sorted_dir+ '/'+swc_fn
          bn.sort_swc(preprocessed_swc_path, sorted_swc_path,GEN_QSUB,data_DIR+'/qsub/sort')
     else:
        failure_file.write(input_swc_path)

failure_file.close()

FEATURE_CALC = 0
if FEATURE_CALC:
  sorted_ANO = sorted_dir+"/sorted.ano"
  bn.genLinkerFile( sorted_dir, sorted_ANO)

  ##batch computing
  feature_file =  sorted_dir+ "/features.nfb"
  bn.batch_compute (sorted_ANO,feature_file)
