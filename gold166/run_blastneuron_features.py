import os

import glob

import pandas as pd

p = path.dirname(path.dirname(path.abspath(__file__)))
sys.path.append(p)
sys.path.append(p+'/blast_neuron')

import blast_neuron_comp as bn





RUN_GOLD = 0

if RUN_GOLD:
    ########################  Resample THE GOLD STANDARDS
    data_DIR = "/data/mat/xiaoxiaol/data/gold166/checked_final_swcs"
    preprocessed_dir = data_DIR +"/preprocessed"
    if  not os.path.exists(preprocessed_dir):
        os.mkdir(preprocessed_dir)


    for input_swc_path in glob.glob(data_DIR+"/*.swc"):
       print input_swc_path
       if(( os.path.getsize(input_swc_path) > 1000) and( os.path.getsize(input_swc_path) < 1024*1024*50)):
            swc_fn = input_swc_path.split('/')[-1]
            preprocessed_swc_path = preprocessed_dir+ '/'+swc_fn
            resample(input_swc_path, preprocessed_swc_path)

    preprocessed_ANO = preprocessed_dir+"/preprocessed.ano"
    bn.genLinkerFile( preprocessed_dir, preprocessed_ANO)

    ##batch computing
    feature_file =  preprocessed_dir+ "/features.nfb"
    bn.batch_compute( preprocessed_ANO,feature_file)


###############################


RUN_RESULTS = 1
if RUN_RESULTS:


data_DIR = "/data/mat/xiaoxiaol/data/gold166/gold166_results_combined"
preprocessed_dir = data_DIR +"/../preprocessed"
if  not os.path.exists(preprocessed_dir):
    os.mkdir(preprocessed_dir)


#1. sort :identify soma, connect broken trees,




#2. resample
for input_swc_path in glob.glob(data_DIR+"/*/*/*.swc"):
   print input_swc_path
   if(( os.path.getsize(input_swc_path) > 1000) and( os.path.getsize(input_swc_path) < 1024*1024*50)):
        swc_fn = "/".join (input_swc_path.split("/")[-3:])
        preprocessed_swc_path = preprocessed_dir+ '/'+swc_fn
        resample(input_swc_path, preprocessed_swc_path)



preprocessed_ANO = preprocessed_dir+"/preprocessed.ano"
genLinkerFile( preprocessed_dir, preprocessed_ANO)

##batch computing
feature_file =  preprocessed_dir+ "/features.nfb"
batch_compute( preprocessed_ANO,feature_file)






