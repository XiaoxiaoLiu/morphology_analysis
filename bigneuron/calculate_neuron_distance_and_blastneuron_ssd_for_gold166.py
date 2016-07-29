

import sys
import os
import platform

if (platform.system() == "Linux"):
    WORK_PATH = "/local1/xiaoxiaol/work"
else:
    WORK_PATH = "/Users/xiaoxiaoliu/work"

p = WORK_PATH + '/src/morphology_analysis'
sys.path.append(p)


import blast_neuron.blast_neuron_comp as bn
import recon_prescreening as rp
import glob
import pandas as pd
import os.path as path
import numpy as np


algorithm_plugin_match_csv = "/data/mat/xiaoxiaol/data/big_neuron/silver/ported_neuron_tracing_spreadsheet.csv"

###########################  preprocessing and organize data ##################################
data_DIR ="/data/mat/xiaoxiaol/data/big_neuron/silver/0401_gold163_all_soma_sort"

lookup_image_id_table_file = data_DIR +"/../image_name_lookup_table.csv"
df_lookup_table = pd.read_csv(lookup_image_id_table_file)



subdirs = [x[0] for x in os.walk(data_DIR)]


#### step 1: run neuron distance (clustering computing), generate txt jobs for xiaoxiao to run
def run_neuron_distance():
    os.system('rm ' +data_DIR+"/qsub_xx/*")

    for recon_dir in subdirs[1:]:
            folder_name = recon_dir.split('/')[-1]
            if 'processed' in  folder_name:
                      print recon_dir

                      files = glob.glob(recon_dir+'/../00_*strict.swc')
                      if len(files)>0:
                           gs_swc_file =files[0]

                      #auto recons
                      swc_files = []# glob.glob(recon_dir+'/*.strict.swc')



                      for i in range(len(swc_files)) :
                          input = swc_files[i]

                          out_fn=recon_dir+"/../processed/"+input.split('/')[-1]+'.dist.r.log'
                          print out_fn
                          if not os.path.exists(out_fn):
                                bn.run_neuron_dist(inputswc_path1=input, inputswc_path2=gs_swc_file, logfile= out_fn,
                                                     GEN_QSUB = 1, qsub_script_dir= data_DIR+"/qsub_xx", id=None)
                      #consensus
                      consensus_swc = recon_dir+'/../consensus3.strict.swc'
                      if path.exists(consensus_swc):
                             out_fn=consensus_swc+'.dist.r.log'
                             bn.run_neuron_dist(inputswc_path1=consensus_swc, inputswc_path2=gs_swc_file, logfile= out_fn,
                                                     GEN_QSUB = 1, qsub_script_dir= data_DIR+"/qsub_xx", id=None)

##

#### step2: collecting neuron distances into a csv file
def run_collecting_neuron_distance():
    output_csv=data_DIR+"/analysis_results/neuron_distance_strict.csv"

    #collect results from log files
    df_neuron_distance = pd.DataFrame(columns=('image_file_name','swc_file', 'gold_swc_file', 'algorithm', 'neuron_distance','neuron_distance_diff','neuron_distance_perc'))
    count =0
    for recon_dir in subdirs[1:]:

            folder_name = recon_dir.split('/')[-1]

            if 'processed' in  folder_name:
                      print recon_dir
                      files = glob.glob(recon_dir+'/../00_*strict.swc')
                      if len(files)>0:
                           gs_swc_file =files[0]
                      image_folder_name =  recon_dir.split('/')[-2]
                      image_id = int(image_folder_name.split(".")[0])
                      image_file_name = df_lookup_table.image_file_name[image_id-1]


                      # auto recons
                      swc_files = glob.glob(recon_dir+'/*.strict.swc')


                      for i in range(len(swc_files)) :
                          input = swc_files[i]
                          print input
                          out_fn=recon_dir+"/../processed/"+input.split('/')[-1]+'.dist.r.log'

                          algorithm = rp.matchFileToAlgorithmName( input.split('/')[-1])

                          if path.isfile(out_fn):
                            nd = bn.read_neuron_dist_log(out_fn)
                            df_neuron_distance.loc[count] = [image_file_name,input, gs_swc_file, algorithm, nd['ave'],nd['diff'],nd['perc']]
                          else:
                            print "Warning: no neuron distance log output for "+out_fn +" :output NAs."
                            df_neuron_distance.loc[count] = [image_file_name,input, gs_swc_file, algorithm, np.nan, np.nan, np.nan]
                          count= count+1

                      #consensus
                      consensus_swc= recon_dir+'/../consensus3.strict.swc'
                      consensus_swc_log = recon_dir+'/../consensus3.strict.swc.dist.r.log'
                      if path.exists(consensus_swc_log):
                            nd = bn.read_neuron_dist_log(consensus_swc_log)
                            df_neuron_distance.loc[count] = [image_file_name,consensus_swc, gs_swc_file, algorithm, nd['ave'],nd['diff'],nd['perc']]
                      else:
                            print "Warning: no neuron distance log output for "+consensus_swc +" :output NAs."
                            df_neuron_distance.loc[count] = [image_file_name,consensus_swc, gs_swc_file, algorithm, np.nan, np.nan, np.nan]
                      count = count+1

    df_neuron_distance['neuron_difference'] = df_neuron_distance['neuron_distance_diff'] *df_neuron_distance['neuron_distance_perc']

    df_neuron_distance.to_csv(output_csv,index=False)
    print "output:"+output_csv




### step 3: run blast neuron plugin ( does not need python scripting)
### all_strict.nfb need to be generated manually using vaa3d plugin
### convert all_strict.nfb to a meaningful table with meta




###$ step4: calculate blastneuron SSD into a csv file



#run_neuron_distance()
run_collecting_neuron_distance()