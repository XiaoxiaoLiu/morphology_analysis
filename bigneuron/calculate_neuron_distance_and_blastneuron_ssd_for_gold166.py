

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
import glob
import pandas as pd
import os.path as path


algorithm_plugin_match_csv = "/data/mat/xiaoxiaol/data/big_neuron/silver/ported_neuron_tracing_spreadsheet.csv"

###########################  preprocessing and organize data ##################################
data_DIR ="/data/mat/xiaoxiaol/data/big_neuron/silver/0401_gold163_all_soma_sort"

lookup_image_id_table_file = data_DIR +"/../image_name_lookup_table.csv"


#### step 1: run neuron distance (clustering computing), generate txt jobs for xiaoxiao to run
os.system('rm ' +data_DIR+"/qsub_xx/*")
subdirs = [x[0] for x in os.walk(data_DIR)]

for recon_dir in subdirs[1:]:
        folder_name = recon_dir.split('/')[-1]
        if 'processed' in  folder_name:
                  print recon_dir
                  folder_name = recon_dir.split('/')[-1]
                  files = glob.glob(recon_dir+'/../00_*strict.swc')
                  if len(files)>0:
                       gs_swc_file =files[0]

                  swc_files = glob.glob(recon_dir+'/*.strict.swc')
                  for i in range(len(swc_files)) :
                      input = swc_files[i]

                      out_fn=recon_dir+"/../processed/"+input.split('/')[-1]+'.dist.r.log'
                      print out_fn
                      if not os.path.exists(out_fn):
                            bn.run_neuron_dist(inputswc_path1=input, inputswc_path2=gs_swc_file, logfile= out_fn,
                                                 GEN_QSUB = 1, qsub_script_dir= data_DIR+"/qsub_xx", id=None)


#### step2: collecting neuron distances into a csv file

output_csv=data_DIR+"/analysis_results/neuron_distance_strict.csv"

#collect results from log files
df_neuron_distance = pd.DataFrame(columns=('image_file_name','swc_file', 'gold_swc_file', 'algorithm', 'neuron_distance','neuron_distance_diff','neuron_distance_perc'))
for recon_dir in subdirs[1:]:
        folder_name = recon_dir.split('/')[-1]
        if 'processed' in  folder_name:
                  print recon_dir
                  folder_name = recon_dir.split('/')[-1]
                  files = glob.glob(recon_dir+'/../00_*strict.swc')
                  if len(files)>0:
                       gs_swc_file =files[0]

                  swc_files = glob.glob(recon_dir+'/*.strict.swc')
                  for i in range(len(swc_files)) :
                      input = swc_files[i]

                      out_fn=recon_dir+"/../processed/"+input.split('/')[-1]+'.dist.r.log'

                      image_file_name =
                      algorithm =
                      if path.isfile(out_fn):
                        nd = bn.read_neuron_dist_log(out_fn)
                        df_neuron_distance.loc[i] = [image_file_name,input, gs_swc_file, algorithm, nd['ave'],nd['diff'],nd['perc']]
                      else:
                        print "Warning: no neuron distance log output for "+tmp +" :output NAs."
                        df_neuron_distance.loc[i] = [image_file_name,input, gs_swc_file, algorithm, np.nan, np.nan, np.nan]


df_neuron_distance['neuron_difference'] = df_neuron_distance['neuron_distance_diff'] *df_neuron_distance['neuron_distance_perc']

df_neuron_distance.to_csv(output_csv,index=False)
print "output:"+output_csv



### step 3: run blast neuron plugin ( does not need python scripting)
### all_strict.nfb need to be generated manually using vaa3d plugin
### convert all_strict.nfb to a meanginful table with meta


###$ step4: calculate blastneuron SSD into a csv file



