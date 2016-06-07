import os
import sys
import platform

if (platform.system() == "Linux"):
    WORK_PATH = "/local1/xiaoxiaol/work"
else:
    WORK_PATH = "/Users/xiaoxiaoliu/work"

p =  WORK_PATH + '/src/morphology_analysis'
sys.path.append(p)


import blast_neuron.blast_neuron_comp as bn



def  run_consensus(input_dir, output_dir):
       # each subfolder contains the reconstructions from one image
       # walk through all folders
        subdirs = [x[0] for x in os.walk(input_dir)]
        count = 0
        for recon_dir in subdirs[1:]:
             print recon_dir
             folder_name = recon_dir.split('/')[-1]
             bn.consensus(input_ano_path=recon_dir+"/*.swc",output_eswc_path=output_dir+"/"+folder_name+"/consensus.eswc",
                          vote_threshold=3, max_cluster_distance = 5,resampling = 0, remove_outlier=1, GEN_QSUB = 1, qsub_script_dir= output_dir+"/qsub", id=None)

        #      #image_file = image_DIR+ '/'+ im[:-7]+'/'+im
        #      logfile= out_dir+"/median_distances.csv.log"
        #      line2 =  QMasterV3D+" -x consensus_swc -f median_swc -i "+ output_eswc_path +"_SelectedNeurons.ano  "+ output_eswc_path +" -o "+  out_dir+"/median_distances.csv > "+logfile
        #
        #
        #      gold_swc = df_image.iloc[0]['gold_swc_file']
        #      gold_swc = out_dir+'/00_'+gold_swc.split('/')[-1]
        #
        #      distance_log_file = output_eswc_path+".weighted.dist.log"
        #      # if os.path.exists(distance_log_file):
        #      #     continue
        #      line3 =  QMasterV3D+" -x neuron_weighted_distance -f  neuron_weighted_distance  -i "+ output_eswc_path +" "+ gold_swc +" -o "+distance_log_file
        #
        return



#### main
data_DIR = "/data/mat/xiaoxiaol/data/big_neuron/simulation/reconstructions"

output_dir = data_DIR
run_consensus(data_DIR, output_dir)











