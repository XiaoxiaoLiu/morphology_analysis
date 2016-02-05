__author__ = 'xiaoxiaol'
__author__ = 'xiaoxiaoliu'



import pandas as pd
import os
import sys
import platform

if (platform.system() == "Linux"):
    WORK_PATH = "/local1/xiaoxiaol/work"
else:
    WORK_PATH = "/Users/xiaoxiaoliu/work"

p =  WORK_PATH + '/src/morphology_analysis'
sys.path.append(p)
import pandas as pd
import numpy as np
import os
import blast_neuron.blast_neuron_comp as bn



def gen_txt_job_script(cmd, job_fn):
    output_dir = os.path.dirname(job_fn)
    if not os.path.exists(output_dir):
        os.system("mkdir -p  " + output_dir)
        print "create output dir: ", output_dir

    FILE = open(job_fn, 'w')

    FILE.write("%s\n" % cmd)

    FILE.close()


data_DIR = "/lustre/atlas2/nro101/proj-shared/BigNeuron/data/gold_163_all_soma_sort"
output_dir = data_DIR


neuron_distance_csv = "/data/mat/xiaoxiaol/data/big_neuron/silver/20160113_merged_gold_gt/neuron_distances_with_gold.csv"
df_nd = pd.read_csv(neuron_distance_csv)


images = np.unique(df_nd['image_file_name'])


dfg = df_nd.groupby('image_file_name')
os.system('rm -r ./txt_jobs')
os.system('mkdir ./txt_jobs')
count=0
for im in images:

     df_image = dfg.get_group(im)

     #df_image.sort(['neuron_distance'], ascending=[1], inplace=True)
     df_image=df_image.sort(['neuron_distance'])


     tmp = df_image.iloc[0]['swc_file']

     im_id = tmp.split('/')[-2]  # 2.v3dpbd

     out_dir = output_dir  + '/' + im_id.split('.')[0]

     gold_swc =  df_image.iloc[0]['gold_swc_file']

     output_gold_swc = out_dir+'/00_'+gold_swc.split('/')[-1]



     linker_fn = out_dir+'/processed/'+im_id+'.ano'

     output_eswc_path = out_dir+'/processed/consensus_p2.eswc'



     #line1 = "./start_vaa3d.sh -x linker_file_generator -f linker -i "+ out_dir+'/processed -o '+ linker_fn +' -p 1 '

     logfile = output_eswc_path+".log"
     line2 = "./start_vaa3d.sh -x consensus_swc -f consensus_swc -i " + linker_fn + " -o " + output_eswc_path + " -p 2  10 > "+logfile

     distance_log_file = output_eswc_path+".dist.log"
     line3 = "./start_vaa3d.sh -x neuron_distance -f neuron_distance -i "+ output_eswc_path +" "+ output_gold_swc +" -o "+distance_log_file

     distance_log_file = output_eswc_path+".weighted.dist.log"
     line4 = "./start_vaa3d.sh -x neuron_weighted_distance -f  neuron_weighted_distance  -i "+ output_eswc_path +" "+ output_gold_swc +" -o "+distance_log_file


     job_fn = './txt_jobs/'+str(count)+'.txt'
     FILE = open(job_fn, 'w')
     #FILE.write("%s;" % line1)
     FILE.write("%s;" % line2)
     FILE.write("%s;" % line3)
     FILE.write("%s\n" % line4)
     FILE.close()

     count = count +1


