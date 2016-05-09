__author__ = 'xiaoxiaol'
__author__ = 'xiaoxiaol'
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




def  gen_txt_jobs(neuron_distance_csv, input_dir, output_dir):
        df_nd = pd.read_csv(neuron_distance_csv)
        images = np.unique( df_nd['image_file_name'])


        dfg = df_nd.groupby('image_file_name')
        os.system('rm -r ./txt_jobs')
        os.system('mkdir ./txt_jobs')

        count = 0
        for im in images:
             df_image = dfg.get_group(im)

             df_image = df_image.sort(['neuron_distance'])
             tmp = df_image.iloc[0]['swc_file']
             im_id = tmp.split('/')[-2]  # 2.v3dpbd
             out_dir = output_dir + '/' + im_id.split('.')[0]
             input_dir = out_dir


             output_eswc_path = out_dir+'/consensus.eswc'
             logfile = output_eswc_path+".log"
             line1 = "./start_vaa3d.sh -x consensus_swc -f consensus_swc -i " +  input_dir +"/processed/*.swc   -o " + output_eswc_path + " -p 6 5 > "+logfile

             #image_file = image_DIR+ '/'+ im[:-7]+'/'+im
             logfile= out_dir+"/median_distances.csv.log"
             line2 = "./start_vaa3d.sh -x consensus_swc -f median_swc -i "+ input_dir +"/processed/*.swc  "+ output_eswc_path +" -o "+  out_dir+"/median_distances.csv > "+logfile


             gold_swc = df_image.iloc[0]['gold_swc_file']
             gold_swc = out_dir+'/00_'+gold_swc.split('/')[-1]

             distance_log_file = output_eswc_path+".weighted.dist.log"
             line3 = "./start_vaa3d.sh -x neuron_weighted_distance -f  neuron_weighted_distance  -i "+ output_eswc_path +" "+ gold_swc +" -o "+distance_log_file


             job_fn = './txt_jobs/'+str(count)+'.txt'
             FILE = open(job_fn, 'w')
             FILE.write("%s;"  % line1)
             FILE.write("%s;"  % line2)
             FILE.write("%s\n" % line3)
             FILE.close()

             count = count +1

        os.system('tar -zcvf ./txt_jobs.tar.gz ./txt_jobs/')

        return




def  gen_qsub_jobs(neuron_distance_csv, input_dir, output_dir):
        df_nd = pd.read_csv(neuron_distance_csv)
        images = np.unique( df_nd['image_file_name'])
        QMasterV3D = "/data/mat/xiaoxiaol/work/bin/bin_vaa3d_for_clusters/vaa3d"

        dfg = df_nd.groupby('image_file_name')
        os.system('rm -r ./qsubs')
        os.system('mkdir ./qsubs')

        count = 0
        for im in images:
             df_image = dfg.get_group(im)

             df_image = df_image.sort(['neuron_distance'])
             tmp = df_image.iloc[0]['swc_file']
             im_id = tmp.split('/')[-2]  # 2.v3dpbd
             out_dir = output_dir + '/' + im_id.split('.')[0]
             input_dir = out_dir


             output_eswc_path = out_dir+'/consensus.eswc'
             output_distance_csv = out_dir+"/median_distances.csv"
        #     if os.path.exists(output_distance_csv):
         #        continue



             gold_swc = df_image.iloc[0]['gold_swc_file']
             gold_swc = out_dir+'/00_'+gold_swc.split('/')[-1]

             distance_log_file = output_eswc_path+".weighted.dist.log"
             line =  "java -jar /data/mat/xiaoxiaol/data/big_neuron/DiademMetric/DiademMetric.jar  -G "+ gold_swc +" -T "+ output_eswc_path +" -D 5 -m true -o "+distance_log_file
             #java -jar DiademMetric.jar -G ./Xiao_Xiao_test1_sn.swc  -T ./test.swc -D 5 -m true


             lines = line1+";"+line2+";"+line3
             bn.run_command_lines(lines, 1,"./qsubs", count)
             count = count +1

        return



#### main
data_DIR = "/data/mat/xiaoxiaol/data/big_neuron/silver/0401_gold163_all_soma_sort"
data_DIR = "/data/mat/xiaoxiaol/data/big_neuron/silver/gold_163_all_soma_sort_0328"
output_dir = data_DIR
neuron_distance_csv = "/data/mat/xiaoxiaol/data/big_neuron/silver/20160113_merged_gold_gt/neuron_distances_with_gold.csv"

gen_qsub_jobs(neuron_distance_csv,data_DIR, output_dir)











