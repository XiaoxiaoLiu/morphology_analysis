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

import os




def gen_txt_job_script(cmd, job_fn):
    output_dir = os.path.dirname(job_fn)
    if not os.path.exists(output_dir):
        os.system("mkdir -p  " + output_dir)
        print "create output dir: ", output_dir

    FILE = open(job_fn, 'w')

    FILE.write("%s\n" % cmd)

    FILE.close()






def generate_txt_files(fn_list, data_DIR, output_dir, rerun_missing_data_dir=None, startJobId = 0):

    df_nd = pd.read_csv(fn_list)
    images = df_nd['image_file_name']


    dfg = df_nd.groupby('image_file_name')

    count = startJobId

    for im in images:
         out_dir = output_dir
         input_dir = data_DIR+'/'+im[:-7]

         output_eswc_path = out_dir+'/'+im+'_consensus.eswc'
         if (rerun_missing_data_dir is not None):
             #test_eswc_path = rerun_missing_data_dir+'/'+im+'_consensus.eswc'
             test_distance_path = rerun_missing_data_dir+'/'+im+'_median_distances.csv'
             #print test_distance_path
             if os.path.exists(test_distance_path):
                 continue

         logfile = output_eswc_path+".log"
         line1 = "./start_vaa3d.sh -x consensus_swc -f consensus_swc -i " +  input_dir +"/*.swc   -o " + output_eswc_path + " -p 6 5 1 > "+logfile

         #image_file = image_DIR+ '/'+ im[:-7]+'/'+im
         output_csv =  out_dir+"/"+im+"_median_distances.csv"
         logfile = output_csv+".log"
         line2 = "./start_vaa3d.sh -x consensus_swc -f median_swc -i "+ output_eswc_path +"_SelectedNeurons.ano  "+ output_eswc_path +" -o "+  output_csv+" > "+logfile



         job_fn = './txt_jobs/'+str(count)+'.txt'
         FILE = open(job_fn, 'w')
         FILE.write("%s;" % line1)
         FILE.write("%s\n" % line2)
         FILE.close()

         count = count +1



    return count





def run_janelia_set(set_type = "set1_original", tag='consensus_0725'):


    os.system('rm -r ./txt_jobs')
    os.system('mkdir ./txt_jobs')


    rerun_missing_data_dir= None
    endJobId = 0

    if set_type =='set1_original':
      data_DIR = "/lustre/atlas2/nro101/proj-shared/BigNeuron/data/Janelia/set1_extract_single/reconstructions_for_img_nopreproprcessing"
      output_dir =  "/lustre/atlas2/nro101/proj-shared/BigNeuron/data/Janelia/set1_extract_single/"+tag

      fn_list = '/data/mat/xiaoxiaol/data/big_neuron/jen1_image_file_name_list.csv'
      #image_DIR="/lustre/atlas2/nro101/proj-shared/BigNeuron/data/Janelia/set1_extract_single/img_nopreproprcessing"

      rerun_missing_data_dir = "/mnt/BigNeuron/data/Janelia/set1_extract_single/"+tag

      endJobId = generate_txt_files(fn_list, data_DIR, output_dir, rerun_missing_data_dir,  0)


    if set_type =='set1_smooth':
      data_DIR = "/lustre/atlas2/nro101/proj-shared/BigNeuron/data/Janelia/set1_extract_single/reconstructions_for_img_anisosmooth"
      output_dir =  "/lustre/atlas2/nro101/proj-shared/BigNeuron/data/Janelia/set1_extract_single/"+tag+"_anisosmooth"

      fn_list = '/data/mat/xiaoxiaol/data/big_neuron/jen1_image_file_name_list.csv'
      #image_DIR="/lustre/atlas2/nro101/proj-shared/BigNeuron/data/Janelia/set1_extract_single/img_nopreproprcessing"

      rerun_missing_data_dir = "/mnt/BigNeuron/data/Janelia/set1_extract_single/"+tag+"_anisosmooth"

      endJobId = generate_txt_files(fn_list, data_DIR, output_dir, rerun_missing_data_dir,  endJobId)


    if  set_type =='set2_original':
      data_DIR = "/lustre/atlas2/nro101/proj-shared/BigNeuron/data/Janelia/set2_accepted_single/reconstructions_for_img_nopreproprcessing"
      output_dir =  "/lustre/atlas2/nro101/proj-shared/BigNeuron/data/Janelia/set2_accepted_single/"+tag
      fn_list =  '/data/mat/xiaoxiaol/data/big_neuron/jen2_image_file_name_list.csv'
      rerun_missing_data_dir = "/mnt/BigNeuron/data/Janelia/set2_accepted_single/"+tag
      endJobId = generate_txt_files(fn_list, data_DIR, output_dir, rerun_missing_data_dir,  endJobId)


    if  set_type =='set2_smooth':
      data_DIR = "/lustre/atlas2/nro101/proj-shared/BigNeuron/data/Janelia/set2_accepted_single/reconstructions_for_img_anisosmooth"
      output_dir =  "/lustre/atlas2/nro101/proj-shared/BigNeuron/data/Janelia/set2_accepted_single/"+tag+"_anisosmooth"
      fn_list = '/data/mat/xiaoxiaol/data/big_neuron/jen2_image_file_name_list.csv'
      rerun_missing_data_dir = "/mnt/BigNeuron/data/Janelia/set2_accepted_single/"+tag+"_anisosmooth"
      endJobId = generate_txt_files(fn_list, data_DIR, output_dir, rerun_missing_data_dir,  endJobId)

    os.system('pwd')
    os.system('tar -zcvf ./txt_jobs.tar.gz ./txt_jobs/')





run_janelia_set(set_type = "set1_original", tag="consensus_0725")
# run_janelia_set(set_type = "set1_smooth")
# run_janelia_set(set_type = "set2_original")
# run_janelia_set(set_type = "set2_smooth")