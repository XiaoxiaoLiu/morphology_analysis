import matplotlib.pyplot as plt
import seaborn as sb
import os
import os.path as path
import numpy as np
import pandas as pd
import platform
import sys


if (platform.system() == "Linux"):
    WORK_PATH = "/local1/xiaoxiaol/work"
else:
    WORK_PATH = "/Users/xiaoxiaoliu/work"

p = WORK_PATH + '/src/morphology_analysis'
sys.path.append(p)

import  bigneuron.recon_prescreening as rp
import  bigneuron.plot_distances as plt_dist



import blast_neuron.blast_neuron_comp as bn


def compute_gold_standard_distances(input_dir,gold_standard_swc):

    subdirs = [x[0] for x in os.walk(input_dir)]
    count = 0
    for recon_dir in subdirs[1:]:
             print recon_dir
             folder_name = recon_dir.split('/')[-1]
    cout = 0
    cout_Cal = 0
    for dirpath, dnames, fnames in os.walk(input_dir):
          for f in fnames:
              if f.endswith(".swc"):
                 swc_f= (os.path.join(dirpath, f))
                 siz=os.path.getsize(swc_f)
                 if (siz<1024*1024):
                    log_file = swc_f + ".r.log"
                    if not os.path.exists(log_file):
                       bn.run_neuron_dist(swc_f, gold_standard_swc,log_file, 0, output_dir+"/nd")
                       cout_Cal=cout_Cal+1
                 else:
                     print "skip big file:",swc_f
                     cout= cout+1
    print cout, " big files skipped"
    print cout_Cal, " finished"



def compute_diadem_score(input_dir,gold_standard_swc):
    return

def collect_distances(input_dir, output_neuron_distance_csv):
    #collect results from log files
    df_neuron_distance = pd.DataFrame(columns=('image_file_name','swc_file', 'gold_swc_file',
                                               'algorithm',
                                               'neuron_distance_12','neuron_distance_21',
                                               'neuron_distance_ave','neuron_distance_diff',
                                               'neuron_distance_perc'))
    i=0
    for dirpath, dnames, fnames in os.walk(input_dir):
           for f in fnames:
              if f.endswith(".r.log"):
                 log_f= (os.path.join(dirpath, f))
                 image = log_f.split('/')[-2]
                 algorithm  = rp.matchFileToAlgorithmName(log_f[0:-6])
                 nd = bn.read_neuron_dist_log(log_f)
                 df_neuron_distance.loc[i] = [image,nd['input_file1'], nd['input_file2'], algorithm, nd['dist_12'],nd['dist_21'],nd['ave'],nd['diff'],
                                             nd['perc']]
                 i=i+1


    df_neuron_distance = df_neuron_distance[df_neuron_distance['neuron_distance_ave'] != -1]  # empty no

    df_neuron_distance.to_csv(output_neuron_distance_csv, index=False)

    return




def plot_neuron_distance_by_groups(neuron_distance_csv, outputDir,algorithms,CASE_BY_CASE_PLOT = 0):
    #neuron_distance_csv = data_DIR + "/neuron_distance.r.csv"
    #outputDir = data_DIR + "/neuron_dist_plots_r"
    df_nd = pd.read_csv(neuron_distance_csv)

    all_images = np.unique(df_nd.image_file_name)
    if not path.exists(outputDir):
        os.mkdir(outputDir)

    dfg = df_nd.groupby('image_file_name')
    sample_size_per_img=[]
    for img in all_images:
        sample_size_per_img.append(dfg.get_group(img).shape[0])

    plt.figure()
    sb.set_context("talk")

    a=sb.barplot(x='image_file_name', y='neuron_distance_ave', data=df_nd,order=all_images)
    a.set_xticklabels(['%s ($n$=%d )'%(all_images[i], sample_size_per_img[i]) for i in range(all_images.size) ])

    plt.xticks(rotation="90")
    plt.ylabel('Average Neuron Distance')
    plt.subplots_adjust(left=0.1, bottom=0.5, top=0.9)
    plt.savefig(outputDir + '/group_by_image_average_neuron_distance.png', format='png')
    plt.show()
    plt.close()

    plt.figure()
    g=sb.lmplot(x="SNR", y="neuron_distance_ave", col="correlated", hue="sigma", data=df_nd,legend_out=False)

    g = (g.set_axis_labels("SNR", "Average Neuron Distance")
           .set( xticks=[1,2,3,5,10,20])
           .fig.subplots_adjust(wspace=.02))

    #g.set_xticklabels([1,2,3,5,10,20],['1','2','3','5','10','noise_free'])
    plt.subplots_adjust(left=0.1, bottom=0.3, top=0.9)
    plt.savefig(outputDir + '/SNR_average_neuron_distance.png', format='png')
    plt.show()
    plt.close()

    return







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
#run_consensus(data_DIR, output_dir)

#compare distances
gold_dir = "/data/mat/xiaoxiaol/data/gold166/checked_final_swcs"
GOLD_CSV = "/data/mat/xiaoxiaol/data/gold166/gold.csv"  #manually generated
gold_feature_csv= "/data/mat/xiaoxiaol/data/gold166/features_with_tags.csv"

output_neuron_distance_csv = output_dir+"/neuron_distance.csv"
input_dir = data_DIR
#compute_gold_standard_distances(input_dir=data_DIR,gold_standard_swc=data_DIR+"/../Cell1_40x_1.0umstep_H15.06.005.07.02_cell1.oif.files.v3dpbd_stamp_2015_06_16_13_31_radiusreestimated.swc")
#collect_distances(input_dir=data_DIR,output_neuron_distance_csv=output_dir+"/neuron_distance.csv")


#plt_dist.plot_neuron_distance(output_neuron_distance_csv,input_dir,None,0)
df_nd = pd.read_csv(output_neuron_distance_csv)
df_new_nd= pd.DataFrame(columns=['image_file_name','SNR','sigma','correlated','swc_file','gold_swc_file','algorithm','neuron_distance_12','neuron_distance_21','neuron_distance_ave','neuron_distance_diff','neuron_distance_perc'])
df_new_nd[df_nd.columns] = df_nd

for i in range(len(df_nd)):
    img = df_nd.iloc[i]['image_file_name']

    if img== "noise_free":
        df_new_nd.loc[i,'SNR']  = 20
        df_new_nd.loc[i,'correlated']=  'No'
        df_new_nd.loc[i,'sigma']= 0
    else:
        df_new_nd.loc[i,'SNR']  = img.split('_')[1]
        if img.split('_')[2] == 'correlated':
             df_new_nd.loc[i,'correlated']=  'Yes'
        else:
             df_new_nd.loc[i,'correlated']=  'No'
        df_new_nd.loc[i,'sigma']= 0
        if img.split('_')[2] == "correlated":
               df_new_nd.loc[i,'sigma'] = img.split('_')[-1]

output_neuron_distance_extend_csv=output_dir+"/neuron_distance_parse.csv"
df_new_nd.to_csv(output_neuron_distance_extend_csv)
plot_neuron_distance_by_groups(output_neuron_distance_extend_csv,input_dir,None,0)



