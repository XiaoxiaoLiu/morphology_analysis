__author__ = 'xiaoxiaoliu'

import os
from os import path, sys
import glob
import pandas as pd
import numpy as np


WORK_PATH = "/Users/xiaoxiaoliu/work"
p =  WORK_PATH + '/src/morphology_analysis'
sys.path.append(p)
import blast_neuron.blast_neuron_comp as bn
import utilities.morph_nfb_2_csv as nfb

data_DIR = WORK_PATH+"/data/20151030_rhea_reconstructions_for_allen300_silver_set"
original_dir = data_DIR + "/auto_recons"
preprocessed_dir = data_DIR +"/79/resampled"
sorted_dir = data_DIR +"/79/sorted"
print sorted_dir
if  not os.path.exists(preprocessed_dir):
      os.mkdir(preprocessed_dir)
if  not os.path.exists(sorted_dir):
      os.mkdir(sorted_dir)


# generate ano files into a single top folder with absoluate paths to swc files reconstructed for one image
def generate_ano_files_two_level_dir(my_dir):
    if  not os.path.exists(my_dir+'/ano'):
       os.mkdir(my_dir+'/ano')

    topdirs = glob.glob(os.path.join(my_dir, '*'))
    # print topdirs
    for subdir in topdirs:
        #print subdir
        for datasetdir in  glob.glob(os.path.join(subdir, '*')):
                 #print datasetdir

                 dataset_path = os.path.abspath(datasetdir)
                 bn.genLinkerFile( dataset_path, my_dir+'/ano/'+subdir.split('/')[-1]+'.'+datasetdir.split('/')[-1]+'.recons.ano')
    return


def recon_table_gen(data_root):
        generate_ano_files_two_level_dir(data_root)  # assuming two-level data archives
        anofiles = glob.glob(os.path.join(data_root+'/ano/', '*.ano'))
        print "there are "+str(len(anofiles))+ " datasets"

        df_silver = pd.DataFrame()

        algorithmList = []
        imageList = []
        swc_list=[]
        for anofile in anofiles:
                f_ano = open(anofile, "r+")
                for line in f_ano:

                    swc_file = line.split('SWCFILE=')[-1]
                    swc_file = swc_file.split('\n')[0]
                    fn = swc_file.split('/')[-1]
                    fn = fn.split('\n')[0]
                    algorithm = fn.split('v3dpbd_')[-1]
                    algorithm = algorithm.split('.')[0]
                    if "app1" in algorithm:   # for patterns like *x245_y234_z234_app1.swc
                          algorithm = "app1"
                    if "app2" in algorithm:
                          algorithm = "app2"
                    if  "spanningtree" in algorithm: # fastmarching_spanningtree is too long
                          algorithm = "spanningtree"
                    image = fn.split('.v3dpbd')[0]
                    image = image.split('sorted_')[-1]  # for sorted_* swc_files
                    algorithmList.append(algorithm)
                    imageList.append(image)
                    swc_list.append(swc_file)

        df_silver['algorithm'] = pd.Series(algorithmList)
        df_silver['image'] = pd.Series(imageList)
        df_silver['swc_file'] =pd.Series(swc_list)
        df_silver.to_csv(data_root + '/recon_table.csv', index=False)
        return




###############################  run on cluster ################

df_shared=pd.read_csv(original_dir+"/shared.csv")

RUN_RESULTS = 0
GEN_QSUB = 0
failure_file = open(data_DIR +"/79_failurecases_via_size.txt","w")
i = 0
if RUN_RESULTS:
      for input_swc_path in df_shared.swc_file:
             i = i+1
             print "\n\n "
             print i
             print "resample and sort : "+ input_swc_path
             if(( os.path.getsize(input_swc_path) > 1000) and( os.path.getsize(input_swc_path) < 1024*1024/4)):
                    if not "tmp_cache_img"  in input_swc_path:   # skip the tmp files
                      swc_fn = "/".join (input_swc_path.split("/")[-3:]) # to keep the subfolder structure

                      # resample
                      preprocessed_swc_path = preprocessed_dir+ '/'+swc_fn
                      bn.resample(input_swc_path, preprocessed_swc_path,3,GEN_QSUB,data_DIR+'/qsub/resample')  # generate QSUB scripts

                      # sort
                      sorted_swc_path = sorted_dir+ '/'+swc_fn
                      bn.sort_swc(preprocessed_swc_path, sorted_swc_path,GEN_QSUB,data_DIR+'/qsub/sort')
             else:
                failure_file.write(input_swc_path+"\n")
      failure_file.close()
      recon_table_gen(sorted_dir)



# run consensus skeleton algorithm
RUN_CONSENSUS = 1
if RUN_CONSENSUS:
    filesano = glob.glob(sorted_dir+'/ano/*.ano')
    output_dir = data_DIR+"/consensus"
    if  not os.path.exists(output_dir):
         os.mkdir(output_dir)
    for file in filesano:
         print file
         #bn.consensus(file,output_dir+'/'+file.split('/')[-1]+'.consensus.eswc',1)
         bn.consensus(file,output_dir+'/'+file.split('/')[-1]+'.consensus_mst.eswc',0)



RUN_VoteMap = 0
if RUN_VoteMap:
    anofiles = glob.glob(os.path.join(sorted_dir+'/ano/', '*.ano'))
    print "there are "+str(len(anofiles))+ " sorted datasets"
    output_dir = data_DIR+"/votemaps"
    if  not os.path.exists(output_dir):
         os.mkdir(output_dir)
    for anofile in anofiles:
         bn.vote_map(anofile,output_dir+'/'+anofile.split('/')[-1]+'.votemap.tif')

RUN_Median = 0
if RUN_Median:
    anofiles = glob.glob(os.path.join(sorted_dir+'/ano/', '*.ano'))
    print "there are "+str(len(anofiles))+ " sorted datasets"
    output_dir = data_DIR+"/medians"
    if  not os.path.exists(output_dir):
         os.mkdir(output_dir)
    for anofile in anofiles:
         bn.median_swc(anofile,output_dir+'/'+anofile.split('/')[-1]+'.median.swc')


COMPARE_NEURON_DISTANCE = 0
if COMPARE_NEURON_DISTANCE:
    df_gold = pd.read_csv(WORK_PATH+"/data/gold79/sorted/gold.csv")
    df_results = pd.read_csv((sorted_dir+'/recon_table.csv'))

    # find common sets
    # images_gold = np.unique(df_gold.image)
    # images_results = np.unique(df_results.image)
    # final_images = np.intersect1d(images_gold, images_results)

    for i in range(df_gold.image.size):
        imageName = df_gold.image[i]
        df_image = df_results[df_results.image == imageName]
        df_gold_image = df_gold[df_gold.image == imageName].iloc[0]

        for j in range(df_image.shape[0]):
            df_swc = df_image.iloc[j]
            bn.run_neuron_dist(df_gold_image.gold_swc_file, df_swc.swc_file, df_swc.swc_file + ".r.log", 0, "nd")

COLLECT_NEURON_DISTANCE = 0
if COLLECT_NEURON_DISTANCE:
    df_gold = pd.read_csv(WORK_PATH+"/data/gold79/sorted/gold.csv")
    df_results = pd.read_csv((sorted_dir+'/recon_table.csv'))
    df_neuron_distance = pd.DataFrame(columns=('swc_file', 'gold_swc_file', 'algorithm', 'neuron_distance'))
    idx = 0
    for i in range(df_gold.image.size):
        imageName =  df_gold.image[i]
        df_image = df_results[df_results.image == imageName]
        df_gold_image = df_gold[df_gold.image == imageName].iloc[0]

        for j in range(df_image.shape[0]):
            df_swc = df_image.iloc[j]

            logfile = df_swc.swc_file + ".r.log"
            if path.isfile(logfile):
                nd = bn.read_neuron_dist_log(logfile)
                algorithm = df_swc.algorithm
                df_neuron_distance.loc[idx] = [df_swc.swc_file, df_gold_image.gold_swc_file,
                                               algorithm, nd['ave']]
                idx = idx + 1
    df_neuron_distance.to_csv(data_DIR + "/neuron_distance.r.csv", index=False)



FEATURE_CALC = 0
if FEATURE_CALC:
          out_sorted_ANO = sorted_dir+"/sorted.ano"
          bn.genLinkerFile( sorted_dir, out_sorted_ANO)

          ##batch computing
          out_feature_file =  sorted_dir+ "/features.nfb"
          bn.batch_compute (out_sorted_ANO,out_feature_file)

FEATURE_CALC = 0
if FEATURE_CALC:
          gold_sorted_dir =WORK_PATH+"/data/gold79/sorted"
          out_sorted_ANO = gold_sorted_dir+"/sorted.ano"
          bn.genLinkerFile(gold_sorted_dir , out_sorted_ANO)

          ##batch computing
          out_feature_file =  gold_sorted_dir + "/features.nfb"
          bn.batch_compute (out_sorted_ANO,out_feature_file)

          nfb.generateALLFeatureCSV_gold166(out_feature_file, gold_sorted_dir +'/features_with_tags.csv')



def SSD(feature_array1, feature_array2):
    diff_v = np.array(feature_array1) - np.array(feature_array2)
    ssd = np.sum(np.abs(diff_v) ** 2)
    return ssd


BLASTNEURON_PLOT=0
if BLASTNEURON_PLOT:
    results_csv = data_DIR + "/79/sorted/features_with_tags.csv"
    df_results = pd.read_csv(results_csv)
    df_gold = pd.read_csv(WORK_PATH+"/data/gold79/sorted/features_with_tags.csv")


    ########## filter out spanning tree, which generates huge numbers
    # df_results = df_results[df_results.algorithm != 'spanningtree']


    # find common sets
    images_gold = np.unique(df_gold.image)
    images_results = np.unique(df_results.image)

    # common set of images
    final_images = np.intersect1d(images_gold, images_results)

    feature_cols = [u'num_nodes', u'soma_surface', u'num_stems', u'num_bifurcations', u'num_branches',
                    u'num_of_tips', u'overall_width', u'overall_height', u'overall_depth', u'average_diameter',
                    u'total_length', u'total_surface', u'total_volume', u'max_euclidean_distance',
                    u'max_path_distance', u'max_branch_order', u'average_contraction', u'average fragmentation',
                    u'parent_daughter_ratio', u'bifurcation_angle_local', u'bifurcation_angle_remote', u'moment1',
                    u'moment2', u'moment3', u'moment4', u'moment5', u'moment6', u'moment7', u'moment8', u'moment9',
                    u'moment10', u'moment11', u'moment12', u'moment13', u'avgR']

    selected_cols = [ u'num_stems', u'num_bifurcations', u'num_branches',u'num_of_tips', u'overall_width', u'overall_height', u'overall_depth',
                    u'total_length', u'max_euclidean_distance',u'max_path_distance', u'max_branch_order', u'average_contraction', u'average fragmentation',u'parent_daughter_ratio', u'bifurcation_angle_local', u'bifurcation_angle_remote', u'moment1',
                    u'moment2', u'moment3', u'moment4', u'moment5', u'moment6', u'moment7', u'moment8', u'moment9',u'moment10', u'moment11', u'moment12', u'moment13', u'avgR']

    my_cols1 = [u'image', u'algorithm', u'swc_file']
    my_cols1.extend(selected_cols)
    df_results_s = df_results[my_cols1]
    del df_results
    my_cols2 = [u'image']
    my_cols2.extend(selected_cols)
    df_gold_s = df_gold[my_cols2]
    del df_gold

    # calcualte std for each image, calculate normalized results
    df_normalized = pd.DataFrame(columns=my_cols1)
    df_gold_normalized = pd.DataFrame(columns=my_cols2)

    for i in range(final_images.size):
        imageName = final_images[i]
        df_image = df_results_s[df_results_s.image == imageName]

        if df_image.shape[0] > 5:  # too few samples
            df_image[selected_cols] = df_image[selected_cols].astype(float)  # some moment values are interpreted as strings
            df_gold_image = df_gold_s[df_gold_s.image == imageName]

            df_normalized_per_image = pd.DataFrame(columns=my_cols1)
            df_gold_normalized_per_image = pd.DataFrame(columns=my_cols2)

            df_normalized_per_image[[u'image', u'algorithm', u'swc_file']] = df_image[[u'image', u'algorithm', u'swc_file']]
            df_gold_normalized_per_image[u'image'] = df_gold_image.image
            print imageName

            for col in selected_cols:
                if (df_image[col].std() == 0.0 ):
                    print "std = 0!", col, " ", imageName
                df_normalized_per_image[col] = (df_image[col] - df_image[col].median() ) / (df_image[col].std() + 0.0000001)
                df_gold_normalized_per_image[col] = ( df_gold_image[col] - df_image[col].median() ) / (df_image[col].std() + 0.0000001)

            # append the results for this image
            df_normalized = df_normalized.append(df_normalized_per_image, ignore_index=True)
            df_gold_normalized= df_gold_normalized.append(df_gold_normalized_per_image, ignore_index=True)


    #calculated sit
    ssd_metrics = []
    for rowIdx in range(df_normalized.shape[0]):
        imageName = df_normalized.image[rowIdx]
        gold_df = df_gold_normalized.loc[df_gold_normalized.image == imageName]

        #normalize restuls with mean and std
        result = df_normalized.iloc[rowIdx]
        sumsquare = SSD(result[selected_cols], gold_df[selected_cols])
        ssd_metrics.append(sumsquare)

    df_normalized['SSD'] = ssd_metrics
    ## reordering columns
    df_normalized.to_csv(data_DIR + '/normalized_bn_dist.csv')

