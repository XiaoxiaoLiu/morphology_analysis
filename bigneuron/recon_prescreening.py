__author__ = 'xiaoxiaol'


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


import utilities.morph_nfb_2_csv as nfb
import blast_neuron.blast_neuron_comp as bn
import glob
import numpy as np

import os.path as path
######################################################
LOOKUP_TABLE_FILE = "/data/mat/xiaoxiaol/data/big_neuron/silver/ported_neuron_tracing_spreadsheet.csv"
##################################




def map_better_algorithm_name(alg,algorithm_plugin_match_csv=LOOKUP_TABLE_FILE):

    df_check_table = pd.read_csv(algorithm_plugin_match_csv)
    keys = df_check_table['algorithm']
    values = df_check_table['better_algorithm_name']
    algorithm_name_mapping = dict(zip(keys, values))
    return algorithm_name_mapping[alg]







def  matchFileToAlgorithmName(file_name,lookup_table_file=LOOKUP_TABLE_FILE):
     if file_name.find("consensus") >-1:
         return 'consensus'

     if file_name.find("v3dpbd") >-1:
             tmp = file_name.split('v3dpbd_')[-1]
     else:
             tmp = file_name.split('v3draw_')[-1]
     al_string = tmp.split('.swc')[0]


     if al_string.startswith('anisodiff.raw_'):
         al_string = al_string.split('anisodiff.raw_')[-1]


     df_lookup = pd.read_csv(lookup_table_file)
     keys = df_lookup['swc_file_label']
     values = df_lookup['algorithm']
     name_mapping = dict(zip(keys, values))

     alg=None
     unique_keys = np.unique(keys)
     for key in unique_keys:
        if key in al_string:
               alg= name_mapping[key]
     return  alg





def SSD(feature_array1, feature_array2):
    ssd = -1
    if  feature_array1.size >0 and feature_array2.size >0:
        diff_v = np.array(feature_array1) - np.array(feature_array2)
        ssd = np.sum(np.abs(diff_v) ** 2)
    return ssd


GL_FEATURE_TAGS = np.array(
['num_nodes', 'soma_surface', 'num_stems', 'num_bifurcations', 'num_branches', 'num_of_tips', 'overall_width',
 'overall_height', 'overall_depth', 'average_diameter', 'total_length', 'total_surface', 'total_volume',
 'max_euclidean_distance', 'max_path_distance', 'max_branch_order', 'average_contraction', 'average fragmentation',
 'parent_daughter_ratio', 'bifurcation_angle_local', 'bifurcation_angle_remote'])
GMI_FEATURE_TAGS = np.array(
['moment1', 'moment2', 'moment3', 'moment4', 'moment5', 'moment6', 'moment7', 'moment8', 'moment9', 'moment10',
 'moment11', 'moment12', 'moment13', 'avgR'])

def readDBFeatures(FEATURE_FILE):
    # TODO: detect nan values
    glf_featureList = []  # each row is a feature vector
    gmi_featureList = []
    swc_file_nameList = []
    with open(FEATURE_FILE, 'r') as  f:
        for fn_line in f:  # ignore the SWCFILE=* line

            swc_file = fn_line[8:].strip()
            swc_file_nameList.append(swc_file)

            line_globalFeature = (f.next()).strip()
            glf = map(float, line_globalFeature.split('\t'))
            glf_featureList.append(glf)

            line_GMI = (f.next()).strip()
            gmi = map(float, line_GMI.split('\t'))
            gmi_featureList.append(gmi)

    return swc_file_nameList, np.array(glf_featureList), np.array(gmi_featureList)


def generateALLFeatureCSV(feature_file, feature_csv_file):
    swc_file_nameList, glFeatures, gmiFeatures = readDBFeatures(feature_file)

    allFeatures = np.append(glFeatures, gmiFeatures, 1)
    allColumns = np.append(GL_FEATURE_TAGS, GMI_FEATURE_TAGS, 0)

    df = pd.DataFrame(allFeatures, columns=allColumns)

    df['swc_file'] = pd.Series(swc_file_nameList, index=df.index)

    allColumns = np.append(np.array(['swc_file']), allColumns, 0)

    df = df[allColumns]

    df.to_csv(feature_csv_file, index=False)

    print 'output all feature csv file to :', feature_csv_file
    return






def generateALLFeatureCSV_gold166(feature_file, feature_csv_file):

    swc_file_nameList, glFeatures, gmiFeatures = readDBFeatures(feature_file)

    allFeatures = np.append(glFeatures, gmiFeatures, 1)
    allColumns = np.append(GL_FEATURE_TAGS, GMI_FEATURE_TAGS, 0)

    df = pd.DataFrame(allFeatures, columns=allColumns)

    df['swc_file'] = pd.Series(swc_file_nameList, index=df.index)

    algorithmList = []
    imageList = []
    for swc_file in swc_file_nameList:
        fn = swc_file.split('/')[-1]

        algorithm = matchFileToAlgorithmName(fn)

        if fn.find("v3dpbd") >-1:
             tmp = fn.split('.v3dpbd')[0] +".v3dpbd"
        else:
             tmp = fn.split('.v3draw')[0] +".v3draw"
        image = tmp.split('sorted_')[-1]  # for sorted_* swc_files
        algorithmList.append(algorithm)
        imageList.append(image)

    df['algorithm'] = pd.Series(algorithmList, index=df.index)
    df['image_file_name'] = pd.Series(imageList, index=df.index)

    allColumns = np.append(np.array(['image_file_name', 'algorithm', 'swc_file']), allColumns, 0)

    df = df[allColumns]

    df.to_csv(feature_csv_file, index=False)

    print 'output all feature csv file to :', feature_csv_file
    return


def cal_bn_features(input_dir,results_feature_csv):
      #results_feature_csv = sorted_dir +'/features_with_tags.csv'
      input_ANO = input_dir+"/input.ano"
      bn.genLinkerFile( input_dir, input_ANO)

      ##batch computing
      feature_file =  input_dir+ "/features.nfb"
      bn.batch_compute (input_ANO,feature_file)

      print "output feature file:"+feature_file
      print "output ano file:"+input_ANO

      generateALLFeatureCSV_gold166(feature_file,results_feature_csv)
      print "output features with tag:"+ results_feature_csv

      return

def cal_blastneuron_distance(results_feature_csv,gold_feature_csv, merged_csv, output_csv, LMEASURE_ONLY = 1):

    df_results = pd.read_csv(results_feature_csv)
    df_gold = pd.read_csv(gold_feature_csv)


    # common set of images
    df_merged = pd.read_csv(merged_csv)
    final_images = np.unique( df_merged.image_file_name)
    print "\n\n Calculating blastneuron ssd scores  for " +str(df_merged.shape[0])+" reconstructions"


    #output_csv =data_DIR + '/normalized_bn_dist.csv'
    feature_cols = [u'num_nodes', u'soma_surface', u'num_stems', u'num_bifurcations', u'num_branches',\
                u'num_of_tips', u'overall_width', u'overall_height', u'overall_depth', u'average_diameter',\
                u'total_length', u'total_surface', u'total_volume', u'max_euclidean_distance',\
                u'max_path_distance', u'max_branch_order', u'average_contraction', u'average fragmentation',\
                u'parent_daughter_ratio', u'bifurcation_angle_local', u'bifurcation_angle_remote', u'moment1',\
                u'moment2', u'moment3', u'moment4', u'moment5', u'moment6', u'moment7', u'moment8', u'moment9',\
                u'moment10', u'moment11', u'moment12', u'moment13', u'avgR']


    # the following three features can easlily have NA values: average_contraction,       average fragmentation ,      parent_daughter_ratio

    selected_allfea_cols = [u'num_nodes', u'num_bifurcations', u'num_branches', u'num_of_tips', u'overall_width',u'overall_height', u'overall_depth',u'total_length',
                     u'bifurcation_angle_remote', u'max_euclidean_distance',u'max_path_distance',
                     u'moment1',u'moment2', u'moment3', u'moment4', u'moment5', u'moment6', u'moment7', u'moment8', u'moment9',u'moment10', u'moment11', u'moment12', u'moment13']
    selected_lmeasure_cols = [ u'num_bifurcations', u'num_branches', u'num_of_tips', u'overall_width',u'overall_height', u'overall_depth',u'total_length',
                     u'bifurcation_angle_remote']
    prefix=""
    if LMEASURE_ONLY:
        selected_cols = selected_lmeasure_cols
        prefix="lm"
    else:
        selected_cols = selected_allfea_cols
        prefix=""

    my_cols1 = [u'image_file_name', u'algorithm', u'swc_file']
    my_cols1.extend(selected_cols)
    df_results_s = df_results[my_cols1]
    del df_results
    my_cols2 = [u'image_file_name']
    my_cols2.extend(selected_cols)
    df_gold_s = df_gold[my_cols2]
    del df_gold

    ### remove na entries
    df_results_s = df_results_s.dropna(axis=0)

    df_results_s[selected_cols]= df_results_s[selected_cols].astype(float)
    df_gold_s[selected_cols]= df_gold_s[selected_cols].astype(float)


    # calcualte std for each image, calculate normalized results
    df_normalized = pd.DataFrame(columns=my_cols1)
    df_gold_normalized = pd.DataFrame(columns=my_cols2)

    df_normalized[[u'image_file_name', u'algorithm', u'swc_file']] = df_results_s[[u'image_file_name', u'algorithm', u'swc_file']]
    df_gold_normalized[u'image_file_name'] = df_gold_s.image_file_name



    for col in selected_cols:
            df_normalized[col] = (df_results_s[col] - df_results_s[col].mean() ) / (df_results_s[col].std() + 0.0000001)
            print df_results_s[col].std()
            print df_results_s[col].mean()
            df_gold_normalized[col] = ( df_gold_s[col] - df_results_s[col].mean() ) / (df_results_s[col].std() + 0.0000001)



    # for i in range(final_images.size):
    #     imageName = final_images[i]
    #
    #     df_image = df_results_s[df_results_s.image_file_name == imageName]
    #
    #     if df_image.shape[0] > 0:  # too few samples
    #         df_image[selected_cols] = df_image[selected_cols].astype(float)  # some moment values are interpreted as strings
    #         df_gold_image = df_gold_s[df_gold_s.image_file_name == imageName]
    #
    #         df_normalized_per_image = pd.DataFrame(columns=my_cols1)
    #         df_gold_normalized_per_image = pd.DataFrame(columns=my_cols2)
    #
    #         df_normalized_per_image[[u'image_file_name', u'algorithm', u'swc_file']] = df_image[[u'image_file_name', u'algorithm', u'swc_file']]
    #         df_gold_normalized_per_image[u'image_file_name'] = df_gold_image.image_file_name
    #
    #         for col in selected_cols:
    #             if (df_image[col].std() == 0.0 ):
    #                 print "warning: std = 0!", col, " ", imageName
    #
    #                 df_normalized_per_image[col] =  (df_image[col] - df_image[col].mean() ) / 1.0
    #                 df_gold_normalized_per_image[col] = ( df_gold_image[col] - df_image[col].mean() ) / 1.0
    #             else:
    #                 df_normalized_per_image[col] = (df_image[col] - df_image[col].mean() ) / (df_image[col].std() + 0.0000001)
    #                 df_gold_normalized_per_image[col] = ( df_gold_image[col] - df_image[col].mean() ) / (df_image[col].std() + 0.0000001)
    #
    #         # append the results for this image
    #         df_normalized = df_normalized.append(df_normalized_per_image, ignore_index=True)
    #         df_gold_normalized= df_gold_normalized.append(df_gold_normalized_per_image, ignore_index=True)

    output_dir = os.path.dirname(output_csv)
    df_gold_normalized.to_csv(output_dir +"/normalized_bnfeatures"+prefix+"_gold.csv",index=False)
    df_normalized.to_csv(output_dir+"/normalized_bnfeatures"+prefix+".csv",index=False)


    mycols= ['image_file_name', 'algorithm', 'swc_file','SSD']
    df_final = pd.DataFrame()
    i = 0
    for rowIdx in range(df_normalized.shape[0]):
        imageName = df_normalized.image_file_name[rowIdx]
        print imageName
        if imageName in final_images:
            gold_df = df_gold_normalized.loc[df_gold_normalized.image_file_name == imageName]

            #normalize restuls with mean and std
            result = df_normalized.iloc[rowIdx]
            print gold_df
            print result
            if gold_df.shape[0] > 0 and result.shape[0] >0 :
                 sumsquare = SSD(result[selected_cols], gold_df[selected_cols])
                 print result[selected_cols]
                 print  "gold: "
                 print gold_df[selected_cols]


                 df_final = df_final.append({'image_file_name': df_normalized.iloc[rowIdx]['image_file_name'],
                                             'algorithm': df_normalized.iloc[rowIdx]['algorithm'],
                                             'swc_file':  df_normalized.iloc[rowIdx]['swc_file'],
                                             'SSD': sumsquare}, ignore_index=True)

                 i=i+1

            else:
                print "warning"
                print imageName
                print result.shape[0]
                print gold_df.shape[0]


    # for rowIdx in range(df_normalized.shape[0]):
    #     imageName = df_normalized.image_file_name[rowIdx]
    #     gold_df = df_gold_normalized.loc[df_gold_normalized.image_file_name == imageName]
    #
    #     #normalize restuls with mean and std
    #     result = df_normalized.iloc[rowIdx]
    #     if gold_df.shape[0] > 0 and result.shape[0] >0 :
    #          sumsquare = SSD(result[selected_cols], gold_df[selected_cols])
    #          ssd_metrics.append(sumsquare)
    #     else:
    #         print "warning"
    #         print result
    #         print gold_df

    #df_final['SSD'] = ssd_metrics
    ## reordering columns
    df_final.to_csv(output_csv,index=False)
    print "output:"+output_csv
    return


def collect_neuron_distance(input_distance_log_file_list, output_distance_csv,lookup_image_id_table_file):
    df_input = pd.read_csv(input_distance_log_file_list, header=None) # txt file contains the list of log files

    df_lookup_table = pd.read_csv(lookup_image_id_table_file)

    df_neuron_distance = pd.DataFrame(columns=('image_file_name','swc_file', 'gold_swc_file',
                                               'algorithm',
                                               'neuron_distance_12','neuron_distance_21',
                                               'neuron_distance_ave','neuron_distance_diff',
                                               'neuron_distance_perc'))
    for i in range(df_input.size):
            logfile_path =  df_input.iloc[i][0]

            image_id = int(logfile_path.split("/")[-3])

            if image_id > df_lookup_table.image_file_name.size:
                  print "error in looking image ids"
            image_file_name = df_lookup_table.image_file_name[image_id-1]


            if path.isfile(logfile_path):
                nd = bn.read_neuron_dist_log(logfile_path)

                algorithm = matchFileToAlgorithmName( logfile_path.split('/')[-1])

                ##return {'input_file1':input_file1, 'input_file2':input_file2,'dist_12':d_12,'dist_21':d_21, 'ave': ave, 'diff': diff, 'perc': perc}
                df_neuron_distance.loc[i] = [image_file_name,nd['input_file1'], nd['input_file2'],
                                             algorithm, nd['dist_12'],nd['dist_21'],nd['ave'],nd['diff'],
                                             nd['perc'] ]
            else:
                print "Warning: no neuron distance log output for "+image_file_name +" :output NAs."
                df_neuron_distance.loc[i]= [image_file_name ,np.nan,np.nan, np.nan, np.nan,np.nan,np.nan, np.nan, np.nan]


    df_neuron_distance['neuron_difference'] = df_neuron_distance['neuron_distance_diff'] *df_neuron_distance['neuron_distance_perc']

    df_neuron_distance.to_csv(output_distance_csv,index=False)
    print "output:"+output_distance_csv
    return






def collect_consensus_distance(input_distance_log_file_list, output_distance_csv,lookup_image_id_table_file):
    print input_distance_log_file_list
    df_input = pd.read_csv(input_distance_log_file_list, header=None) # txt file contains the list of log files

    df_lookup_table = pd.read_csv(lookup_image_id_table_file)

    df_neuron_distance = pd.DataFrame(columns=('image_file_name','consensus_swc_file', 'gold_swc_file',
                                               'algorithm',
                                               'weighted_neuron_distance_12','weighted_neuron_distance_21',
                                               'weighted_neuron_distance_ave','neuron_distance_diff',
                                               'neuron_distance_perc', 'max_distance'))
    for i in range(df_input.size):
            logfile_path =  df_input.iloc[i][0]


            image_id = int(logfile_path.split("/")[-2])

            if image_id > df_lookup_table.image_file_name.size:
                  print "error in looking image ids"
            image_file_name = df_lookup_table.image_file_name[image_id-1]


            if path.isfile(logfile_path):
                nd = bn.read_weighted_neuron_dist_log(logfile_path)

                algorithm ='consensus'
                #{'w_dis_12': wd12, 'w_dis_21': wd21, 'w_ave': w_ave, 'diff': diff, 'perc': perc, 'max_dist':max_dist}
                df_neuron_distance.loc[i] = [image_file_name,nd['input_file1'], nd['input_file2'],
                                             algorithm, nd['w_dis_12'],nd['w_dis_21'],nd['w_ave'],nd['diff'],
                                             nd['perc'],nd['max_dist'] ]
            else:
                print "Warning: no neuron distance log output for "+image_file_name +" :output NAs."
                df_neuron_distance.loc[i]= [image_file_name ,np.nan,np.nan, 'consensus', np.nan, np.nan, np.nan]


    df_neuron_distance['neuron_difference'] = df_neuron_distance['neuron_distance_diff'] *df_neuron_distance['neuron_distance_perc']

    df_neuron_distance.to_csv(output_distance_csv,index=False)
    print "output:"+output_distance_csv
    return



def cal_neuron_dist(input_csv_file,output_csv,overwrite_existing = 1,GEN_QSUB = 0 ):

    df_input = pd.read_csv(input_csv_file)
    df_already_have = pd.DataFrame(columns = df_input.columns)

    # if (not overwrite_existing):
    #     if os.path.isfile(old_output_csv):
    #         # only run new data
    #         df_old = pd.read_csv(old_output_csv)
    #         df_already_have = pd.merge(df_input, df_old, on='swc_file')
    #         print "there are already "+ str(df_already_have['swc_file'].size) +"  swcs calculated"

    output_dir = os.path.dirname(output_csv)
    print " Calculate neuron distances:"
    for i in range(df_input.image_file_name.size):
             print "swc file :" + str(i)
             swc_f = df_input.iloc[i].swc_file
             log_file = df_input.iloc[i].swc_file + ".r.log"
             #if not swc_f in list(df_already_have['swc_file'])
             if overwrite_existing   or   not os.path.isfile(log_file) :
                   bn.run_neuron_dist(swc_f, df_input.iloc[i].gold_swc_file,log_file, GEN_QSUB, output_dir+"/nd")

    #collect results from log files
    df_neuron_distance = pd.DataFrame(columns=('image_file_name','swc_file', 'gold_swc_file',
                                               'algorithm',
                                               'neuron_distance_12','neuron_distance_21',
                                               'neuron_distance_ave','neuron_distance_diff',
                                               'neuron_distance_perc'))


    for i in range(df_input.image_file_name.size):
            tmp = df_input.iloc[i].swc_file
            logfile = tmp + ".r.log"
            if path.isfile(logfile):
                nd = bn.read_neuron_dist_log(logfile)
                df_neuron_distance.loc[i] = [df_input.iloc[i].image_file_name,df_input.iloc[i].swc_file, df_input.iloc[i].gold_swc_file, df_input.iloc[i].algorithm, nd['dist_12'],nd['dist_21'],nd['ave'],nd['diff'],
                                             nd['perc']]
            else:
                print "Warning: no neuron distance log output for "+tmp +" :output NAs."
                df_neuron_distance.loc[i]= [df_input.iloc[i].image_file_name,df_input.iloc[i].swc_file, df_input.iloc[i].gold_swc_file, df_input.iloc[i].algorithm, np.nan, np.nan, np.nan, np.nan, np.nan]


    df_neuron_distance['neuron_difference'] = df_neuron_distance['neuron_distance_diff'] *df_neuron_distance['neuron_distance_perc']

    df_neuron_distance.to_csv(output_csv,index=False)
    print "output:"+output_csv
    return






def cal_neuron_dist_deprecated(input_csv_file,output_csv,overwrite_existing = 1,GEN_QSUB = 0 ):

    df_input = pd.read_csv(input_csv_file)
    df_already_have = pd.DataFrame(columns = df_input.columns)

    # if (not overwrite_existing):
    #     if os.path.isfile(old_output_csv):
    #         # only run new data
    #         df_old = pd.read_csv(old_output_csv)
    #         df_already_have = pd.merge(df_input, df_old, on='swc_file')
    #         print "there are already "+ str(df_already_have['swc_file'].size) +"  swcs calculated"

    output_dir = os.path.dirname(output_csv)
    print " Calculate neuron distances:"
    for i in range(df_input.image_file_name.size):
             print "swc file :" + str(i)
             swc_f = df_input.iloc[i].swc_file
             log_file = df_input.iloc[i].swc_file + ".r.log"
             #if not swc_f in list(df_already_have['swc_file'])
             if overwrite_existing   or   not os.path.isfile(log_file) :
                   bn.run_neuron_dist(swc_f, df_input.iloc[i].gold_swc_file,log_file, GEN_QSUB, output_dir+"/nd")

    #collect results from log files
    df_neuron_distance = pd.DataFrame(columns=('image_file_name','swc_file', 'gold_swc_file', 'algorithm', 'neuron_distance','neuron_distance_diff','neuron_distance_perc'))
    for i in range(df_input.image_file_name.size):
            tmp = df_input.iloc[i].swc_file
            logfile = tmp + ".r.log"
            if path.isfile(logfile):
                nd = bn.read_neuron_dist_log(logfile)
                df_neuron_distance.loc[i] = [df_input.iloc[i].image_file_name,df_input.iloc[i].swc_file, df_input.iloc[i].gold_swc_file, df_input.iloc[i].algorithm, nd['ave'],nd['diff'],nd['perc']]
            else:
                print "Warning: no neuron distance log output for "+tmp +" :output NAs."
                df_neuron_distance.loc[i]= [df_input.iloc[i].image_file_name,df_input.iloc[i].swc_file, df_input.iloc[i].gold_swc_file, df_input.iloc[i].algorithm, np.nan, np.nan, np.nan]


    df_neuron_distance['neuron_difference'] = df_neuron_distance['neuron_distance_diff'] *df_neuron_distance['neuron_distance_perc']

    df_neuron_distance.to_csv(output_csv,index=False)
    print "output:"+output_csv
    return


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


# generate ano files into a single top folder with absoluate paths to swc files reconstructed for one image
def generate_ano_files_dir(my_dir):
    if  not os.path.exists(my_dir+'/ano'):
       os.mkdir(my_dir+'/ano')

    topdirs = glob.glob(os.path.join(my_dir, '*'))
    # print topdirs
    for subdir in topdirs:
         dataset_path = os.path.abspath(subdir)
         bn.genLinkerFile( dataset_path, my_dir+'/ano/'+subdir.split('/')[-1]+'.recons.ano')



def map_image_name(tmp_feature_csv,lookup_image_id_table_file, output_feature_csv):
        df_in = pd.read_csv(tmp_feature_csv)
        df_lookup_table = pd.read_csv(lookup_image_id_table_file)
        df_output = df_in
        for i in range(df_in.image_file_name.size):
            #print df_in.iloc[i].image_file_name
            image_id = int(df_in.iloc[i].image_file_name.split(".")[0])
            if image_id > df_lookup_table.image_file_name.size:
                  print "error in looking image ids"
            image_file_name = df_lookup_table.image_file_name[image_id-1]
            df_output.loc[i,'image_file_name']= image_file_name

        df_output.to_csv(output_feature_csv,index=False)
        return




def recon_table_gen(data_root, lookup_image_id_table_file=None, output_csv_file=None):
        #generate_ano_files_two_level_dir(data_root)  # assuming two-level data archives
        generate_ano_files_dir(data_root)
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
                    swc_file = swc_file[:-1] # get rid of \r or \n
                    fn = swc_file.split('/')[-1]

                    algorithm = matchFileToAlgorithmName(fn,LOOKUP_TABLE_FILE)


                    if fn.find("v3dpbd") >-1:
                        tmp = fn.split('.v3dpbd')[0]
                    else:
                        tmp = fn.split('.v3draw')[0]
                    image = tmp.split('sorted_')[-1]  # for sorted_* swc_files

                    image_file_name = image
                    if os.path.exists(lookup_image_id_table_file):
                        df_lookup_table = pd.read_csv(lookup_image_id_table_file)
                        image_id = int(image)
                        if image_id > df_lookup_table.image_file_name.size:
                            print "error!"
                        image_file_name = df_lookup_table.image_file_name[image_id-1]

                    algorithmList.append(algorithm)
                    imageList.append(image_file_name)
                    swc_list.append(swc_file)

        df_silver['algorithm'] = pd.Series(algorithmList)
        df_silver['image_file_name'] = pd.Series(imageList)
        df_silver['swc_file'] = pd.Series(swc_list)
        if not output_csv_file:
            output_csv_file = data_root + '/recon_table.csv'
        df_silver.to_csv(output_csv_file, index=False)
        return


def merge_gold_silver(GOLD_CSV,SILVER_CSV, output_merged_csv):
    df_gold = pd.read_csv(GOLD_CSV)
    df_silver = pd.read_csv(SILVER_CSV)
    #print df_silver.columns
    #print df_gold.columns

    df_share = pd.merge(df_silver,df_gold,on="image_file_name")


    # df_share = pd.DataFrame([],columns = df_silver.columns)
    #
    # j=0
    # for i in range(df_silver.shape[0]):
    #      if df_silver.image_file_name[i] in df_gold.image_file_name.values:
    #          df_share.loc[j] = df_silver.iloc[i].values
    #          j=j+1

    df_share.to_csv(output_merged_csv, index=False)
    return


def gen_gold_feature_csv(gold_dir,output_gold_csv_file,output_gold_feature_csv):
    #sorted_GMR_57C10_AD_01-1xLwt_attp40_4stop1-m-A02-20111101_2_F3-left_optic_lobe.v3draw.extract_6.v3dpbd.ano_stamp_2015_06_17_12_23.swc
    gold_files = glob.glob(os.path.join(gold_dir, '*.swc'))
    df_gold = pd.DataFrame()
    images = []
    gold_swc_files=[]
    for file in gold_files:
        if file.find(".v3dpbd") > -1:
             image = file.split('.v3dpbd')[0]+".v3dpbd"
        else:
             image = file.split('.v3draw')[0]+".v3draw"

        image = image.split('sorted_')[-1]
        images.append(image)
        gold_swc_files.append(file)

    df_gold['image_file_name'] = pd.Series(images)
    df_gold['gold_swc_file'] = pd.Series(gold_swc_files)
    df_gold.to_csv(output_gold_csv_file, index=False)

    # generate ano file for feature calcuation
    out_sorted_ANO = gold_dir+"/sorted.ano"
    bn.genLinkerFile(gold_dir , out_sorted_ANO)

    out_feature_file =  gold_dir + "/features.nfb"
    bn.batch_compute (out_sorted_ANO,out_feature_file)
    generateALLFeatureCSV_gold166(out_feature_file, output_gold_feature_csv)

    return




def resample_and_sort(data_dir,resampled_dir,sorted_dir, GEN_QSUB = 0,overwrite_sorted = 1,filelist=None):
    failure_file = open(data_dir +"/failurecases_via_size.txt","w")
    i = 0
    subfolder = 1
    if filelist == None:
        filelist= glob.glob(data_dir+"/*/*.swc")
        if len(filelist)<1 :
            filelist= glob.glob(data_dir+"/*.swc")
            subfolder = 0
         
    for input_swc_path in filelist:
         i = i+1
         print "\n\n "
         print i

         #if( os.path.getsize(input_swc_path) < 1024*1024)):
         if not "tmp_cache_img"  in input_swc_path:   # skip the tmp files
              if subfolder >0:
                   swc_fn = "/".join (input_swc_path.split("/")[-2:]) # to keep the subfolder structure
              else:
                   swc_fn = "/".join (input_swc_path.split("/")[-1:])


              sorted_swc_path = sorted_dir+ '/'+swc_fn
              resampled_swc_path = resampled_dir+ '/'+swc_fn
              if not os.path.isfile(resampled_swc_path) or overwrite_sorted:  # if already generated
                 # resample
                 resampled_swc_path = resampled_dir+ '/'+swc_fn
                 print "resample : "+ input_swc_path
                 bn.resample(input_swc_path, resampled_swc_path,3,GEN_QSUB,data_dir+'/qsub/resample')  # generate QSUB scripts

                 # sort
                 sorted_swc_path = sorted_dir+ '/'+swc_fn
                 #bn.sort_swc(preprocessed_swc_path, sorted_swc_path,GEN_QSUB,data_dir+'/qsub/sort')
         #else:
          #  failure_file.write(input_swc_path+" "+str(os.path.getsize(input_swc_path))+"\n ")

    #failure_file.close()
    #print "failure cases with file sizes that are too big or too small are logged at:"+ data_dir+"/failurecases_via_size.txt"
    print "done resampling "#and sorting"
    return



def summerize_running_time(time_csv, algorithm_plugin_match_csv, lookup_image_id_table_file, output_csv):

    # time_csv="/data/mat/xiaoxiaol/data/reconstructions_2015_1214/auto_recons/running_time.csv"
    # output_csv= "/data/mat/xiaoxiaol/data/reconstructions_2015_1214/auto_recons/running_time_algorithm.csv"
    # algorithm_plugin_match_csv ="/data/mat/xiaoxiaol/data/reconstructions_2015_1214/ported_neuron_tracing_spreadsheet.csv"


    df_time = pd.read_csv(time_csv)
    #image,plugin,running_time


    df_check_table = pd.read_csv(algorithm_plugin_match_csv)
    df_lookup_table = pd.read_csv(lookup_image_id_table_file)


    keys = df_check_table['plugin_function']
    values = df_check_table['algorithm']
    #print np.unique(df_check_table.algorithm)
    dictionary = dict(zip(keys, values))


    for i in range(df_time.shape[0]):
        plugin= df_time.iloc[i].plugin
        #print plugin
        alg = dictionary.get(plugin)
        df_time.ix[i,'algorithm'] = alg
        image_id = int(df_time.iloc[i]['image'].split(".")[0])
        image_file_name = df_lookup_table.image_file_name[image_id-1]
        df_time.ix[i,'image_file_name']= image_file_name


    df_time.to_csv(output_csv, index=False)
    return

#####################################################


def main():

    return



if __name__ == "__main__":
    main()



