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



def cal_bn_features(sorted_dir,results_feature_csv):
      #results_feature_csv = sorted_dir +'/features_with_tags.csv'
      sorted_ANO = sorted_dir+"/sorted.ano"
      bn.genLinkerFile( sorted_dir, sorted_ANO)

      ##batch computing
      feature_file =  sorted_dir+ "/features.nfb"
      bn.batch_compute (sorted_ANO,feature_file)

      print "output feature file:"+feature_file
      print "output ano file:"+sorted_ANO

      nfb.generateALLFeatureCSV_gold166(feature_file,results_feature_csv)
      print "output features with tag:"+ results_feature_csv

      return




def SSD(feature_array1, feature_array2):
    diff_v = np.array(feature_array1) - np.array(feature_array2)
    ssd = np.sum(np.abs(diff_v) ** 2)
    return ssd



def cal_blastneuron_distance(results_feature_csv,gold_feature_csv, merged_csv, output_csv):

    df_results = pd.read_csv(results_feature_csv)
    df_gold = pd.read_csv(gold_feature_csv)

    # common set of images
    df_merged = pd.read_csv(merged_csv)
    final_images = np.unique( df_merged.image_file_name)

    #output_csv =data_DIR + '/normalized_bn_dist.csv'
    feature_cols = [u'num_nodes', u'soma_surface', u'num_stems', u'num_bifurcations', u'num_branches',\
                u'num_of_tips', u'overall_width', u'overall_height', u'overall_depth', u'average_diameter',\
                u'total_length', u'total_surface', u'total_volume', u'max_euclidean_distance',\
                u'max_path_distance', u'max_branch_order', u'average_contraction', u'average fragmentation',\
                u'parent_daughter_ratio', u'bifurcation_angle_local', u'bifurcation_angle_remote', u'moment1',\
                u'moment2', u'moment3', u'moment4', u'moment5', u'moment6', u'moment7', u'moment8', u'moment9',\
                u'moment10', u'moment11', u'moment12', u'moment13', u'avgR']

    selected_cols = [u'num_nodes', u'num_bifurcations', u'num_branches', u'num_of_tips', u'overall_width',u'overall_height', u'overall_depth',u'total_length', u'average fragmentation',
            u'bifurcation_angle_local', u'bifurcation_angle_remote', u'moment1',u'moment2', u'moment3', u'moment4', u'moment5', u'moment6', u'moment7', u'moment8', u'moment9',u'moment10', u'moment11', u'moment12', u'moment13']


    my_cols1 = [u'image_file_name', u'algorithm', u'swc_file']
    my_cols1.extend(selected_cols)
    df_results_s = df_results[my_cols1]
    del df_results
    my_cols2 = [u'image_file_name']
    my_cols2.extend(selected_cols)
    df_gold_s = df_gold[my_cols2]
    del df_gold

    # calcualte std for each image, calculate normalized results
    df_normalized = pd.DataFrame(columns=my_cols1)
    df_gold_normalized = pd.DataFrame(columns=my_cols2)

    for i in range(final_images.size):
        imageName = final_images[i]

        df_image = df_results_s[df_results_s.image_file_name == imageName]

        if df_image.shape[0] > 5:  # too few samples
            df_image[selected_cols] = df_image[selected_cols].astype(float)  # some moment values are interpreted as strings
            df_gold_image = df_gold_s[df_gold_s.image_file_name == imageName]

            df_normalized_per_image = pd.DataFrame(columns=my_cols1)
            df_gold_normalized_per_image = pd.DataFrame(columns=my_cols2)

            df_normalized_per_image[[u'image_file_name', u'algorithm', u'swc_file']] = df_image[[u'image_file_name', u'algorithm', u'swc_file']]
            df_gold_normalized_per_image[u'image_file_name'] = df_gold_image.image_file_name

            for col in selected_cols:
                if (df_image[col].std() == 0.0 ):
                    print "warning: std = 0!", col, " ", imageName
                    df_normalized_per_image[col] =  (df_image[col] - df_image[col].median() ) / 1.0
                    df_gold_normalized_per_image[col] = ( df_gold_image[col] - df_image[col].median() ) / 1.0
                else:
                    df_normalized_per_image[col] = (df_image[col] - df_image[col].median() ) / (df_image[col].std() + 0.0000001)
                    df_gold_normalized_per_image[col] = ( df_gold_image[col] - df_image[col].median() ) / (df_image[col].std() + 0.0000001)

            # append the results for this image
            df_normalized = df_normalized.append(df_normalized_per_image, ignore_index=True)
            df_gold_normalized= df_gold_normalized.append(df_gold_normalized_per_image, ignore_index=True)

    output_dir = os.path.dirname(output_csv)
    df_gold_normalized.to_csv(output_dir +"/normalized_bnfeatures_gold.csv")
    df_normalized.to_csv(output_dir+"/normalized_bnfeatures.csv")


    ssd_metrics = []
    for rowIdx in range(df_normalized.shape[0]):
        imageName = df_normalized.image_file_name[rowIdx]
        gold_df = df_gold_normalized.loc[df_gold_normalized.image_file_name == imageName]

        #normalize restuls with mean and std
        result = df_normalized.iloc[rowIdx]
        if gold_df.shape[0] > 0 and result.shape[0] >0 :
             sumsquare = SSD(result[selected_cols], gold_df[selected_cols])
             ssd_metrics.append(sumsquare)
        else:
            print "warning"
            print result
            print gold_df

    df_normalized['SSD'] = ssd_metrics
    ## reordering columns
    df_normalized.to_csv(output_csv,index=False)
    print "output:"+output_csv
    return


def cal_neuron_dist(input_csv_file,output_csv,overwrite_existing = 1, old_output_csv=None):

    df_input = pd.read_csv(input_csv_file)
    df_already_have = pd.DataFrame()

    if (not overwrite_existing):
        if os.path.isfile(old_output_csv):
            # only run new data
            df_old = pd.read_csv(old_output_csv)
            df_already_have = pd.merge(df_input, df_old, on='swc_file')
            print "there are already "+ str(df_already_have['swc_file'].size) +"  swcs calculated"

    print " Calculate neuron distances:"
    for i in range(df_input.image_file_name.size):
             print i
             swc_f = df_input.iloc[i].swc_file
             if not swc_f in list(df_already_have['swc_file']):
                   bn.run_neuron_dist(swc_f, df_input.iloc[i].gold_swc_file,df_input.iloc[i].swc_file + ".r.log", 0, "nd")

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

                    if fn.find("v3dpbd") >-1:
                        tmp = fn.split('v3dpbd_')[-1]
                    else:
                        tmp = fn.split('v3draw_')[-1]
                    algorithm = tmp.split('.')[0]

                    if "app1" in algorithm:   # for patterns like *x245_y234_z234_app1.swc
                          algorithm = "app1"
                    if "app2" in algorithm:
                          algorithm = "app2"
                    if  "fastmarching_spanningtree" in algorithm: # fastmarching_spanningtree is too long
                          algorithm = "spanningtree"

                    if "tubularity_model_S" in algorithm:
                         algorithm = "RegMST"


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
    nfb.generateALLFeatureCSV_gold166(out_feature_file, output_gold_feature_csv)

    return


def resample_and_sort(data_dir,preprocessed_dir,sorted_dir, GEN_QSUB = 0,overwrite_sorted = 1,filelist=None):
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
         print "resample and sort : "+ input_swc_path
         if(( os.path.getsize(input_swc_path) > 1000) and ( os.path.getsize(input_swc_path) < 1024*1024)):
            if not "tmp_cache_img"  in input_swc_path:   # skip the tmp files
              if subfolder >0:
                   swc_fn = "/".join (input_swc_path.split("/")[-2:]) # to keep the subfolder structure
              else:
                   swc_fn = "/".join (input_swc_path.split("/")[-1:])


              sorted_swc_path = sorted_dir+ '/'+swc_fn
              if not os.path.isfile(sorted_swc_path) or overwrite_sorted:  # if already generated
                 # resample
                 preprocessed_swc_path = preprocessed_dir+ '/'+swc_fn
                 bn.resample(input_swc_path, preprocessed_swc_path,3,GEN_QSUB,data_dir+'/qsub/resample')  # generate QSUB scripts

                # sort
                 sorted_swc_path = sorted_dir+ '/'+swc_fn
                 bn.sort_swc(preprocessed_swc_path, sorted_swc_path,GEN_QSUB,data_dir+'/qsub/sort')
         else:
            failure_file.write(input_swc_path+" "+str(os.path.getsize(input_swc_path))+"\n ")

    failure_file.close()
    print "failure cases with file sizes that are too big or too small are logged at:"+ data_dir+"/failurecases_via_size.txt"
    print "done resampling and sorting"
    return


#####################################################


def main():

    return



if __name__ == "__main__":
    main()



