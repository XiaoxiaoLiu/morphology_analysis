__author__ = 'xiaoxiaol'


import pandas as pd
import numpy as np
import os
import sys


WORK_PATH = "/Users/xiaoxiaoliu/work"
p =  WORK_PATH + '/src/morphology_analysis'
sys.path.append(p)
import blast_neuron.blast_neuron_comp as bn
import glob



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
#####################################################

data_DIR = WORK_PATH+"/data/20151030_rhea_reconstructions_for_allen300_silver_set"
original_dir = data_DIR + "/auto_recons"
preprocessed_dir = data_DIR +"/resampled"
sorted_dir = data_DIR +"/sorted"


#one-time operation
GEN_GOLD_CSV = 0
gold_dir = WORK_PATH+"/data/gold79/sorted"
if GEN_GOLD_CSV:
    #sorted_GMR_57C10_AD_01-1xLwt_attp40_4stop1-m-A02-20111101_2_F3-left_optic_lobe.v3draw.extract_6.v3dpbd.ano_stamp_2015_06_17_12_23.swc
    gold_files = glob.glob(os.path.join(gold_dir, '*.swc'))
    df_gold = pd.DataFrame()
    images = []
    gold_swc_files=[]
    for file in gold_files:
        image = file.split('.v3dpbd')[0]
        image = image.split('sorted_')[-1]
        images.append(image)
        gold_swc_files.append(file)

    df_gold['image'] = pd.Series(images)
    df_gold['gold_swc_file'] = pd.Series(gold_swc_files)
    df_gold.to_csv(gold_dir + '/gold.csv', index=False)




######################
###### tmp
#sorted_dir = original_dir
######

sorted_dir = WORK_PATH+"/data/20151030_rhea_reconstructions_for_allen300_silver_set/79/sorted"
recon_table_gen(sorted_dir)

#merge to get the common set between gold and silver



SELECT_SWC_WITH_GOLD = 0
if  SELECT_SWC_WITH_GOLD:
    df_gold = pd.read_csv(WORK_PATH+"/data/gold79/sorted/gold.csv")
    df_silver=pd.read_csv((original_dir+'/recon_table.csv'))

    df_share = pd.DataFrame([],columns = df_silver.columns)

    j=0
    for i in range(df_silver.shape[0]):
         if df_silver.image[i] in df_gold.image.values:
             df_share.loc[j] = df_silver.iloc[i].values
             j=j+1

    df_share.to_csv(original_dir+"/shared.csv", index=False)

# compare consensus to gold standard
# generate the gold standard




