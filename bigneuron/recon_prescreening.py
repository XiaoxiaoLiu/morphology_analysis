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


import blast_neuron.blast_neuron_comp as bn
import glob




def cal_neuron_dist(merged_csv_file):
    df_merged = pd.read_csv(merged_csv_file)
    for i in range(df_merged.image_file_name.size):
        bn.run_neuron_dist(df_merged.iloc[i].swc_file, df_merged.iloc[i].gold_swc_file,df_merged.iloc[i].swc_file + ".r.log", 0, "nd")
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
                    if  "spanningtree" in algorithm: # fastmarching_spanningtree is too long
                          algorithm = "spanningtree"

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
                            print "error! line 87 at recon_prescreening.py"
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
    print df_silver.columns
    print df_gold.columns

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


def gen_gold_csv(gold_dir,output_gold_csv_file):
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
    return


def resample_and_sort(data_dir,preprocessed_dir,sorted_dir, GEN_QSUB = 0):
    failure_file = open(data_dir +"/failurecases_via_size.txt","w")
    i = 0
    print data_dir
    for input_swc_path in glob.glob(data_dir+"/*/*.swc"):
         i = i+1
         print "\n\n "
         print i
         print "resample and sort : "+ input_swc_path
         if(( os.path.getsize(input_swc_path) > 1000) and ( os.path.getsize(input_swc_path) < 1024*1024/2)):
            if not "tmp_cache_img"  in input_swc_path:   # skip the tmp files
              swc_fn = "/".join (input_swc_path.split("/")[-2:]) # to keep the subfolder structure

              # resample
              preprocessed_swc_path = preprocessed_dir+ '/'+swc_fn
              bn.resample(input_swc_path, preprocessed_swc_path,3,GEN_QSUB,data_dir+'/qsub/resample')  # generate QSUB scripts

              # sort
              sorted_swc_path = sorted_dir+ '/'+swc_fn
              bn.sort_swc(preprocessed_swc_path, sorted_swc_path,GEN_QSUB,data_dir+'/qsub/sort')
         else:
            failure_file.write(input_swc_path+" "+os.path.getsize(input_swc_path)+"\n ")

    failure_file.close()
    print "failure cases with file sizes that are too big or too small are logged at:"+ data_dir+"/failurecases_via_size.txt"
    print "done resampling and sorting"
    return


#####################################################


def main():
    data_DIR ="/data/mat/xiaoxiaol/data/20151130_rhea_reconstructions_for_allen300_silver_set"
    original_dir = data_DIR +"/auto_recons"
    preprocessed_dir = data_DIR+ "/resampled"
    sorted_dir = data_DIR +"/sorted"


    ###### generate gold standard 166 table
    gold_dir = "/data/mat/xiaoxiaol/data/gold166/checked_final_swcs/preprocessed"
    GOLD_CSV = "/data/mat/xiaoxiaol/data/gold166/gold.csv"
    gen_gold_csv(gold_dir, GOLD_CSV)



    ######  resample and sort
    resample_and_sort(original_dir,preprocessed_dir,sorted_dir)


    ###### genearte sliver data table
    SILVER_CSV = data_DIR+'/recon_table.csv'
    lookup_image_id_table_file = data_DIR +"/image_name_lookup_table.csv"
    recon_table_gen(sorted_dir,lookup_image_id_table_file,SILVER_CSV)


    #####merge to get the common set between gold and silver
    merged_csv_file = data_DIR+'/shared.csv'
    merge_gold_silver(GOLD_CSV,SILVER_CSV,merged_csv_file)

    cal_neuron_dist(merged_csv_file)
    return



if __name__ == "__main__":
    main()



