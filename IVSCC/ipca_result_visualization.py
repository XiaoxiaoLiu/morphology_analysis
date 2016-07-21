__author__ = 'xiaoxiaoliu'



import pandas as pd
import numpy as np
import os


def generateLinkerFileFromCSV(result_dir, csvfile, column_name=None, strip_path=True, fpath='.',swc_postfix=None):
    df = pd.read_csv(csvfile, dtype=str)
    if (column_name == None):
        swc_files = df['swc_file']
        with open(result_dir + '/all.ano', 'w') as outf:
            for afile in swc_files:
                filename = afile
                if strip_path:
                    filename = afile.split('/')[-1]
                line = 'SWCFILE=' + filename + '\n'
                outf.write(line)
            outf.close()
        return

    types = df[column_name]

    _, idx = np.unique(types, return_index=True)
    unique_types = types[np.sort(idx)]
    print unique_types


    if swc_postfix is None:
         swc_postfix="_pia.swc"
    count = 1
    for atype in unique_types:
        print atype
        print count
        idxs = np.nonzero(types == atype)[0]
        swc_files = df['swc_file_name']
        with open(result_dir + '/' + str(count) + '.ano', 'w') as outf:
            for afile in swc_files[idxs]:
                filename = afile[:-4]
                if strip_path:
                    filename = afile.split('/')[-1]
                line = 'SWCFILE=' + fpath +"/"+filename + swc_postfix+'\n'
                outf.write(line)
            outf.close()
        count = count+1
    return



# Data PATHs


data_dir = "/home/xiaoxiaol/work/data/lims2/ivscc_0607/CK"

default_all_feature_merged_file = data_dir + '/Hclust.spiny_apical_basal_feature.csv'
output_dir=data_dir

# df_c = pd.read_csv(data_dir + '/../ivscc_0607_new_aligned_with_meta.csv')
# df_m = pd.read_csv(default_all_feature_merged_file)
#
# df_merged = pd.merge(df_c,df_m,how='inner',on=['specimen_name'])
# merged_file = output_dir +"/selected_feature_merged.csv"
# df_merged.to_csv(merged_file)

subfolder = default_all_feature_merged_file .split('/')[-1].split('.csv')[0]
os.system("mkdir  "+output_dir+'/'+subfolder)

#'AP.euc'
generateLinkerFileFromCSV(output_dir+'/'+subfolder, default_all_feature_merged_file,'clusterID.all',strip_path=False,fpath=data_dir + "/../pia_swc", swc_postfix=".swc")
