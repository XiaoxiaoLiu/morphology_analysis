__author__ = 'xiaoxiaoliu'




import pandas as pd

import numpy as np





def generateLinkerFileFromCSV(result_dir, csvfile, column_name=None, strip_path=True, fpath='.',swc_postfix=None):
    df = pd.read_csv(csvfile)
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
    print np.unique(types)
    if swc_postfix is None:
         swc_postfix=".swc"
    for atype in np.unique(types):
        idxs = np.nonzero(types == atype)[0]
        swc_files = df['swc_file_name']
        with open(result_dir + '/' + str(atype) + '.ano', 'w') as outf:
            for afile in swc_files[idxs]:
                filename = afile[:-4]
                if strip_path:
                    filename = afile.split('/')[-1]
                line = 'SWCFILE=' + fpath +"/"+filename + swc_postfix+'\n'
                outf.write(line)
            outf.close()
    return



# Data PATHs


data_dir = "/home/xiaoxiaol/work/data/lims2/ivscc_0519/CK"

default_all_feature_merged_file = data_dir + '/Final.specimen.selected.all.ss.06032016.csv'
output_dir=data_dir

df_c = pd.read_csv(data_dir + '/../ivscc_0519_features_with_meta.csv')
df_m = pd.read_csv(default_all_feature_merged_file)

df_merged = pd.merge(df_c,df_m,how='inner',on=['specimen_name'])
merged_file = output_dir +"/selected_feature_merged.csv"
df_merged.to_csv(merged_file)

generateLinkerFileFromCSV(output_dir+'/all', merged_file,'clusterID.all',strip_path=False,fpath=data_dir + "/../PCA_aligned",swc_postfix="_pia.swc")
generateLinkerFileFromCSV(output_dir+'/selected', merged_file,'clusterID',strip_path=False,fpath=data_dir + "/../PCA_aligned",swc_postfix="_pia.swc")