__author__ = 'xiaoxiaoliu'




import pandas as pd

import numpy as np





def generateLinkerFileFromCSV(result_dir, csvfile, column_name=None, strip_path=True, fpath='.'):
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
    for atype in np.unique(types):
        idxs = np.nonzero(types == atype)[0]
        swc_files = df['swc_file_name']
        with open(result_dir + '/' + atype + '.ano', 'w') as outf:
            for afile in swc_files[idxs]:
                filename = afile
                if strip_path:
                    filename = afile.split('/')[-1]
                line = 'SWCFILE=' + fpath +"/"+filename + '\n'
                outf.write(line)
            outf.close()
    return



# Data PATHs


data_dir = "/data/mat/xiaoxiaol/data/lims2/pw_aligned_1223"

default_all_feature_merged_file = data_dir + '/0107_new_features.csv'
output_dir=data_dir+'/clustering_results/ipca/'

df_c = pd.read_csv(output_dir + '/Final.Ephys.seed0.Regular.P0.05.GrpbyPC.FSbyttest.IncPC.PCA.0.056.csv')
df_m = pd.read_csv(default_all_feature_merged_file)

df_merged = pd.merge(df_c,df_m,how='inner',on=['specimen_id'])
merged_file = output_dir +"/merged.csv"
df_merged.to_csv(merged_file)

generateLinkerFileFromCSV(output_dir, merged_file,'clusterID',strip_path=False,fpath=data_dir + "/keith_swc_22dec")