
__author__ = 'xiaoxiaoliu'
# Module PATHs


import pandas as pd
import platform
import sys
import numpy as np

if (platform.system() == "Linux"):
        WORK_PATH = "/local1/xiaoxiaol/work"
else:
        WORK_PATH = "/Users/xiaoxiaoliu/work"



p =  WORK_PATH + '/src/cell-type-analysis'
sys.path.append(p)

sys.path.append(p + '/utilities')
import morph_nfb_2_csv as nfb

sys.path.append(p + '/blast_neuron')
import blast_neuron_comp as bn







def generateLinkerFileFromCSV(result_dir, csvfile, column_name=None, strip_path=True):
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
    for atype in np.unique(types):
        idxs = np.nonzero(types == atype)[0]
        swc_files = df['swc_file']
        with open(result_dir + '/' + atype + '.ano', 'w') as outf:
            for afile in swc_files[idxs]:
                filename = afile
                if strip_path:
                    filename = afile.split('/')[-1]
                line = 'SWCFILE=' + filename + '\n'
                outf.write(line)
            outf.close()
    return



# Data PATHs
data_dir = WORK_PATH + '/data/lims2/1027_pw_aligned'

output_dir=data_dir+'/clustering_results/ipca/'

df_c = pd.read_csv('~/work/data/lims2/1027_pw_aligned/clustering_results/ipca/ipca_11_clusterid.csv')
df_m = pd.read_csv('~/work/data/lims2/1027_pw_aligned/meta_merged_allFeatures.csv')

df_merged = pd.merge(df_c,df_m,how='inner',on=['specimen_id'])
merged_file = output_dir +"/merged.csv"
df_merged.to_csv(merged_file)

generateLinkerFileFromCSV(output_dir, merged_file,'clusterID',False)