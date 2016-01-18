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


data_DIR = "/data/mat/xiaoxiaol/data/BBP/raw_converted_swcs"
csv_file = data_DIR+"/meta.csv"


df_bbp = pd.read_csv(csv_file)

generateLinkerFileFromCSV(data_DIR, csv_file, column_name='type',strip_path=False,
                          fpath=data_DIR+'/swc_files')