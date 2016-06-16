# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 14:51:29 2015

@author: xiaoxiaoliu
"""

#merge spreadsheet to generate ano file

import pandas as pd



#df1 = pd.read_csv('/data/mat/xiaoxiaol/data/IVSCC_testing/list.csv')
df1 = pd.read_csv('/data/mat/xiaoxiaol/data/IVSCC_testing/')


#result_dir = '/data/mat/xiaoxiaol/data/IVSCC_testing/linker_files'
result_dir = '/data/mat/xiaoxiaol/data/IVSCC_testing/linker_files_all'


for idx in range(df1.shape[0]):
    specimen_id = df1['specimen_id'][idx]
    with open( result_dir + '/'+str(specimen_id)+'.ano','w') as outf:
        #swc_file = '../swc_pixel/'+ df1['orca_path'][idx].split('/')[-1]
        df1['orca_path'][idx].split('/')[-1]
        swc_file = '../swc_enhanced_1strun/'+T301_min_xy_'+str(specimen_id)+'.jp2.tif_region_APP2.swc'

        mip_file = '../mip_tif/T301_min_xy_'+str(specimen_id)+'.jp2.tif'
        line1='GRAYIMG='+mip_file+'\n'
        line2='SWCFILE='+swc_file+'\n'
        outf.write(line1+line2)
        outf.close()