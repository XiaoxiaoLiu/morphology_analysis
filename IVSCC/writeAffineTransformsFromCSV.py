# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 18:19:18 2015

@author: xiaoxiaol


[ trv_00 trv_01 trv_02 trv_09
  trv_03 trv_04 trv_05 trv_10
  trv_06 trv_07 trv_08 trv_11 ]
"""
 
 
# table column names
# Index([u'specimen_id', u'specimen_name', u'id', u'tvr_00', u'tvr_01', u'tvr_02', u'tvr_03', u'tvr_04', u'tvr_05', u'tvr_06',
# u'tvr_07', u'tvr_08', u'tvr_09', u'tvr_10', u'tvr_11', u'trv_00', u'trv_01', u'trv_02', u'trv_03', u'trv_04', u'trv_05', 
 #u'trv_06', u'trv_07', u'trv_08', u'trv_09', u'trv_10', u'trv_11', u'metric', u'scale_x', u'scale_y', u'scale_z', u'rotation_x',
# u'rotation_y', u'rotation_z', u'skew_x', u'skew_y', u'skew_z', u'created_at', u'updated_at', u'orca_path'], dtype='object')



import pandas as pd
import os



###########################
data_DIR="/home/xiaoxiaol/work/data/lims2/nr_june_25_filter_aligned/transforms"  ## where to store the transform.txt
csv_file = "/home/xiaoxiaol/work/data/lims2/nr_june_25_filter_aligned/june_25_alignement_transform.csv"  ## where the transform parameters are, obtained from lims
##########################




if not os.path.exists(data_DIR):
        os.makedirs(data_DIR)


df = pd.read_csv(csv_file)
data_table = df.values
num_samples, num_cols = data_table.shape


# transform = np.array([ [ float(a3d.find('tvr-00').text), float(a3d.find('tvr-01').text), float(a3d.find('tvr-02').text), float(a3d.find('tvr-09').text) ],
#                           [ float(a3d.find('tvr-03').text), float(a3d.find('tvr-04').text), float(a3d.find('tvr-05').text), float(a3d.find('tvr-10').text) ],
#                           [ float(a3d.find('tvr-06').text), float(a3d.find('tvr-07').text), float(a3d.find('tvr-08').text), float(a3d.find('tvr-11').text) ],

for i in range(num_samples):
     orca_path = data_table[i][num_cols-1]
     fn= orca_path.split('/')[-1]
     text_file = open(data_DIR+"/%s.txt" % fn.split('.swc')[0]  , "w")

     SCALE = 1000   # from  mm to microns
    
     text_file.write( "%f %f %f %f \n" % (SCALE* df['tvr_00'][i], SCALE* df['tvr_01'][i],SCALE* df['tvr_02'][i],SCALE* df['tvr_09'][i]))

     text_file.write( "%f %f %f %f \n" % (SCALE* df['tvr_03'][i], SCALE* df['tvr_04'][i],SCALE* df['tvr_05'][i],SCALE* df['tvr_10'][i]))
    
     text_file.write( "%f %f %f %f \n" % (SCALE* df['tvr_06'][i], SCALE* df['tvr_07'][i],SCALE* df['tvr_08'][i],SCALE* df['tvr_11'][i]))

     text_file.close()
     
     
     