# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 18:13:25 2015

@author: xiaoxiaol
"""

import pandas as pd
import os



def applyTransformBySpecimenName(swc_path, transform_path,  output_path): 
    
    cmd = ' /local1/xiaoxiaol/work/v3d/v3d_external/bin/vaa3d -x affine_transform_swc -f apply_transform  -i '+ swc_path +  ' -o ' + output_path +' -p ' + transform_path
    
    print cmd
    os.system(cmd)



###########################
transform_DIR="/home/xiaoxiaol/work/data/lims2/nr_june_25_filter_aligned/transforms"  ## where to store the transform.txt
output_DIR="/home/xiaoxiaol/work/data/lims2/nr_june_25_filter_aligned/original" 
origin_data_DIR= "/home/xiaoxiaol/work/data/lims2/nr_june_25_filter/original"
csv_file = "/home/xiaoxiaol/work/data/lims2/nr_june_25_filter_aligned/june_25_alignement_transform.csv"  ## where the transform parameters are, obtained from lims
##########################


if not os.path.exists(output_DIR):
        os.makedirs(output_DIR)

df = pd.read_csv(csv_file)
data_table = df.values
num_samples, num_cols = data_table.shape


for i in range(num_samples):
    orca_path = data_table[i][num_cols-1]
    fn = orca_path.split('/')[-1]
    transform_fn =  transform_DIR+"/%s.txt" % fn.split('.swc')[0]
    output_fn = output_DIR+"/"+fn
    applyTransformBySpecimenName(origin_data_DIR+'/'+fn, transform_fn, output_fn)
     
