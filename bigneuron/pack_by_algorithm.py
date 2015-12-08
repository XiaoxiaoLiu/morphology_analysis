__author__ = 'xiaoxiaol'
import pandas as pd
import os
import numpy as np



data_dir = "/data/mat/xiaoxiaol/data/reconstructions_2015_1207/auto_recons"
df_r=pd.read_csv(data_dir+'/list.csv')

lookup_image_id_table_file = data_dir +"/../image_name_lookup_table.csv"
df_lookup_table = pd.read_csv(lookup_image_id_table_file)
out_folder = "org_by_alg"


for swc_file  in df_r.swc_file:

    fn = swc_file.split('/')[-1]

    if fn.find("v3dpbd_") >-1:
        tmp = fn.split('v3dpbd_')[-1]
    else:
        tmp = fn.split('v3draw_')[-1]
    algorithm = tmp.split('.')[0]

    if "app1" in algorithm:   # for patterns like *x245_y234_z234_app1.swc
          algorithm = "app1"
    if "app2" in algorithm:
          algorithm = "app2"
    if  "fastmarching_spanningtree" in algorithm: # fastmarching_spanningtree is too long
          algorithm = "spanningtree"
    if "tubularity_model_S" in algorithm:
         algorithm = "RegMST"





    image_id = int(fn.split(".v3d")[0])
    if image_id > df_lookup_table.image_file_name.size:
          print "error in looking image ids"
    image_file_name = df_lookup_table.image_file_name[image_id-1]
    #print image_file_name

    outfilename = image_file_name+".swc"
    outdir= data_dir+"/"+out_folder+"/"+algorithm
    os.system("mkdir -p "+outdir)
    os.system("cp "+data_dir+"/"+swc_file +" "+ outdir+"/"+outfilename)

    #print data_dir+"/"+swc_file
    #print data_dir+"/pack/"+algorithm+"/"+outfilename


