__author__ = 'xiaoxiaoliu'
__author__ = 'xiaoxiaoliu'


import pandas as pd
import numpy as np
import os
import sys

import matplotlib.pyplot as plt
WORK_PATH = "/Users/xiaoxiaoliu/work"
p =  WORK_PATH + '/src/morphology_analysis'
sys.path.append(p)
import blast_neuron.blast_neuron_comp as bn
import glob

import libtiff3d
from scipy.stats.stats import pearsonr


def read_tif(tif_file):
    t = libtiff3d.TIFF3D.open(tif_file)
    im = t.read_image()   # z_dim, y_dim, x_dim
    return im

gold_image_dir = WORK_PATH+"/data/gold79/origin_data"
gold_images =  glob.glob(gold_image_dir+'/*/*/*.tif')

votemaps_dir = WORK_PATH+"/data/20151030_rhea_reconstructions_for_allen300_silver_set/votemaps"
vote_images  = glob.glob(votemaps_dir+'/*.tif')
images=[]
gold_image_files=[]
vote_image_files=[]
pearsons=[]
dim_x=[]
dim_y=[]
dim_z=[]
pval=[]
for vote_image_file in vote_images:
    Iv= read_tif(vote_image_file)
    image_name = vote_image_file.split('/')[-1]
    image_name = image_name.split('recons.ano')[0]
    image_name = ".".join(image_name.split('.')[1:])

    for gold_image_file in gold_images:
        if image_name in gold_image_file:
             Ig=read_tif(gold_image_file)

             siz_x = min(Ig.shape[0], Iv.shape[0])
             siz_y = min(Ig.shape[1], Iv.shape[1])
             siz_z = min(Ig.shape[2], Iv.shape[2])

             Ig_match = Ig[0:siz_x, 0:siz_y,0:siz_z]
             Iv_match = Iv[0:siz_x, 0:siz_y,0:siz_z]

             if (Ig.shape[0] < Iv.shape[0]) or (Ig.shape[1] < Iv.shape[1]) or (Ig.shape[2] < Iv.shape[2]):
                print "wrong dim"


             print Ig_match.shape
             print Iv_match.shape

             # plt.figure()
             # plt.imshow(Iv_match[Iv_match.shape[0]/2,:,:])
             # plt.show()
             # plt.figure()
             # plt.imshow(Ig_match[Ig_match.shape[0]/2,:,:])
             # plt.show()


             pr = pearsonr(Iv_match.flatten(), Ig_match.flatten())

             pvalue = 0 #? place holder for pvalue

             pearsons.append(pr[0])
             pval.append(pvalue)
             dim_x.append(Ig.shape[2])
             dim_y.append(Ig.shape[1])
             dim_z.append(Ig.shape[0])
             images.append(image_name)
             gold_image_files.append(gold_image_file)
             vote_image_files.append(vote_image_file)



df=pd.DataFrame()
df['image'] = pd.Series(images)
df['gold_image_file'] = pd.Series(gold_image_files)
df['votemap_image_file'] = pd.Series(vote_image_files)
df['perasonr'] = pd.Series(pearsons)
df['pval']=pd.Series(pval)

df['dim_x'] = pd.Series(dim_x)
df['dim_y'] = pd.Series(dim_y)
df['dim_z'] = pd.Series(dim_z)
df.to_csv(votemaps_dir+"/pearsonr.csv", index=False)

