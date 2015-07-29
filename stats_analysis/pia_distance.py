# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 00:16:10 2015

@author: xiaoxiaoliu
"""
#generate  pia distance

# read from swc files and find out the root node and grab the z cooordinates

data_DIR = '/Users/xiaoxiaoliu/work/data/lims2/nr_june_25_filter_aligned/original'
fn='/Users/xiaoxiaoliu/work/data/lims2/nr_june_25_filter_aligned/allFeatures_withid.csv'
df=pd.read_csv(fn)

i=0
dis=[]
for m_path in  df.orca_path:
    swc_fn = data_DIR + '/'+m_path.split('/')[-1]
    with open(swc_fn, 'r') as FILE:
        for line in FILE: ##n,type,x,y,z,radius,parent
          if (  (line.strip()).split(' ')[-1] == '-1'):
              #root identified
              dis.append(line.split(' ')[4])
              print dis
              break
    i= i+1;

df['pia_distance'] = dis
df.to_csv('/Users/xiaoxiaoliu/work/data/lims2/nr_june_25_filter_aligned/allFeatures_withid_piaDistance.csv')