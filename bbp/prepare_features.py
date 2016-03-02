__author__ = 'xiaoxiaoliu'


import pandas as pd


df_axon = pd.read_csv('/Users/xiaoxiaoliu/work/data/BBP/BBP_data_uncurated_SWC/sorted/axon_all_features.csv')
df_dendrite = pd.read_csv('/Users/xiaoxiaoliu/work/data/BBP/BBP_data_uncurated_SWC/sorted/merged_dendrites_tags.csv')

for i in range(len(df_axon)):
     a=df_axon.iloc[i]['file_name']
     fn=a.split('/')[-1]
     df_axon.set_value(i,'file_name',fn.split('.')[0])

df_merged_axon = pd.merge(df_axon, df_dendrite,on='file_name')
df_merged_axon.to_csv('/Users/xiaoxiaoliu/work/data/BBP/BBP_data_uncurated_SWC/sorted/merged_all_features_with_axon_with_tags.csv')
