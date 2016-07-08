__author__ = 'xiaoxiaol'
import matplotlib.pyplot as plt
import seaborn as sb
import os
import os.path as path
import numpy as np
import pandas as pd

data_DIR = "/data/mat/xiaoxiaol/data/big_neuron/silver/0401_gold163_all_soma_sort"
#df_metrics = pd.read_csv(data_DIR+'/image_profiling_p0.05_original_gs.csv')
df_metrics = pd.read_csv(data_DIR+'/radius_estimation_profiling-strict.csv')

#
df_meta = pd.read_csv(data_DIR+'/image_name_lookup_table_with_limited_meta.csv')
df_data= pd.merge(df_metrics,df_meta, on='image_id' )
df_data['image_id']=df_data['image_id'].astype('int')
df_data["meta"] =df_data["species"]+ "-"+df_data["lab"]+ " ("+df_data["image_id"].map(str)+")"
#df_data.to_csv(data_DIR+'/image_profiling_p0.05_original_gs_with_meta.csv', index=False)
df_data.to_csv(data_DIR+'/image_profiling_p0.05_reestimated_radius_gs_with_meta.csv', index=False)
#df_data = pd.read_csv(data_DIR+'/image_profiling_p0.05_original_gs_with_meta.csv')

df_data.sort(['SNR'], ascending=[1], inplace=True)
df_data.loc[df_data['CNR'] >100, 'CNR'] = 100
#sort_by_cnr = np.argsort(df_metrics['CNR'])

sb.set_context("poster")


f, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(30, 15), sharex=True)

sb.barplot(data=df_data,x='image_id', y='CNR', order = df_data.image_id,    ax=ax1)
#   plt.xticks(range(len(df_data)), df_data['image_id'], rotation='vertical')

sb.barplot(data=df_data,x='image_id', y='SNR', order = df_data.image_id,    ax=ax2)
#   plt.xticks(range(len(df_data)), df_data['image_id'], rotation='vertical')

sb.barplot(data=df_data,x='image_id', y='mean_tubularity', order = df_data.image_id,    ax=ax3)
#plt.xticks(range(len(df_data)), df_data['meta'], rotation='vertical')

ax1.set_xlabel('')
ax2.set_xlabel('')
#sb.barplot(data=df_data,x='image_id', y='dynamic_range', order = df_data.image_id,    ax=ax4)
plt.xticks(range(len(df_data)), df_data['meta'], rotation='vertical')

sb.despine(bottom=True)
#plt.setp(f.axes, yticks=[])
plt.tight_layout(h_pad=0)
plt.xlabel('images')

#plt.plot(range(len(df_metrics)),df_metrics['CNR'])
#plt.savefig(data_DIR+'/image_profile_metrics_meta.png')
plt.savefig(data_DIR+'/image_profile_metrics_meta_radius_reestimated.png')
plt.show()
#plt.close()
