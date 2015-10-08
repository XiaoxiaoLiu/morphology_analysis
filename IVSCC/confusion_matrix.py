__author__ = 'xiaoxiaoliu'

import pandas as pd
from sklearn.metrics import confusion_matrix
import matplotlib.pylab as pl

data_dir = "/Users/xiaoxiaoliu/work/data/lims2/0923_pw_aligned"

cluster1_ids =  pd.read_csv(data_dir+'/clustering_results/ward_all_ol_clipped/cluster_id.csv')
cluster2_ids =  pd.read_csv(data_dir+'/clustering_results/ward_mrmr_ol_clipped/cluster_id.csv')

# matching by specimen_name
df_merged = pd.merge(cluster1_ids, cluster2_ids, how='inner', on=['specimen_name'])

# order
#sorted_cluster2_ids = cluster2_ids.sort(['specimen_name'])


cm = confusion_matrix(df_merged.cluster_id_x, df_merged.cluster_id_y)
pl.matshow(cm)
pl.colorbar()
pl.show()
