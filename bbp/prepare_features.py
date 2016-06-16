
import pandas as pd
import numpy as np
import IVSCC.Clustering_morphology_on_detailed_features as mc

from sklearn.metrics import confusion_matrix
import matplotlib.pylab as plt


work_path = "/data/mat/xiaoxiaol"




## old features
# df_axon = pd.read_csv(work_path+'/data/BBP/BBP_data_uncurated_SWC/sorted/axon_all_features.csv')
# df_dendrite = pd.read_csv(work_path+'/data/BBP/BBP_data_uncurated_SWC/sorted/merged_dendrites_tags.csv')
#
# for i in range(len(df_axon)):
#      a=df_axon.iloc[i]['file_name']
#      fn=a.split('/')[-1]
#      df_axon.set_value(i,'file_name',fn.split('.')[0])
#
# df_merged_axon = pd.merge(df_axon, df_dendrite,on='file_name')
# df_merged_axon.to_csv('/Users/xiaoxiaoliu/work/data/BBP/BBP_data_uncurated_SWC/sorted/merged_all_features_with_axon_with_tags.csv')


def plot_confusion_matrix(cm, xlabel, ylabel, xnames, ynames,  title='Confusion matrix', cmap=plt.cm.Blues):
    plt.grid(False)
    plt.imshow(cm, interpolation = 'none',cmap=cmap)
    plt.title(title)
    plt.colorbar()
    tick_marksx = np.arange(len(xnames))
    tick_marksy = np.arange(len(ynames))
    plt.xticks(tick_marksx, xnames)
    plt.yticks(tick_marksy, ynames)
    #plt.tight_layout()
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)





data_DIR = work_path+'/data/BBP/BBP_data_uncurated_SWC'
df_cloud_2d = pd.read_csv(data_DIR +'/bbp_features_0411.cloud.xy.csv')
df_meta = pd.read_csv(data_DIR+'/rich_meta_type.csv')

if 0:
    csv_file = data_DIR+"/rich_meta_type.csv"
    df_bbp = pd.read_csv(csv_file)
    for i in range (df_bbp.shape[0]):
        df_bbp.ix[i,'layer'] = df_bbp.iloc[i]['m-type'].split('_')[0]

    df_bbp.to_csv(data_DIR+'/final_meta.csv')
    #generateLinkerFileFromCSV(data_DIR, csv_file, column_name='m-type',strip_path=False,
     #                         fpath=data_DIR+'/sorted')
## process feature csv
if 0:
     df_f =df_cloud_2d

     #df_f= df_f.drop(['specimen_name', 'specimen_id', 'dendrite_type', 'region_info','ignore'],axis=1)

     for i in range (df_f.shape[0]):
        tmp = df_f.iloc[i]['filename'].split('/')[-1]
        df_f.ix[i,'swc_file_name'] = (tmp.split('.')[0]).upper() +'.swc'
        df_f.ix[i,'specimen_name'] = (tmp.split('.')[0]).upper() +'.swc'


     #df_f.to_csv(data_DIR+'/bbp_features_0411.cloud.xy_processed.csv', index=False)
     print df_f.shape
     print df_bbp.shape
     df_2=pd.merge(df_f, df_bbp, on="swc_file_name")
     print df_2.shape
     df_2.to_csv(data_DIR+'/bbp_features_0411.cloud.xy.with_meta.csv', index=False)

     dfg=df_2.groupby(['dendrite_type'])


     df_i =dfg.get_group('inhibitory')
     df_i.to_csv(data_DIR+'/bbp_inhibitory_cloud_xy_features.csv', index=False)

if 1:
     #run hierachical clustering
     cloud_xy_features = [ u'cloud_first_compartment_moment_x',
                           u'cloud_first_compartment_moment_y',
                           u'cloud_height',
                           u'cloud_second_compartment_moment_x',
                           u'cloud_second_compartment_moment_y',
                           u'cloud_width']
     df_i = pd.read_csv(data_DIR+'/bbp_inhibitory_cloud_xy_features.csv')


     linkage, df_zscore= mc.ward_cluster(df_i, cloud_xy_features, 38,
                               data_DIR + '/cloud2D_inh_WARD', None, 0, 'bbp')
     mc.silhouette_clusternumber(linkage, df_zscore,  data_DIR + '/cloud2D_inh_WARD')
     #gen confusion matrix

    #assign_ids=np.array(df_i['m-type'])



#http://stackoverflow.com/questions/5821125/how-to-plot-confusion-matrix-with-string-axis-rather-than-integer-in-python

df_i = pd.read_csv(data_DIR+'/bbp_inhibitory_cloud_xy_features.csv')


df_cluster_ids =  pd.read_csv( data_DIR + '/cloud2D_inh_WARD/cluster_id.csv')
cluster_ids = np.unique(df_cluster_ids['cluster_id'])
print cluster_ids
metric ='m-type'
shapes = np.unique(df_i[metric])
shapes_lut = dict(zip(shapes, range(1,len(shapes)+1)))  # map creline type to color
shape_id = df_i[metric].map(shapes_lut)
print shapes

cm = confusion_matrix(shape_id,df_cluster_ids['cluster_id'])
# print cm.shape
# #np.set_printoptions(precision=2)
# print('Confusion matrix, without normalization')
# print(cm)
# plt.figure()
# plot_confusion_matrix(cm,'cluster id','shape types',cluster_ids, shapes)
#
# plt.show()


# Normalize the confusion matrix by row (i.e by the number of samples
# in each class)
cm_normalized = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
print('Normalized confusion matrix')
print(cm_normalized)
plt.figure()
plot_confusion_matrix(cm_normalized, 'cluster id','shape types',cluster_ids, shapes, "Normaized confusion matrix")

plt.show()






