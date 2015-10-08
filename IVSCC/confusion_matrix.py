__author__ = 'xiaoxiaoliu'

import pandas as pd
from sklearn.metrics import confusion_matrix
import matplotlib.pylab as pl




import numpy as np
import scipy.stats as stats
import seaborn as sns
def cluster_specific_features(df_all, assign_ids, feature_names, output_csv_fn):
    #student t to get cluster specific features

    labels=[]
    clusters = np.unique(assign_ids)
    num_cluster = len(clusters)
    df_pvalues =  pd.DataFrame(index = feature_names, columns = clusters)

    for cluster_id in clusters:

        ids_a = np.nonzero(assign_ids == cluster_id)[0]  # starting from  0
        ids_b = np.nonzero(assign_ids != cluster_id)[0]  # starting from  0
        labels.append(cluster_id + "("+ str(len(ids_a))+")" )
        for feature in feature_names:
            a = df_all.iloc[ids_a][feature]
            b = df_all.iloc[ids_b][feature]

            t_stats,pval = stats.ttest_ind(a,b,equal_var=False)
            df_pvalues.ix[feature,cluster_id] = -np.log10(pval)


    df_pvalues.to_csv(output_csv_fn)

    ### visulaize
    df_pvalues.index.name = "Features"
    df_pvalues.columns.name ="Types"
    d=df_pvalues[df_pvalues.columns].astype(float)
    g = sns.heatmap(data=d,linewidths=0.1)
     #               cmap =sns.color_palette("coolwarm",7, as_cmap=True))

    g.set_xticklabels(labels)
    pl.yticks(rotation=0)
    pl.xticks(rotation=90)
    pl.subplots_adjust(left=0.3, right=0.9, top=0.9, bottom=0.3)
    pl.title('-log10(P value)')
    filename = output_csv_fn + '.png'
    #pl.show()
    pl.savefig(filename, dpi=300)
    pl.close()


    return df_pvalues



data_dir = "/Users/xiaoxiaoliu/work/data/lims2/0923_pw_aligned"
types_id = pd.read_csv(data_dir+'/staci_tag_id.csv')




gl_feature_names = np.array(
        ['num_nodes', 'soma_surface', 'num_stems', 'num_bifurcations', 'num_branches', 'num_of_tips',
         'overall_width', 'overall_height', 'overall_depth', 'average_diameter', 'total_length',
         'total_surface', 'total_volume', 'max_euclidean_distance', 'max_path_distance', 'max_branch_order',
         'average_contraction', 'average fragmentation', 'parent_daughter_ratio', 'bifurcation_angle_local',
         'bifurcation_angle_remote'])

gmi_feature_names = np.array(
    ['moment1', 'moment2', 'moment3', 'moment4', 'moment5', 'moment6', 'moment7', 'moment8',
     'moment9', 'moment10', 'moment11', 'moment12', 'moment13'])  ### removed ave_R
all_feature_names = np.append(gl_feature_names, gmi_feature_names)

df_all = pd.read_csv(data_dir+'/meta_merged_allFeatures.csv')
cluster_specific_features(df_all, types_id['types'], all_feature_names, data_dir+'/fea_pvalues_staci_tags.csv')


#http://stackoverflow.com/questions/5821125/how-to-plot-confusion-matrix-with-string-axis-rather-than-integer-in-python
unique_layers, a = np.unique(df_all.layer_corrected, return_inverse=True)
unique_types, b = np.unique(df_all.types, return_inverse=True)
cm = confusion_matrix(a,b)
print cm.shape
print a.shape
print b.shape

#
# #ax = fig.add_subplot(111)
# pl.matshow(cm)
# #ax.imshow(cm,cmap=pl.cm.jet)
# pl.xticks(range(len(unique_types)),unique_types,rotation=90)
# pl.yticks(range(len(unique_layers)),unique_layers)
# pl.colorbar()
# pl.show()



# cluster1_ids =  pd.read_csv(data_dir+'/clustering_results/ward_all_ol_clipped/cluster_id.csv')
# cluster2_ids =  pd.read_csv(data_dir+'/clustering_results/ward_mrmr_ol_clipped/cluster_id.csv')
#
# # matching by specimen_name
# df_merged = pd.merge(cluster1_ids, cluster2_ids, how='inner', on=['specimen_name'])
#
# # order
# #sorted_cluster2_ids = cluster2_ids.sort(['specimen_name'])
#
#
# cm = confusion_matrix(df_merged.cluster_id_x, df_merged.cluster_id_y)
# pl.matshow(cm)
# pl.colorbar()
# pl.show()



