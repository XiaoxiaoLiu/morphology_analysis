from sklearn.cluster import AffinityPropagation
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs

##############################################################################


########################################## data dir
data_DIR = WORK_PATH + "/data/lims2/0923_pw_aligned"
#########################################################
#data_linker_file = data_DIR + '/original/mylinker.ano'
#preprocessed_data_linker_file = data_DIR + '/preprocessed/mylinker.ano'
FEATURE_FILE = data_DIR + '/preprocessed/prep_features.nfb'

gl_feature_names = np.array(
    ['num_nodes', 'soma_surface', 'num_stems', 'num_bifurcations', 'num_branches', 'num_of_tips',
     'overall_width', 'overall_height', 'overall_depth', 'average_diameter', 'total_length',
     'total_surface', 'total_volume', 'max_euclidean_distance', 'max_path_distance', 'max_branch_order',
     'average_contraction', 'average fragmentation', 'parent_daughter_ratio', 'bifurcation_angle_local',
     'bifurcation_angle_remote'])

gmi_feature_names = np.array(['moment1', 'moment2', 'moment3', 'moment4', 'moment5', 'moment6', 'moment7', 'moment8',
                              'moment9', 'moment10', 'moment11', 'moment12', 'moment13'])  ### removed ave_R

selected_features = ['max_euclidean_distance', 'num_stems', 'num_bifurcations', 'average_contraction',
                     'parent_daughter_ratio']

all_feature_names = np.append(gl_feature_names, gmi_feature_names)



all_feature_merged_file = data_DIR + '/preprocessed/features_with_db_tags.csv'
merged = pd.read_csv(all_feature_merged_file)

merged[all_feature_names]= merged[all_feature_names].astype(float)



df_zscores = zscore_features(merged, all_feature_names, None, 1)

X=df_zscores.as_matrix()
##############################################################################
# Compute Affinity Propagation
af = AffinityPropagation().fit(X)
cluster_centers_indices = af.cluster_centers_indices_
labels = af.labels_

n_clusters_ = len(cluster_centers_indices)

print('Estimated number of clusters: %d' % n_clusters_)
#print("Homogeneity: %0.3f" % metrics.homogeneity_score(labels_true, labels))
#print("Completeness: %0.3f" % metrics.completeness_score(labels_true, labels))
#print("V-measure: %0.3f" % metrics.v_measure_score(labels_true, labels))
#print("Adjusted Rand Index: %0.3f"
#      % metrics.adjusted_rand_score(labels_true, labels))
#print("Adjusted Mutual Information: %0.3f"
#      % metrics.adjusted_mutual_info_score(labels_true, labels))
#print("Silhouette Coefficient: %0.3f"
#      % metrics.silhouette_score(X, labels, metric='sqeuclidean'))
