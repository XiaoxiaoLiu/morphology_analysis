import pandas as pd
import platform


if (platform.system() == "Linux"):
    WORK_PATH = "/local1/xiaoxiaol/work"
else:
    WORK_PATH = "/Users/xiaoxiaoliu/work"
###############################################################################
#data_DIR = '/data/mat/xiaoxiaol/data/lims2/0903_filtered_ephys_qc'
data_DIR = WORK_PATH +'/data/lims2/0923_pw_aligned'
# /original stores the downloaded swc files
original_dir = data_DIR + "/pw_aligned"
preprocessed_dir = data_DIR + "/preprocessed"

###############################################################################


#gl_feature_names = np.array(
        # ['total_length', 'soma_surface', 'num_stems', 'num_bifurcations', 'num_branches', 'num_of_tips',
        #  'overall_width', 'overall_height', 'overall_depth', 'average_diameter', 'num_nodes',
        #  'total_surface', 'total_volume', 'max_euclidean_distance', 'max_path_distance', 'max_branch_order',
        #  'average_contraction', 'average fragmentation', 'parent_daughter_ratio', 'bifurcation_angle_local',
        #  'bifurcation_angle_remote','height_width_ratio','average_branch_length','length_surface_ratio'])
        #

fnames = [u'Unnamed: 0', u'specimen_id', u'specimen_name', u'dendrite_type',
       u'cre_line', u'layer', u'swc_file', u'num_nodes', u'soma_surface',
       u'num_stems', u'num_bifurcations', u'num_branches', u'num_of_tips',
       u'overall_depth', u'overall_width', u'overall_height',
       u'average_diameter', u'total_length', u'total_surface', u'total_volume',
       u'max_euclidean_distance', u'max_path_distance', u'max_branch_order',
       u'average_contraction', u'average fragmentation',
       u'parent_daughter_ratio', u'bifurcation_angle_local',
       u'bifurcation_angle_remote', u'moment1', u'moment2', u'moment3',
       u'moment4', u'moment5', u'moment6', u'moment7', u'moment8', u'moment9',
       u'moment10', u'moment11', u'moment12', u'moment13', u'avgR']

df_f = pd.read_csv(preprocessed_dir + '/features_with_db_tags.csv')

df_f.columns=fnames

df_f['height_width_ratio'] = df_f['overall_height']/df_f['overall_width']
df_f['average_branch_length']=df_f['total_length']/df_f['num_branches']
df_f ['length_surface_ratio'] = df_f ['total_length']/df_f['total_surface']

df_f.to_csv(preprocessed_dir + '/features_with_db_tags_added.csv')

