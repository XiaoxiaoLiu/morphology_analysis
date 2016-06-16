import numpy as np
import pandas as pd


# program path on this machine
#===================================================================
WORK_PATH="/Users/xiaoxiaoliu/work"
#WORK_PATH="/home/xiaoxiaol/work"



# data dir
data_DIR= "/Volumes/mat/xiaoxiaol/data/lims2/nr_june_25_filter_aligned/apical"


FEATURE_FILE = data_DIR + '/preprocessed/prep_features.nfb'
gl_feature_names= np.array(['num_nodes', 'soma_surface', 'num_stems','num_bifurcations', 'num_branches', 'num_of_tips',  'overall_width', 'overall_height',  'overall_depth', 'average_diameter',    'total_length', 'total_surface', 'total_volume', 'max_euclidean_distance',       'max_path_distance', 'max_branch_order',  'average_contraction', 'average fragmentation', 'parent_daughter_ratio', 'bifurcation_angle_local', 'bifurcation_angle_remote'])
gmi_feature_names = np.array(['moment1', 'moment2', 'moment3','moment4', 'moment5', 'moment6',  'moment7', 'moment8',  'moment9', 'moment10',    'moment11', 'moment12', 'moment13', 'avgR'])

selected_features=['max_euclidean_distance','num_stems','num_bifurcations','average_contraction','parent_daughter_ratio']



#==============================================================================


##################################################################################################

csv_file = data_DIR+'/glfeatures.csv'
    # TODO: detect nan values

glf_featureList = []  # each row is a feature vector
gmi_featureList = []
fn_list=[]
with open (FEATURE_FILE,'r') as  f:
    for fn_line in f: # ignore the SWCFILE=* line
        fn_list.append(fn_line[8:-1])
        line_globalFeature = (f.next()).strip()
        glf = map(float,line_globalFeature.split('\t'))
        glf_featureList.append(glf)

        line_GMI = (f.next()).strip()
        gmi = map(float,line_GMI.split('\t'))
        gmi_featureList.append(gmi)



gl_matrix=np.array(glf_featureList)
gmi_matrix= np.array(gmi_featureList)



df = pd.DataFrame(data = gl_matrix, columns = ['num_nodes', 'soma_surface', 'num_stems','num_bifurcations', 'num_branches', 'num_of_tips',  'overall_width', 'overall_height',  'overall_depth', 'average_diameter',    'total_length', 'total_surface', 'total_volume', 'max_euclidean_distance',       'max_path_distance', 'max_branch_order',  'average_contraction', 'average fragmentation', 'parent_daughter_ratio', 'bifurcation_angle_local', 'bifurcation_angle_remote'] )
df['swc_file'] = pd.Series(fn_list)


df.to_csv(csv_file)


