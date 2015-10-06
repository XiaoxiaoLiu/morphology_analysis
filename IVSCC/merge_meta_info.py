__author__ = 'xiaoxiaoliu'

import platform
import pandas as pd
import numpy as np




def generateLinkerFileFromCSV(result_dir, csvfile, column_name, strip_path=True):
    df = pd.read_csv(csvfile)
    types = df[column_name]
    for atype in np.unique(types):
        idxs = np.nonzero(types == atype)[0]
        swc_files = df['swc_file']
        with open(result_dir + '/' + atype + '.ano', 'w') as outf:
            for afile in swc_files[idxs]:
                filename = afile
                if strip_path:
                    filename = afile.split('/')[-1]
                line = 'SWCFILE=' + filename + '\n'
                outf.write(line)
            outf.close()
    return


gl_feature_names = np.array(
        ['num_nodes', 'soma_surface', 'num_stems', 'num_bifurcations', 'num_branches', 'num_of_tips',
         'overall_width', 'overall_height', 'overall_depth', 'average_diameter', 'total_length',
         'total_surface', 'total_volume', 'max_euclidean_distance', 'max_path_distance', 'max_branch_order',
         'average_contraction', 'average fragmentation', 'parent_daughter_ratio', 'bifurcation_angle_local',
         'bifurcation_angle_remote'])

# remove scales
gl_feature_names_inv = np.array(
    ['num_nodes', 'soma_surface', 'num_stems', 'num_bifurcations', 'num_branches', 'num_of_tips',
     'average_diameter', 'total_length',
     'total_surface', 'total_volume', 'max_euclidean_distance', 'max_path_distance', 'max_branch_order',
     'average_contraction', 'average fragmentation', 'parent_daughter_ratio', 'bifurcation_angle_local',
     'bifurcation_angle_remote'])

gmi_feature_names = np.array(
    ['moment1', 'moment2', 'moment3', 'moment4', 'moment5', 'moment6', 'moment7', 'moment8',
     'moment9', 'moment10', 'moment11', 'moment12', 'moment13'])  ### removed ave_R




if (platform.system() == "Linux"):
    WORK_PATH = "/local1/xiaoxiaol/work"
else:
    WORK_PATH = "/Users/xiaoxiaoliu/work"

data_DIR = WORK_PATH + "/data/lims2/0923_pw_aligned"
all_feature_file = data_DIR + '/preprocessed/features_with_db_tags_or.csv'
output_dir = data_DIR + '/clustering_results'
Meta_CSV_FILE = data_DIR + '/IVSCC_qual_calls_XiaoXiao_150cells_092915.csv'





# require the following col names in the merged spread sheet
col_names = ['specimen_id','specimen_name','cre_line','layer_corrected','dendrite_type','swc_file','types']
all_feature_names = np.append(gl_feature_names, gmi_feature_names)
col_names.extend(all_feature_names)



df_complete = pd.read_csv(all_feature_file)
df_meta = pd.read_csv(Meta_CSV_FILE)

merged = pd.merge(df_complete,df_meta,how='inner',on=['specimen_name'])

merged = merged[col_names]
merged[all_feature_names] = merged[all_feature_names].astype(float)


output_merged_csv = data_DIR+'/meta_merged_allFeatures.csv'
merged.to_csv(output_merged_csv,index=False)

generateLinkerFileFromCSV(output_dir, output_merged_csv,'cre_line',False)

generateLinkerFileFromCSV(output_dir, output_merged_csv,'types',False)
