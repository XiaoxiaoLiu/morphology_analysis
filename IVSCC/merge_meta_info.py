__author__ = 'xiaoxiaoliu'

import platform
import pandas as pd
import numpy as np
import os
def copySnapshots(df_in, snapshots_dir, output_dir):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    swc_files = df_in['swc_file']
    if len(swc_files) > 0:
        for afile in swc_files:
            filename = snapshots_dir + '/' + afile.split('/')[-1] + '.BMP'
            if os.path.exists(filename):
                os.system("cp  " + filename + "  " + output_dir + "/\n")
    return



import glob
from PIL import Image
def assemble_screenshots(input_dir, output_image_file_name, size):
    files = glob.glob(input_dir + "/*.BMP")

    if len(files) < 1 :
        print "no BMP files found"
        return
    cols = len(files)
    rows = 1
    # if cols >100:
    #     assemble_image = Image.new("RGB", (size * cols,size*rows))
    #     cols = 10
    #     rows = np.ceil( len(files)/10.0)
    assemble_image = Image.new("RGB", (size * cols,size*rows))
    y=0
    for infile in files:
        if infile:
            im = Image.open(infile)
            im.thumbnail((size, size), Image.ANTIALIAS)
            assemble_image.paste(im, (y, 0))
            y += size
    print output_image_file_name
    assemble_image.save(output_image_file_name)

    return




def generateLinkerFileFromCSV(result_dir, csvfile, column_name=None, strip_path=True):
    df = pd.read_csv(csvfile)
    if (column_name == None):
        swc_files = df['swc_file']
        with open(result_dir + '/all.ano', 'w') as outf:
            for afile in swc_files:
                filename = afile
                if strip_path:
                    filename = afile.split('/')[-1]
                line = 'SWCFILE=' + filename + '\n'
                outf.write(line)
            outf.close()
        return

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
         'bifurcation_angle_remote','height_width_ratio','average_branch_length','length_surface_ratio'])



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
all_feature_file = data_DIR + '/preprocessed/features_with_db_tags_added.csv'





output_dir = data_DIR + '/clustering_results'
Meta_CSV_FILE = data_DIR + '/IVSCC_qual_calls_XiaoXiao_150cells_092915-UPDATED3.csv'
to_remove_filename = data_DIR +'/wrong_scale_removed.csv'



# require the following col names in the merged spread sheet
col_names = ['specimen_id','specimen_name','cre_line','layer','dendrite_type','swc_file','types']
all_feature_names = np.append(gl_feature_names, gmi_feature_names)
col_names.extend(all_feature_names)


# remove outliers identified in
df_remove  = pd.read_csv(to_remove_filename)
df_complete = pd.read_csv(all_feature_file)
df_complete_filter = df_complete[~df_complete['specimen_name'].isin(df_remove['specimen_name'])]
# !!!!! df_complete_filter drop  dendrite_type and layer

df_meta = pd.read_csv(Meta_CSV_FILE)


merged = pd.merge(df_complete_filter,df_meta,how='inner',on=['specimen_name'])

merged = merged[col_names]
merged[all_feature_names] = merged[all_feature_names].astype(float)


output_merged_csv = data_DIR+'/meta_merged_allFeatures.csv'
merged.to_csv(output_merged_csv,index=False)

#generateLinkerFileFromCSV(output_dir, output_merged_csv,'cre_line',False)

generateLinkerFileFromCSV(output_dir, output_merged_csv,None,False)


output_dir= data_DIR + '/figures/staci_types'
generateLinkerFileFromCSV(output_dir, output_merged_csv,'types',False)
swc_screenshot_folder =  data_DIR + "/figures/pw_aligned_bmps"
types = np.unique(merged['types'])
for type in types:
    print type
    df_type = merged[merged.types == type]
    type_format_string = '_'.join(type.split(" "))
    copySnapshots(df_type, swc_screenshot_folder, output_dir + '/' +type_format_string )
    assemble_screenshots(output_dir + '/' + type_format_string, output_dir + '/' +type_format_string+ '_assemble.png', 128)