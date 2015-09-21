import numpy as np
import matplotlib.pylab as pl
from scipy import stats

import pandas as pd
import seaborn as sns

from scipy.spatial import distance
from scipy.cluster import hierarchy



# program path on this machine
#===================================================================
WORK_PATH="/Users/xiaoxiaoliu/work"
MRMR= WORK_PATH+"/src/mrmr_c_src/mrmr"
V3D="qv3d"

########################################## data dir
data_DIR= WORK_PATH + "/data/lims2/0903_filtered_ephys_qc"
#########################################################
data_linker_file =  data_DIR+'/original/mylinker.ano'
preprocessed_data_linker_file = data_DIR+'/preprocessed/mylinker.ano'
FEATURE_FILE = data_DIR + '/preprocessed/prep_features.nfb'

gl_feature_names= np.array(['num_nodes', 'soma_surface', 'num_stems','num_bifurcations', 'num_branches', 'num_of_tips',  'overall_width', 'overall_height',  'overall_depth', 'average_diameter',    'total_length', 'total_surface', 'total_volume', 'max_euclidean_distance',       'max_path_distance', 'max_branch_order',  'average_contraction', 'average fragmentation', 'parent_daughter_ratio', 'bifurcation_angle_local', 'bifurcation_angle_remote'])

gmi_feature_names = np.array(['moment1', 'moment2', 'moment3','moment4', 'moment5', 'moment6',  'moment7', 'moment8',  'moment9', 'moment10',    'moment11', 'moment12', 'moment13'])

selected_features=['max_euclidean_distance','num_stems','num_bifurcations','average_contraction','parent_daughter_ratio']

all_feature_names = np.append(gl_feature_names,gmi_feature_names)
#===================================================================
def zscore(features):
    zscores = stats.zscore(features,0)
    return zscores


def normalizeFeatures(features):
    meanFeatures = np.mean(features,0)
    stdFeatures = np.std(features, 0)
    if np.count_nonzero(stdFeatures)< len(stdFeatures):
          print "zero detected"
          print stdFeatures
    normalized = (features - meanFeatures)/stdFeatures
    return normalized


def plotFeatureVector(featureArray,fig_title):
    normalized = normalizeFeatures(featureArray)
    pl.figure()
    pl.imshow(normalized,interpolation='none')
    pl.colorbar()
    pl.title(fig_title)
    pl.xlabel('feature ID')
    pl.ylabel('neuron ID')
    pl.show()
    pl.savefig(data_DIR+'/'+fig_title+'.png')


def zscores_matrix(featureArray, out_file, REMOVE_OUTLIER=1):
    zscoreMatrix=[]
    normalized= zscore(featureArray)
    #normalized = normalizeFeatures(featureArray)

    # remove outliers!!!
    normalized[normalized < -3]  =-3
    normalized[normalized > 3] = 3

    df = pd.DataFrame(normalized)
    df.to_csv(out_file)
    return normalized

def distance_matrix(featureArray,out_distanceMatrix_file, REMOVE_OUTLIER=1):
    distanceMatrix =[]
    normalized= zscore(featureArray)
    #normalized = normalizeFeatures(featureArray)

    # remove outliers!!!
    normalized[normalized < -3]  =-3
    normalized[normalized > 3] = 3

    for i in range(len(normalized)):
        queryFeature = normalized[i] # each row is a feature vector
        #scores = np.exp(-np.sum(abs(normalized-queryFeature)**2,1)/100)
        scores = np.sum(np.abs(normalized-queryFeature)**2,1)
        distanceMatrix.append(scores)

    df = pd.DataFrame(distanceMatrix)
    df.to_csv(out_distanceMatrix_file)
    return distanceMatrix



def generateLinkerFileFromCSV(result_dir, csvfile, column_name):
	df = pd.read_csv(csvfile)
	types = df[column_name]
	for atype in np.unique(types):
             idxs = np.nonzero(types==atype)[0]
             swc_files = df['swc_file']
             with open(result_dir+'/'+atype+'.ano','w') as outf:
        	   for afile in swc_files[idxs]:
                       filename = afile.split('/')[-1]
                       line='SWCFILE='+filename+'\n'
                       outf.write(line)
                   outf.close()


##############  heatmap plot: hierachical clustering  ########

def heatmap_plot_distancematrix(distanceMatrix,merged, output_dir, title):
    pl.figure()

    # Create a custom palette for creline colors
    cre_lines = np.unique(merged['cre_line'])
    cre_line_pal = sns.color_palette("hls",len(cre_lines))
    cre_line_lut = dict(zip(cre_lines, cre_line_pal)) # map creline type to color
    creline_colors = merged['cre_line'].map(cre_line_lut)

    # Create a custom palette for dendrite_type colors
    dendrite_types= np.unique(merged['dendrite_type'])
    dendrite_type_pal = sns.color_palette("hls",len(dendrite_types))
    dendrite_type_lut = dict(zip(dendrite_types, dendrite_type_pal)) 
    dendritetype_colors = merged['dendrite_type'].map(dendrite_type_lut)

    # Create a custom colormap for the heatmap values
    cmap = sns.diverging_palette(240,10, as_cmap=True)


    # PLOT
    g = sns.clustermap(distanceMatrix, method='ward',metric='euclidean',linewidths=0.0, row_colors = dendritetype_colors,col_colors=creline_colors,  cmap=cmap,xticklabels=False, yticklabels=False)

    # Legend for row and col colors
    for label in dendrite_types:
            g.ax_row_dendrogram.bar(0, 0, color=dendrite_type_lut[label], label=label, linewidth=0)
            g.ax_row_dendrogram.legend(loc="center", ncol=1)

    for label in cre_lines:
            g.ax_col_dendrogram.bar(0, 0, color=cre_line_lut[label], label=label, linewidth=0)
            g.ax_col_dendrogram.legend(loc="center", ncol=3)

    pl.title(title)
    #pl.show()
    pl.savefig(output_dir+ '/distance_matrix_heatmap'+title+'.pdf', dpi=600)
    pl.close()

def heatmap_plot_zscore(zscore_features,merged, output_dir, title):
    pl.figure()

    # Create a custom palette for creline colors
    cre_lines = np.unique(merged['cre_line'])
    cre_line_pal = sns.color_palette("hls",len(cre_lines))
    cre_line_lut = dict(zip(cre_lines, cre_line_pal)) # map creline type to color
    creline_colors = merged['cre_line'].map(cre_line_lut)

    # Create a custom palette for dendrite_type colors
    dendrite_types= np.unique(merged['dendrite_type'])
    dendrite_type_pal = sns.color_palette("hls",len(dendrite_types))
    dendrite_type_lut = dict(zip(dendrite_types, dendrite_type_pal)) 
    dendritetype_colors = merged['dendrite_type'].map(dendrite_type_lut)

    # Create a custom colormap for the heatmap values
    cmap = sns.diverging_palette(240,10, as_cmap=True)

    r_linkage = hierarchy.linkage( distance.pdist(zscore_features), method='ward', metric = 'euclidean')
    c_linkage = hierarchy.linkage( distance.pdist(zscore_features.T), method='ward', metric = 'euclidean')

    # PLOT
    g = sns.clustermap(zscore_features, row_linkage = r_linkage, col_linkage = c_linkage, method='ward',metric='euclidean',linewidths=0.0, row_colors = dendritetype_colors,col_colors=creline_colors,  cmap=cmap,xticklabels=False, yticklabels=False)

    # Legend for row and col colors
    for label in dendrite_types:
            g.ax_row_dendrogram.bar(0, 0, color=dendrite_type_lut[label], label=label, linewidth=0)
            g.ax_row_dendrogram.legend(loc="center", ncol=1)

    for label in cre_lines:
            g.ax_col_dendrogram.bar(0, 0, color=cre_line_lut[label], label=label, linewidth=0)
            g.ax_col_dendrogram.legend(loc="center", ncol=3)

    pl.title(title)
    #pl.show()
    pl.savefig(output_dir+ '/zscore_feature_heatmap'+title+'.pdf', dpi=600)
    pl.close()


##################################################################################################

# merge all info
#df_type = pd.read_csv(data_DIR+'/../custom_report-IVSCC_classification-April_2015.csv')
#merged = pd.merge(df_complete,df_type,how='inner',on=['specimen_name'])
#merged.to_csv(data_DIR+'/merged_allFeatures.csv')

# To qualitative look though crelines
#generateLinkerFileFromCSV(data_DIR+'/original',data_DIR +'/merged_allFeatures.csv','cre_line')

all_feature_merged_file = data_DIR +'/preprocessed/features_with_db_tags.csv'
generateLinkerFileFromCSV(data_DIR+'/preprocessed',all_feature_merged_file,'cre_line')
merged = pd.read_csv(all_feature_merged_file)
allFeatures = merged[all_feature_names].astype(float)

# score matrix saves zscore distance

# all features
distanceMatrix_csv_file = data_DIR+"/preprocessed/morph_features_distance_matrix.csv"
distanceMatrix = distance_matrix(allFeatures,distanceMatrix_csv_file,1)
print("score matrix is saved to : "+distanceMatrix_csv_file + "\n")

#heatmap_plot(distanceMatrix, merged,data_DIR, '34 Morph Features')
zscores = zscores_matrix(allFeatures, data_DIR+"/preprocessed/morph_features_zscores_matrix.csv", 1)
heatmap_plot_zscore(zscores, merged,data_DIR, '34 Morph Features')

# gmi features
gmiFeatures = merged[gmi_feature_names].astype(float)
distanceMatrix_csv_file = data_DIR+"/preprocessed/morph_gmi_distance_matrix.csv"
distanceMatrix = distance_matrix(gmiFeatures,distanceMatrix_csv_file,1)
print("score matrix is saved to : "+distanceMatrix_csv_file + "\n")

#heatmap_plot(distanceMatrix, merged,data_DIR ,'13 GMI Features')

# gl features
glFeatures = merged[gl_feature_names].astype(float)
distanceMatrix_csv_file = data_DIR+"/preprocessed/morph_gl_distance_matrix.csv"
distanceMatrix = distance_matrix(glFeatures,distanceMatrix_csv_file,1)
print("score matrix is saved to : "+distanceMatrix_csv_file + "\n")

#heatmap_plot(distanceMatrix, merged,data_DIR ,'21 L-M Features')
