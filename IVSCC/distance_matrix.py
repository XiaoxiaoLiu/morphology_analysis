import numpy as np
import matplotlib.pylab as pl
from scipy import stats

import pandas as pd
import seaborn as sns

from scipy.spatial import distance
from scipy.cluster import hierarchy


# program path on this machine
# ===================================================================
WORK_PATH = "/local1/xiaoxiaol/work"

########################################## data dir
data_DIR = WORK_PATH + "/data/lims2/0903_filtered_ephys_qc"
#########################################################
data_linker_file = data_DIR + '/original/mylinker.ano'
preprocessed_data_linker_file = data_DIR + '/preprocessed/mylinker.ano'
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
# ===================================================================
def zscore(features):
    zscores = stats.zscore(features, 0)
    return zscores


def normalizeFeatures(features):
    meanFeatures = np.mean(features, 0)
    stdFeatures = np.std(features, 0)
    if np.count_nonzero(stdFeatures) < len(stdFeatures):
        print "zero detected"
        print stdFeatures
    normalized = (features - meanFeatures) / stdFeatures
    return normalized


def plotFeatureVector(featureArray, fig_title):
    normalized = normalizeFeatures(featureArray)
    pl.figure()
    pl.imshow(normalized, interpolation='none')
    pl.colorbar()
    pl.title(fig_title)
    pl.xlabel('feature ID')
    pl.ylabel('neuron ID')
    pl.show()
    pl.savefig(data_DIR + '/' + fig_title + '.png')
    return


OUTLIER_THRESHOLD = 5

def zscores_matrix(featureArray, out_file, REMOVE_OUTLIER=1):
    normalized = zscore(featureArray)
    if REMOVE_OUTLIER:
        normalized[normalized < -OUTLIER_THRESHOLD] = -OUTLIER_THRESHOLD
        normalized[normalized > OUTLIER_THRESHOLD] = OUTLIER_THRESHOLD
    df = pd.DataFrame(normalized)
    df.to_csv(out_file)
    return normalized


def zscore_features(merged, feature_names, out_file, REMOVE_OUTLIER=1):
    featureArray = merged[feature_names].astype(float)
    normalized = zscore(featureArray)
    if REMOVE_OUTLIER:
        normalized[normalized < -OUTLIER_THRESHOLD] = -OUTLIER_THRESHOLD
        normalized[normalized > OUTLIER_THRESHOLD] = OUTLIER_THRESHOLD
    df = pd.DataFrame(normalized)
    df.columns = feature_names
    if out_file:
        df.to_csv(out_file)
        print("save to " + out_file )
    return df


def distance_matrix(featureArray, out_distanceMatrix_file, REMOVE_OUTLIER=1):
    distanceMatrix = []
    normalized = zscore(featureArray)
    #normalized = normalizeFeatures(featureArray)

    # remove outliers!!!
    if REMOVE_OUTLIER:
        normalized[normalized < -OUTLIER_THRESHOLD] = -OUTLIER_THRESHOLD
        normalized[normalized > OUTLIER_THRESHOLD] = OUTLIER_THRESHOLD

    for i in range(len(normalized)):
        queryFeature = normalized[i]  # each row is a feature vector
        scores = np.exp(-np.sum(abs(normalized-queryFeature)**2,1)/100)  #similarity
        #scores = np.sum(np.abs(normalized - queryFeature) ** 2, 1)  # distance
        distanceMatrix.append(scores)

    df = pd.DataFrame(distanceMatrix)
    df.to_csv(out_distanceMatrix_file)
    return distanceMatrix


def generateLinkerFileFromDF(df_in, output_ano_file, strip_path= False):
    swc_files = df_in['swc_file']
    if len(swc_files) > 0:
        with open(output_ano_file, 'w') as outf:
            for afile in swc_files:
                filename = afile
                if strip_path:
                   filename = afile.split('/')[-1]
                line = 'SWCFILE=' + filename + '\n'
                outf.write(line)
            outf.close()
    return


def generateLinkerFileFromCSV(result_dir, csvfile, column_name, strip_path = True):
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


##############  heatmap plot: hierachical clustering  ########

def heatmap_plot_distancematrix(distanceMatrix, merged, output_dir, title):
    pl.figure()

    # Create a custom palette for creline colors
    cre_lines = np.unique(merged['cre_line'])
    cre_line_pal = sns.color_palette("hls", len(cre_lines))
    cre_line_lut = dict(zip(cre_lines, cre_line_pal))  # map creline type to color
    creline_colors = merged['cre_line'].map(cre_line_lut)

    # Create a custom palette for dendrite_type colors
    dendrite_types = np.unique(merged['dendrite_type'])
    dendrite_type_pal = sns.color_palette("hls", len(dendrite_types))
    dendrite_type_lut = dict(zip(dendrite_types, dendrite_type_pal))
    dendritetype_colors = merged['dendrite_type'].map(dendrite_type_lut)

    # Create a custom colormap for the heatmap values
    cmap = sns.diverging_palette(240, 10, as_cmap=True)


    g = sns.clustermap(distanceMatrix, method='ward', metric='euclidean', linewidths=0.0,
                       row_colors=dendritetype_colors, col_colors=creline_colors, cmap=cmap, xticklabels=False,
                       yticklabels=False)


    # Legend for row and col colors
    print dendrite_types
    for label in dendrite_types:
        g.ax_row_dendrogram.bar(0, 0, color=dendrite_type_lut[label], label=label, linewidth=0)
        g.ax_row_dendrogram.legend(loc="center", ncol=1)

    for label in cre_lines:
        g.ax_col_dendrogram.bar(0, 0, color=cre_line_lut[label], label=label, linewidth=0)
        g.ax_col_dendrogram.legend(loc="center", ncol=3)

    pl.title(title)
    #pl.show()

    filename = output_dir + '/similarity_heatmap' + title + '.pdf'
    pl.savefig(filename, dpi=600)
    print("save similarity matrix heatmap figure to :" + filename)
    pl.close()
    return g


def heatmap_plot_zscore(zscore_features, merged, output_dir, title):
    pl.figure()

    # Create a custom palette for creline colors
    cre_lines = np.unique(merged['cre_line'])
    cre_line_pal = sns.color_palette("hls", len(cre_lines))
    cre_line_lut = dict(zip(cre_lines, cre_line_pal))  # map creline type to color
    creline_colors = merged['cre_line'].map(cre_line_lut)

    # Create a custom palette for dendrite_type colors
    dendrite_types = np.unique(merged['dendrite_type'])
    dendrite_type_pal = sns.color_palette("hls", len(dendrite_types))
    dendrite_type_lut = dict(zip(dendrite_types, dendrite_type_pal))
    dendritetype_colors = merged['dendrite_type'].map(dendrite_type_lut)

    # Create a custom colormap for the heatmap values
    cmap = sns.diverging_palette(240, 10, as_cmap=True)

    r_linkage = hierarchy.linkage(zscore_features, method='ward', metric='euclidean')
    c_linkage = hierarchy.linkage(zscore_features.T, method='ward', metric='euclidean')

    # PLOT
    g = sns.clustermap(zscore_features, row_linkage=r_linkage, method='ward', metric='euclidean',
                       linewidths=0.0, row_colors=dendritetype_colors, cmap=cmap,
                       xticklabels=True, yticklabels = merged['cre_line'])
    g.ax_heatmap.set_yticklabels( merged['cre_line'],fontsize='small')


    assignment = hierarchy.fcluster(r_linkage, 2, criterion="maxclust")

    # Legend for row and col colors
    for label in dendrite_types:
        g.ax_row_dendrogram.bar(0, 0, color=dendrite_type_lut[label], label=label, linewidth=0)
        g.ax_row_dendrogram.legend(loc="center", ncol=1)

    #for label in cre_lines:
    #   g.ax_col_dendrogram.bar(0, 0, color=cre_line_lut[label], label=label, linewidth=0)
    #   g.ax_col_dendrogram.legend(loc="center", ncol=3)

    pl.title(title)
    pl.show()
    filename = output_dir + '/zscore_feature_heatmap' + title + '.pdf'
    pl.savefig(filename, dpi=600)
    print("save zscore matrix heatmap figure to :" + filename)
    pl.close()
    return c_linkage



def output_clusters(assign_ids, df_all,feature_names,  output_dir):
    clusters = np.unique(assign_ids)
    num_cluster = len(clusters)

    df_zscores =  zscore_features(df_all, feature_names, None, REMOVE_OUTLIER=1)

    print("there are %d cluster" %num_cluster)
    for i in clusters:
        ids = np.nonzero(assign_ids == i)[0]  # starting from  0
        df_cluster = df_all.loc[ids]

        csv_file = output_dir + '/ward_cluster_'+str(i) +'.csv'
        df_cluster.to_csv(csv_file)

        df_zscore_cluster = df_zscores.loc[ids]
        csv_file2 = output_dir + '/ward_cluster_zscore_'+str(i) +'.csv'
        df_zscore_cluster.to_csv(csv_file2)

        ano_file = output_dir + '/ward_cluster_'+str(i) +'.ano'
        generateLinkerFileFromDF(df_cluster,ano_file, True)
    return
##################################################################################################

# merge all info
#df_type = pd.read_csv(data_DIR+'/../custom_report-IVSCC_classification-April_2015.csv')
#merged = pd.merge(df_complete,df_type,how='inner',on=['specimen_name'])
#merged.to_csv(data_DIR+'/merged_allFeatures.csv')

# To qualitative look though crelines
#generateLinkerFileFromCSV(data_DIR+'/original',data_DIR +'/merged_allFeatures.csv','cre_line')

all_feature_merged_file = data_DIR + '/preprocessed/features_with_db_tags.csv'
generateLinkerFileFromCSV(data_DIR + '/preprocessed', all_feature_merged_file, 'cre_line')
merged = pd.read_csv(all_feature_merged_file)
allFeatures = merged[all_feature_names].astype(float)


cre_lines = np.unique(merged['cre_line'])
num_of_crelines = len(cre_lines)

PLOT_DISTANCE_MATRIX = 1
if PLOT_DISTANCE_MATRIX:
    # all features
    distanceMatrix_csv_file = data_DIR + "/preprocessed/morph_features_similarity_matrix.csv"
    distanceMatrix = distance_matrix(allFeatures, distanceMatrix_csv_file, 1)
    print("score matrix is saved to : " + distanceMatrix_csv_file + "\n")
    out = heatmap_plot_distancematrix(distanceMatrix, merged, data_DIR, '34 Morph Features')
    linkage = out.dendrogram_row.calculated_linkage

    # choose the number of clusteres by the crelines
    assignments = hierarchy.fcluster(linkage,num_of_crelines,criterion="maxclust")
    #hierarchy.dendrogram(linkage)

    ## put assignments into ano files and csv files
    output_clusters(assignments,merged, all_feature_names,data_DIR+'/preprocessed')


if 0:

    gmiFeatures = merged[gmi_feature_names].astype(float)
    distanceMatrix_csv_file = data_DIR + "/preprocessed/morph_gmi_similarity_matrix.csv"
    distanceMatrix = distance_matrix(gmiFeatures, distanceMatrix_csv_file, 1)
    print("score matrix is saved to : " + distanceMatrix_csv_file + "\n")
    heatmap_plot_distancematrix(distanceMatrix, merged, data_DIR, '13 GMI Features')

    # gl features
    glFeatures = merged[gl_feature_names].astype(float)
    distanceMatrix_csv_file = data_DIR + "/preprocessed/morph_gl_similarity_matrix.csv"
    distanceMatrix = distance_matrix(glFeatures, distanceMatrix_csv_file, 1)
    print("score matrix is saved to : " + distanceMatrix_csv_file + "\n")
    heatmap_plot_distancematrix(distanceMatrix, merged, data_DIR, '21 L-M Features')



PLOT_ZSCORE_HEATMAP = 0
if PLOT_ZSCORE_HEATMAP:
    # score matrix saves zscore distance
    prefix = data_DIR + "/preprocessed/morph_features_zscores"

    fn = prefix + '_all.csv'
    df_zscores = zscore_features(merged, all_feature_names, fn, 1)
    a = heatmap_plot_zscore(df_zscores, merged, data_DIR, 'all_Features')

if 0:
    fn = prefix + '_gmi.csv'
    df_zscores = zscore_features(merged, gmi_feature_names, fn, 1)
    a = heatmap_plot_zscore(df_zscores, merged, data_DIR, 'gmi_Features')

    fn = prefix + '_gl.csv'
    df_zscores = zscore_features(merged, gl_feature_names, fn, 1)
    a = heatmap_plot_zscore(df_zscores, merged, data_DIR, 'gl_Features')
