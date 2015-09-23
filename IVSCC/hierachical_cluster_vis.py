import numpy as np
import matplotlib.pylab as pl
from scipy import stats

import pandas as pd
import seaborn as sns
import os

from scipy.cluster import hierarchy
import platform

from scipy.stats.stats import pearsonr

# program path on this machine
# ===================================================================
if (platform.system() == "Linux"):
   WORK_PATH = "/local1/xiaoxiaol/work"
else:
   WORK_PATH = "/Users/xiaoxiaoliu/work"

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
# ===================================================================


# zscore
ZSCORE_OUTLIER_THRESHOLD = 3.5
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




def zscores_matrix(featureArray, out_file, REMOVE_OUTLIER=1):
    normalized = zscore(featureArray)
    if REMOVE_OUTLIER:
         num_outliers = np.count_nonzero(normalized < -ZSCORE_OUTLIER_THRESHOLD) + np.count_nonzero(normalized > ZSCORE_OUTLIER_THRESHOLD)
         print(" found %d outliers" %num_outliers)
         if num_outliers >0 :
            normalized[normalized < -ZSCORE_OUTLIER_THRESHOLD] = -ZSCORE_OUTLIER_THRESHOLD
            normalized[normalized > ZSCORE_OUTLIER_THRESHOLD] = ZSCORE_OUTLIER_THRESHOLD
    df = pd.DataFrame(normalized)
    df.to_csv(out_file)
    return normalized


def zscore_features(merged, feature_names, out_file, REMOVE_OUTLIER=1):
    featureArray = merged[feature_names].astype(float)
    normalized = zscore(featureArray)
    if REMOVE_OUTLIER:
         num_outliers = np.count_nonzero(normalized < -ZSCORE_OUTLIER_THRESHOLD) + np.count_nonzero(normalized > ZSCORE_OUTLIER_THRESHOLD)
         print(" found %d outliers" %num_outliers)
         if num_outliers >0 :
            normalized[normalized < -ZSCORE_OUTLIER_THRESHOLD] = -ZSCORE_OUTLIER_THRESHOLD
            normalized[normalized > ZSCORE_OUTLIER_THRESHOLD] = ZSCORE_OUTLIER_THRESHOLD
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
         num_outliers = np.count_nonzero(normalized < -ZSCORE_OUTLIER_THRESHOLD) + np.count_nonzero(normalized > ZSCORE_OUTLIER_THRESHOLD)
         print(" found %d outliers" %num_outliers)
         if num_outliers >0 :
            normalized[normalized < -ZSCORE_OUTLIER_THRESHOLD] = -ZSCORE_OUTLIER_THRESHOLD
            normalized[normalized > ZSCORE_OUTLIER_THRESHOLD] = ZSCORE_OUTLIER_THRESHOLD

    for i in range(len(normalized)):
        queryFeature = normalized[i]  # each row is a feature vector
        scores = np.exp(-np.sum(abs(normalized-queryFeature)**2,1)/100)  #similarity
        #scores = np.sum(np.abs(normalized - queryFeature) ** 2, 1)  # distance
        distanceMatrix.append(scores)

    df = pd.DataFrame(distanceMatrix)
    df.to_csv(out_distanceMatrix_file)
    print("score sim matrix is saved to : " + out_distanceMatrix_file + "\n")
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

def heatmap_plot_distancematrix(distanceMatrix, merged, output_dir, title=None):
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
    if title:
        pl.title(title)
    # Legend for row and col colors
    print dendrite_types
    for label in dendrite_types:
        g.ax_row_dendrogram.bar(0, 0, color=dendrite_type_lut[label], label=label, linewidth=0)
        g.ax_row_dendrogram.legend(loc="center", ncol=1)

    for label in cre_lines:
        g.ax_col_dendrogram.bar(0, 0, color=cre_line_lut[label], label=label, linewidth=0)
        g.ax_col_dendrogram.legend(loc="center", ncol=3)

    pl.title('Similarities')

    filename = output_dir + '/similarity_heatmap.png'
    pl.savefig(filename, dpi=300)
    print("save similarity matrix heatmap figure to :" + filename)
    pl.close()
    return g


def heatmap_plot_zscore(zscore_features, merged, output_dir, title=None):
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
                       xticklabels=True, yticklabels =False)
    if title:
        pl.title(title)
    # TODO : adjust creline tag size
    # print type(g.data)
    #print g.data.columns
    #crelines = g.data['cre_line']
    #g.ax_heatmap.set_yticklabels(crelines, fontsize=3)

    assignment = hierarchy.fcluster(r_linkage, 2, criterion="maxclust")

    # Legend for row and col colors
    for label in dendrite_types:
        g.ax_row_dendrogram.bar(0, 0, color=dendrite_type_lut[label], label=label, linewidth=0)
        g.ax_row_dendrogram.legend(loc="center", ncol=1)

    #for label in cre_lines:
    #   g.ax_col_dendrogram.bar(0, 0, color=cre_line_lut[label], label=label, linewidth=0)
    #   g.ax_col_dendrogram.legend(loc="center", ncol=3)


    #pl.show()
    pl.title('zscore')
    filename = output_dir + '/zscore_feature_heatmap.png'
    pl.savefig(filename, dpi=300)
    print("save zscore matrix heatmap figure to :" + filename)
    pl.close()
    return g



def output_clusters(assign_ids, df_all,feature_names,  output_dir):
    clusters = np.unique(assign_ids)
    num_cluster = len(clusters)

    df_zscores = zscore_features(df_all, feature_names, None, REMOVE_OUTLIER=1)

    print("there are %d cluster" %num_cluster)
    df_cluster = pd.DataFrame()
    df_zscore_cluster = pd.DataFrame()
    for i in clusters:
        ids = np.nonzero(assign_ids == i)[0]  # starting from  0
        df_cluster = df_all.loc[ids]

        csv_file = output_dir + '/ward_cluster_'+str(i) +'.csv'
        df_cluster.to_csv(csv_file)

        df_zscore_cluster = df_zscores.loc[ids]
        csv_file2 = output_dir + '/ward_cluster_zscore_'+str(i) +'.csv'
        df_zscore_cluster.to_csv(csv_file2)

        ano_file = output_dir + '/ward_cluster_'+str(i) +'.ano'
        generateLinkerFileFromDF(df_cluster,ano_file, False)

        print("there are %d neurons in cluster %d" %(df_cluster.shape[0], i))
    return



#############################################################################################
def ward_cluster(df_all, feature_names, max_cluster_num, output_dir, fig_title=None):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    feature_array = df_all[feature_names].astype(float)

    simMatrix = distance_matrix(feature_array, output_dir + "/morph_features_similarity_matrix.csv", 1)

    # visualize heatmap using ward on similarity matrix
    out = heatmap_plot_distancematrix(simMatrix, df_all, output_dir, fig_title)
    linkage = out.dendrogram_row.calculated_linkage

    assignments = hierarchy.fcluster(linkage,max_cluster_num,criterion="maxclust")
    #hierarchy.dendrogram(linkage)

    ## put assignments into ano files and csv files
    output_clusters(assignments, df_all, feature_names,output_dir)


    ##### zscores  featuer plots
    df_zscores = zscore_features(merged, feature_names, output_dir + '/zscore.csv', 1)
    a = heatmap_plot_zscore(df_zscores, merged, output_dir, fig_title)


    return

#############################################################################################
def remove_correlated_features(df_all, feature_names, coef_threshold = 0.99):
    num_features = len(feature_names)
    removed_names=[]
    for i in range(num_features):
        if not feature_names[i] in removed_names:
            a = df_all[feature_names[i]].astype(float)

            for j in range(i+1, num_features):
                if not feature_names[j] in removed_names:
                    b = df_all[feature_names[j]].astype(float)
                    corrcoef = pearsonr(a,b)
                    if (corrcoef[0] > coef_threshold):
                        removed_names.append(feature_names[j])
                        print("highly correlated:[" +feature_names[i]+", "+feature_names[j]+" ]")

    subset_features_names = feature_names.tolist()
    for i in range(len(removed_names)):
        if removed_names[i] in subset_features_names:
              print ("remove "+removed_names[i])
              subset_features_names.remove(removed_names[i])

    return np.asarray(subset_features_names)

##################################################################################################
######################## WARD  Clustering  main function           ###############################
##################################################################################################

# merge all info, waiting to get cell_shape tags....
#df_type = pd.read_csv(data_DIR+'/../custom_report-IVSCC_classification-April_2015.csv')
#merged = pd.merge(df_complete,df_type,how='inner',on=['specimen_name'])
#merged.to_csv(data_DIR+'/merged_allFeatures.csv')

# To qualitative look though crelines
#generateLinkerFileFromCSV(data_DIR+'/original',data_DIR +'/merged_allFeatures.csv','cre_line')

all_feature_merged_file = data_DIR + '/preprocessed/features_with_db_tags.csv'
generateLinkerFileFromCSV(data_DIR + '/preprocessed', all_feature_merged_file, 'cre_line')
merged = pd.read_csv(all_feature_merged_file)

merged[all_feature_names]= merged[all_feature_names].astype(float)

cre_lines = np.unique(merged['cre_line'])
num_of_crelines = len(cre_lines)

output_dir = data_DIR+'/clustering_results'
if  not os.path.exists(output_dir):
      os.mkdir(output_dir)



PLOT_DISTANCE_MATRIX = 1
if PLOT_DISTANCE_MATRIX:
    ward_cluster(merged, all_feature_names, num_of_crelines, output_dir +'/ward_all_features')
    ward_cluster(merged, gmi_feature_names, num_of_crelines, output_dir +'/ward_GMI_features')


#################################################################################################
#############   feature selection

redundancy_removed_features_names = remove_correlated_features(merged, all_feature_names,0.98)
ward_cluster(merged, redundancy_removed_features_names, num_of_crelines,output_dir +'/ward_rr_all_features')
