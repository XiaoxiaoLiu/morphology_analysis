__author__ = 'xiaoxiaol'


import numpy as np
import pylab as pl
import scipy
import pandas as pd
import seaborn as sns
import os
import sys, getopt
from scipy.cluster import hierarchy
import platform
from scipy.stats.stats import pearsonr
import scipy.stats as stats
from PIL import Image
import glob
from sklearn.metrics import silhouette_samples, silhouette_score
import math

from sklearn.cluster import AffinityPropagation
from sklearn import metrics
from itertools import cycle


####################################
ZSCORE_OUTLIER_THRESHOLD = 5
####################################

sns.set_context("poster")


def zscore(features, remove_outlier=0):
    zscores = scipy.stats.zscore(features, 0)
    # zscores = normalizeFeatures(features)
    return zscores


# def normalizeFeatures(features):
# meanFeatures = np.median(features, 0)
# stdFeatures = np.std(features, 0)
#     if np.count_nonzero(stdFeatures) < len(stdFeatures):
#         print "zero detected"
#         print stdFeatures
#     normalized = (features - meanFeatures) / stdFeatures
#     return normalized



#### need to be updated
# def distance_matrix(df_all, feature_names, out_distanceMatrix_file, REMOVE_OUTLIER=0):
#     feature_array = df_all[feature_names].astype(float)
#     distanceMatrix = []
#     normalized = zscore(feature_array)
#     #normalized = normalizeFeatures(feature_array)
#
#     if num_outliers > 0:
#         if not REMOVE_OUTLIER:  # only clp
#             normalized[normalized < -ZSCORE_OUTLIER_THRESHOLD] = -ZSCORE_OUTLIER_THRESHOLD
#             normalized[normalized > ZSCORE_OUTLIER_THRESHOLD] = ZSCORE_OUTLIER_THRESHOLD
#
#     for i in range(len(normalized)):
#         queryFeature = normalized[i]  # each row is a feature vector
#         scores = np.exp(-np.sum(abs(normalized - queryFeature) ** 2, 1) / 100)  #similarity
#         #scores = np.sum(np.abs(normalized - queryFeature) ** 2, 1)  # distance
#         distanceMatrix.append(scores)
#
#     df_dist = pd.DataFrame(distanceMatrix)
#     df_dist.to_csv(out_distanceMatrix_file, index=False)
#     print("score sim matrix is saved to : " + out_distanceMatrix_file + "\n")
#     return df_dist


def copySnapshots(df_in, snapshots_dir, output_dir):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    swc_files = df_in['swc_file_name']
    if len(swc_files) > 0:
        for afile in swc_files:
            filename = snapshots_dir + '/' + afile.split('/')[-1] + '.BMP'
            if os.path.exists(filename):
                os.system("cp  " + filename + "  " + output_dir + "/\n")
    return




def assemble_screenshots(input_dir, output_image_file_name, size):
    files = glob.glob(input_dir + "/*.BMP")

    assemble_image = Image.new("RGB", (size * len(files),size))

    y = 0
    for infile in files:
        im = Image.open(infile)
        im.thumbnail((size, size), Image.ANTIALIAS)
        assemble_image.paste(im, (y, 0))
        y += size

    assemble_image.save(output_image_file_name)

    return


def generateLinkerFileFromDF(df_in, output_ano_file, strip_path=False, swc_path=None):
    swc_files = df_in['swc_file_name']
    if len(swc_files) > 0:
        with open(output_ano_file, 'w') as outf:
            for afile in swc_files:

                if swc_path is not None:
                    filename = swc_path + '/'+afile
                else:
                    filename = afile

                if strip_path:
                    filename = afile.split('/')[-1]
                line = 'SWCFILE=' + filename + '\n'
                outf.write(line)
            outf.close()
    return




##############  heatmap plot: hierachical clustering  ########
#
# def heatmap_plot_distancematrix(df_distanceMatrix, merged, output_dir, title=None):
#     pl.figure()
#
#     # Create a custom palette for creline colors
#     cre_lines = np.unique(merged['cre_line'])
#     cre_line_pal = sns.color_palette("hls", len(cre_lines))
#     cre_line_lut = dict(zip(cre_lines, cre_line_pal))  # map creline type to color
#     creline_colors = merged['cre_line'].map(cre_line_lut)
#
#     # Create a custom palette for dendrite_type colors  thre colors
#     dendrite_types = np.unique(merged['dendrite_type'])
#     dendrite_type_pal = sns.color_palette(['white','gray','black'])
#     dendrite_type_lut = dict(zip(dendrite_types, dendrite_type_pal))
#     dendritetype_colors = merged['dendrite_type'].map(dendrite_type_lut)
#
#     # Create a custom colormap for the heatmap values
#     #cmap = sns.diverging_palette(240, 10, as_cmap=True)
#
#     g = sns.clustermap(df_distanceMatrix, method='ward', metric='euclidean', linewidths=0.0,
#                        row_colors=dendritetype_colors, col_colors=creline_colors, cmap=cmap, xticklabels=False,
#                        yticklabels=False)
#     if title:
#         pl.title(title)
#     # Legend for row and col colors
#
#     for label in dendrite_types:
#         pl.bar(0, 0, color=dendrite_type_lut[label], label=label, linewidth=0)
#         pl.legend(loc="center", ncol=1)
#
#     for label in cre_lines:
#         g.ax_col_dendrogram.bar(0, 0, color=cre_line_lut[label], label=label, linewidth=0)
#         g.ax_col_dendrogram.legend(loc="center", ncol=3)
#
#     pl.title('Similarities')
#
#     filename = output_dir + '/similarity_heatmap.png'
#     pl.savefig(filename, dpi=300)
#     print("save similarity matrix heatmap figure to :" + filename)
#     pl.close()
    return g


def plot_confusion_matrix(cm, xlabel, ylabel, xnames, ynames,  title='Confusion matrix', cmap=pl.cm.Blues):
    pl.grid(False)
    pl.imshow(cm, interpolation = 'none',cmap=cmap)
    pl.title(title)
    pl.colorbar()
    tick_marksx = np.arange(len(xnames))
    tick_marksy = np.arange(len(ynames))
    pl.xticks(tick_marksx, xnames)
    pl.yticks(tick_marksy, ynames)
    pl.tight_layout()
    pl.ylabel(ylabel)
    pl.xlabel(xlabel)


def heatmap_plot_zscore_ivscc(df_zscore_features, df_all, output_dir, title=None):

    # Create a custom palette for dendrite_type colors
    dendrite_types = [np.nan, 'aspiny', 'sparsely spiny', 'spiny']
    # dendrite_type_pal = sns.color_palette("coolwarm", len(dendrite_types))
    dendrite_type_pal = sns.color_palette(["gray","black","purple","red"])
    dendrite_type_lut = dict(zip(dendrite_types, dendrite_type_pal))
    dendrite_type_colors = df_all['dendrite_type'].map(dendrite_type_lut)


    # Create a custom palette for creline colors
    cre_lines = np.unique(df_all['cre_line'])
    print cre_lines
    cre_lines = ['Pvalb-IRES-Cre','Sst-IRES-Cre','Gad2-IRES-Cre', 'Htr3a-Cre_NO152',
                 'Nr5a1-Cre', 'Ntsr1-Cre','Rbp4-Cre_KL100' ,'Rorb-IRES2-Cre-D', 'Scnn1a-Tg2-Cre',
                 'Scnn1a-Tg3-Cre','Slc17a6-IRES-Cre','Cux2-CreERT2']

    cre_line_pal = sns.color_palette("BrBG", len(cre_lines))

    cre_line_lut = dict(zip(cre_lines, cre_line_pal))  # map creline type to color
    cre_line_colors = df_all['cre_line'].map(cre_line_lut)


    # layers = np.unique(df_all['layer'])
    # layer_pal = sns.light_palette("green", len(layers))
    # layer_lut = dict(zip(layers, layer_pal))
    # layer_colors = df_all['layer'].map(layer_lut)

    # # only if types are available
    # types = np.unique(df_all['types'])
    # #reorder
    # types = ['NGC','multipolar','symm', 'bitufted','bipolar','tripod', 'Martinotti','cortico-cortical', 'cortico-thal','non-tufted', 'short-thick-tufted', 'tufted','thick-tufted']
    # type_pal = sns.color_palette("coolwarm", len(types))#  sns.diverging_palette(220, 20, n=len(types))# sns.color_palette("husl", len(types))
    # type_lut = dict(zip(types, type_pal))
    # type_colors = df_all['types'].map(type_lut)


    # Create a custom colormap for the heatmap values
    #cmap = sns.diverging_palette(240, 10, as_cmap=True)

    linkage = hierarchy.linkage(df_zscore_features, method='ward', metric='euclidean')

    data = df_zscore_features.transpose()
    row_linkage = hierarchy.linkage(data, method='ward', metric='euclidean')
    feature_order = hierarchy.leaves_list(row_linkage)

    #print data.index
    matchIndex = [data.index[x] for x in feature_order]
    #print matchIndex
    data = data.reindex(matchIndex)


    g = sns.clustermap(data, row_cluster = False, col_linkage=linkage, method='ward', metric='euclidean',
                       linewidths = 0.0,col_colors = [cre_line_colors,dendrite_type_colors],
                       cmap = sns.cubehelix_palette(light=1, as_cmap=True),figsize=(40,20))
    #g.ax_heatmap.xaxis.set_xticklabels()
    pl.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=90 )
    pl.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    pl.subplots_adjust(left=0.1, bottom=0.5, right=0.9, top=0.95)  # !!!!!

    #pl.tight_layout( fig, h_pad=20.0, w_pad=20.0)


    if title:
        pl.title(title)
    location ="best"
    num_cols=1
    # Legend for row and col colors

    for label in cre_lines:
         g.ax_row_dendrogram.bar(0, 0, color=cre_line_lut[label], label=label, linewidth=0.0)
         g.ax_row_dendrogram.legend(loc=location, ncol=num_cols,borderpad=0)

    for i in range(3):
        g.ax_row_dendrogram.bar(0, 0, color = "white", label=" ", linewidth=0)
        g.ax_row_dendrogram.legend(loc=location, ncol=num_cols, borderpad=0.0)

    # for label in layers:
    #      pl.bar(0, 0, color=layer_lut[label], label=label, linewidth=1)
    #      pl.legend(loc="left", ncol=2,borderpad=0.5)

    #
    # for label in types:
    #      g.ax_row_dendrogram.bar(0, 0, color=type_lut[label], label=label,linewidth=0)
    #      g.ax_row_dendrogram.legend(loc=location, ncol=num_cols,borderpad=0.0)
    #
    #
    # g.ax_row_dendrogram.bar(0, 0, color = "white", label=" ", linewidth=0)
    # g.ax_row_dendrogram.legend(loc=location, ncol=num_cols, borderpad=0.0)


    for label in dendrite_types:
        g.ax_row_dendrogram.bar(0, 0, color = dendrite_type_lut[label], label=label, linewidth=0)
        g.ax_row_dendrogram.legend(loc=location, ncol= num_cols, borderpad=0.0)


    filename = output_dir + '/zscore_feature_heatmap.png'
    pl.savefig(filename, dpi=300)
    #pl.show()
    print("save zscore matrix heatmap figure to :" + filename)
    pl.close()
    return linkage


def heatmap_plot_zscore_bbp(df_zscore_features, df_all, output_dir, title=None):

    print "heatmap plot"
    metric ='m-type'
    mtypes = np.unique(df_all[metric])
    print mtypes
    mtypes_pal = sns.color_palette("hls", len(mtypes))

    mtypes_lut = dict(zip(mtypes, mtypes_pal))  # map creline type to color
    mtypes_colors = df_all[metric].map(mtypes_lut)


    layers = np.unique(df_all['layer'])
    layer_pal = sns.light_palette("green", len(layers))
    layers_lut = dict(zip(layers, layer_pal))
    layer_colors = df_all['layer'].map(layers_lut)


    # Create a custom colormap for the heatmap values
    #cmap = sns.diverging_palette(240, 10, as_cmap=True)

    linkage = hierarchy.linkage(df_zscore_features, method='ward', metric='euclidean')

    data = df_zscore_features.transpose()
    row_linkage = hierarchy.linkage(data, method='ward', metric='euclidean')
    feature_order = hierarchy.leaves_list(row_linkage)

    #print data.index
    matchIndex = [data.index[x] for x in feature_order]
    #print matchIndex
    data = data.reindex(matchIndex)


    g = sns.clustermap(data, row_cluster = False, col_linkage=linkage, method='ward', metric='euclidean',
                       linewidths = 0.0,col_colors = [mtypes_colors,layer_colors],
                       cmap = sns.cubehelix_palette(light=1, as_cmap=True),figsize=(40,20))
    #g.ax_heatmap.xaxis.set_xticklabels()
    pl.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=90 )
    pl.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    pl.subplots_adjust(left=0.1, bottom=0.5, right=0.9, top=0.95)  # !!!!!

    #pl.tight_layout( fig, h_pad=20.0, w_pad=20.0)


    if title:
        pl.title(title)
    location ="best"
    num_cols=1
    # Legend for row and col colors

    for label in mtypes:
         g.ax_row_dendrogram.bar(0, 0, color=mtypes_lut[label], label=label, linewidth=0.0)
         g.ax_row_dendrogram.legend(loc=location, ncol=num_cols,borderpad=0)

    for i in range(3):
        g.ax_row_dendrogram.bar(0, 0, color = "white", label=" ", linewidth=0)
        g.ax_row_dendrogram.legend(loc=location, ncol=num_cols, borderpad=0.0)

    for label in layers:
         g.ax_row_dendrogram.bar(0, 0, color=layers_lut[label], label=label, linewidth=0.0)
         g.ax_row_dendrogram.legend(loc=location, ncol=num_cols,borderpad=0)

    filename = output_dir + '/zscore_feature_heatmap.png'
    pl.savefig(filename, dpi=300)
    #pl.show()
    print("save zscore matrix heatmap figure to :" + filename)
    pl.close()
    return linkage


##########################   feature selection   ########################
def remove_correlated_features(df_all, feature_names, coef_threshold=0.98):
    num_features = len(feature_names)
    removed_names = []
    for i in range(num_features):
        if not feature_names[i] in removed_names:
            a = df_all[feature_names[i]].astype(float)

            for j in range(i + 1, num_features):
                if not feature_names[j] in removed_names:
                    b = df_all[feature_names[j]].astype(float)
                    corrcoef = pearsonr(a, b)
                    if (corrcoef[0] > coef_threshold):
                        removed_names.append(feature_names[j])
                        print("highly correlated:[" + feature_names[i] + ", " + feature_names[j] + " ]")

    subset_features_names = feature_names.tolist()
    for i in range(len(removed_names)):
        if removed_names[i] in subset_features_names:
            print ("remove " + removed_names[i])
            subset_features_names.remove(removed_names[i])

    return np.asarray(subset_features_names)


#######################################  cluster evaluations ##################
def delta(ck, cl):
    values = np.ones([len(ck), len(cl)]) * 10000

    for i in range(0, len(ck)):
        for j in range(0, len(cl)):
            values[i, j] = np.linalg.norm(ck[i] - cl[j])

    return np.min(values)


def big_delta(ci):
    values = np.zeros([len(ci), len(ci)])

    for i in range(0, len(ci)):
        for j in range(0, len(ci)):
            values[i, j] = np.linalg.norm(ci[i] - ci[j])

    return np.max(values)


def dunn(k_list):
    """ Dunn index [CVI]

    Parameters
    ----------
    k_list : list of np.arrays
        A list containing a numpy array for each cluster |c| = number of clusters
        c[K] is np.array([N, p]) (N : number of samples in cluster K, p : sample dimension)
    """
    deltas = np.ones([len(k_list), len(k_list)]) * 1000000
    big_deltas = np.zeros([len(k_list), 1])
    l_range = range(0, len(k_list))

    for k in l_range:
        for l in (l_range[0:k] + l_range[k + 1:]):
            deltas[k, l] = delta(k_list[k], k_list[l])

        big_deltas[k] = big_delta(k_list[k])

    di = np.min(deltas) / np.max(big_deltas)
    return di


###############################  cluster specific features #####

def cluster_specific_features(df_all, assign_ids, feature_names, output_csv_fn):
    #student t to get cluster specific features

    labels=[]
    clusters = np.unique(assign_ids)
    num_cluster = len(clusters)
    df_pvalues =  pd.DataFrame(index = feature_names, columns = clusters)

    for cluster_id in clusters:

        ids_a = np.nonzero(assign_ids == cluster_id)[0]  # starting from  0
        ids_b = np.nonzero(assign_ids != cluster_id)[0]  # starting from  0
        labels.append("C"+str(cluster_id) + "("+ str(len(ids_a))+")" )
        for feature in feature_names:
            a = df_all.iloc[ids_a][feature]
            b = df_all.iloc[ids_b][feature]
            t_stats,pval = stats.ttest_ind(a,b,equal_var=False)
            df_pvalues.loc[feature,cluster_id] = -np.log10(pval)


    df_pvalues.to_csv(output_csv_fn)



    ### visulaize
    df_pvalues.index.name = "Features"
    df_pvalues.columns.name ="Clusters"
    d=df_pvalues[df_pvalues.columns].astype(float)
    g = sns.heatmap(data=d,linewidths=0.1)
     #               cmap =sns.color_palette("coolwarm",7, as_cmap=True))

    g.set_xticklabels(labels)
    pl.yticks(rotation=0)
    pl.xticks(rotation=90)
    pl.subplots_adjust(left=0.5, right=0.9, top=0.9, bottom=0.1)
    pl.title('-log10(P value)')
    filename = output_csv_fn + '.png'
    pl.savefig(filename, dpi=300)
    #pl.show()
    pl.close()


    return df_pvalues


#############################################################################################
def get_zscore_features(df_all, feature_names, out_file, REMOVE_OUTLIER=0,
                        zscore_threshold=ZSCORE_OUTLIER_THRESHOLD):  # if remove_outlier ==0 , just clip at threshold
    featureArray = df_all[feature_names].astype(float)
    featureArray.fillna(0,inplace=True)  ### might introduce some bias

    normalized = zscore(featureArray)
    # normalized = featureArray
    # normalized[~np.isnan(featureArray)] = zscore(featureArray[~np.isnan(featureArray)])

    num_outliers = np.count_nonzero(normalized < -zscore_threshold) + np.count_nonzero(
        normalized > zscore_threshold)
    print("Found %d  |z score| > %f in zscore matrix :" % (num_outliers, zscore_threshold) )

    df_all_modified = df_all
    df_outliers = pd.DataFrame()
    if num_outliers > 0:
        if not REMOVE_OUTLIER:  # just clip
            normalized[normalized < -zscore_threshold] = -zscore_threshold
            normalized[normalized > zscore_threshold] = zscore_threshold
        # else:
        #     outliers_l = np.nonzero(normalized < -zscore_threshold)
        #     outliers_h = np.nonzero(normalized > zscore_threshold)
        #     outlier_index = np.unique((np.append(outliers_l[0], outliers_h[0])))
        #
        #     # remove outlier rows
        #     df_all_modified = df_all_modified.drop(df_all_modified.index[outlier_index])
        #     normalized = np.delete(normalized, outlier_index, 0)
        #
        #     # re-zscoring and clipping
        #     # m_featureArray = df_all_modified[feature_names].astype(float)
        #     # normalized = zscore(m_featureArray)
        #     # normalized[normalized < -zscore_threshold] = -zscore_threshold
        #     # normalized[normalized > zscore_threshold] = zscore_threshold
        #
        #
        #     print("Removed %d outlier neurons" % len(outlier_index))
        #
        #     df_outliers = df_all.iloc[outlier_index]

    df_z = pd.DataFrame(normalized)
    df_z.columns = feature_names
    df_z.index = df_all['swc_file_name']

    if out_file:
        df_z.to_csv(out_file, index=True)
        print("save to " + out_file )

    if (df_z.shape[0] != df_all_modified.shape[0]):
        print ("error:  the sample size of the zscore and the original table does not match!")

    return df_z, df_all_modified, df_outliers


#############################################################################################
def output_single_cluster_results(df_cluster, output_dir, output_prefix, snapshots_dir=None, swc_path = None):
    csv_file = output_dir + '/' + output_prefix + '.csv'
    df_cluster.to_csv(csv_file, index=False)

    ano_file = output_dir + '/' + output_prefix + '.ano'
    generateLinkerFileFromDF(df_cluster, ano_file, False, swc_path)
    # copy bmp vaa3d snapshots images over

    if (snapshots_dir):
        copySnapshots(df_cluster, snapshots_dir, output_dir + '/' + output_prefix)
        assemble_screenshots(output_dir + '/' + output_prefix, output_dir + '/' + output_prefix + '_assemble.png', 128)
    else:
        print "no bmp copying from:", snapshots_dir
    return


def output_clusters(assign_ids, df_zscores, df_all, feature_names, output_dir, snapshots_dir=None):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    df_assign_id = pd.DataFrame()
    df_assign_id['specimen_name'] = df_all['specimen_name']
    df_assign_id['cluster_id'] = assign_ids
    df_assign_id.to_csv(output_dir + "/cluster_id.csv", index=False)

    clusters = np.unique(assign_ids)
    num_cluster = len(clusters)

    cluster_list = []  # for dunn index calculation
    print("There are %d clusters in total" % num_cluster)

    df_cluster = pd.DataFrame()
    df_zscore_cluster = pd.DataFrame()

    for i in clusters:
        ids = np.nonzero(assign_ids == i)[0]  # starting from  0
        df_cluster = df_all.iloc[ids]
        print("  %d neurons in cluster %d" % (df_cluster.shape[0], i))
        output_single_cluster_results(df_cluster, output_dir, "/cluster_" + str(i), snapshots_dir)

        df_zscore_cluster = df_zscores.iloc[ids]
        csv_file2 = output_dir + '/cluster_zscore_' + str(i) + '.csv'
        df_zscore_cluster.to_csv(csv_file2, index=False)

        cluster_list.append(df_zscore_cluster.values)

    ## pick the cluster specific feature and plot histogram
    cluster_specific_features(df_all, assign_ids, feature_names, output_dir+'/pvalues.csv')
    return cluster_list


####### ward  hierachichal clustering  ###########
def ward_cluster(df_all, feature_names, max_cluster_num, output_dir, snapshots_dir= None, RemoveOutliers = 0, datasetType='ivscc'):
  print("\n\n\n  ***************  ward computation, max_cluster = %d  *************:" % max_cluster_num)

  if not os.path.exists(output_dir):
    os.mkdir(output_dir)
  else:
    os.system("rm -r  " + output_dir + '/*')

  #### similarity plots
  # df_simMatrix = distance_matrix(df_all, feature_names, output_dir + "/morph_features_similarity_matrix.csv", 1)
  # # visualize heatmap using ward on similarity matrix
  # out = heatmap_plot_distancematrix(df_simMatrix, df_all, output_dir, "Similarity")
  # linkage = out.dendrogram_row.calculated_linkage


  ##### zscores  featuer plots
  df_zscores, df_all_outlier_removed, df_outliers = get_zscore_features(df_all, feature_names,
      output_dir + '/zscore.csv', RemoveOutliers)
  if (df_outliers.shape[0] > 0 ):
    output_single_cluster_results(df_outliers, output_dir, "outliers", snapshots_dir)

  if datasetType =='ivscc':
      linkage = heatmap_plot_zscore_ivscc(df_zscores, df_all_outlier_removed, output_dir, "feature zscores")
  if datasetType =='bbp':
      linkage = heatmap_plot_zscore_bbp(df_zscores, df_all_outlier_removed, output_dir, "feature zscores")

  assignments = hierarchy.fcluster(linkage, max_cluster_num, criterion="maxclust")
  #hierarchy.dendrogram(linkage)

  ## put assignments into ano files and csv files
  clusters_list = output_clusters(assignments, df_zscores, df_all_outlier_removed, feature_names, output_dir, snapshots_dir)
  dunn_index = dunn(clusters_list)
  print("dunn index is %f" % dunn_index)
  return linkage,df_zscores


def silhouette_clusternumber(linkage,df_zscores,output_dir ="."):
    #Silhouette analysis for determining the number of clusters

    print("Silhouettee analysis:")
    scores=[]
    for n_clusters in range(2,30):
         assignments = hierarchy.fcluster(linkage, n_clusters, criterion="maxclust")
         silhouette_avg = silhouette_score(df_zscores, assignments)
         print("For n_clusters =", n_clusters,"The average silhouette_score is :", silhouette_avg)
         scores.append(silhouette_avg)
    # plot sihouettee and cut
    pl.figure()
    pl.plot(range(2,30),scores,"*-")
    pl.xlabel("cluster number")
    pl.ylabel("average sihouettee coefficient")
    pl.savefig(output_dir+'/sihouettee_clusternumber.pdf')
    #pl.show()
    pl.close()
    return


def dunnindex_clusternumber(linkage,df_zscores, output_dir ="."):
     index_list=[]
     for n_clusters in range(2,30):
         assignments = hierarchy.fcluster(linkage, n_clusters, criterion="maxclust")
         df_assign_id = pd.DataFrame()

         df_assign_id['cluster_id'] = assignments

         clusters = np.unique(assignments)
         num_cluster = len(clusters)

         cluster_list = []  # for dunn index calculation

         df_cluster = pd.DataFrame()
         df_zscore_cluster = pd.DataFrame()
         for i in clusters:
            ids = np.nonzero(assignments == i)[0]  # starting from  0
            df_zscore_cluster = df_zscores.iloc[ids]
            cluster_list.append(df_zscore_cluster.values)

         dunn_index = dunn(cluster_list)
         index_list.append(dunn_index)
     pl.figure()
     pl.plot(range(2,30),index_list,"*-")
     pl.xlabel("cluster number")
     pl.ylabel("dunn index")
     pl.savefig(output_dir+'/dunnindex_clusternumber.pdf')
     #pl.show()
     return





def affinity_propagation(df_all, feature_names, output_dir, snapshots_dir=None, RemoveOutliers=0):
  ######  Affinity Propogation ##############

    print("\n\n\n ***************  affinity propogation computation ****************:")


    redundancy_removed_features_names = remove_correlated_features(df_all, feature_names, 0.95)
    print(" The %d features that are not closely correlated are %s" % (
        len(redundancy_removed_features_names), redundancy_removed_features_names))


    if not os.path.exists(output_dir):
      os.mkdir(output_dir)
    else:
      os.system("rm -r  " + output_dir + '/*')


        # Compute Affinity Propagation
    df_zscores, df_all_outlier_removed, df_outliers = get_zscore_features(df_all, redundancy_removed_features_names, None, RemoveOutliers)
    if (df_outliers.shape[0] > 0 ):
      output_single_cluster_results(df_outliers, output_dir, "outliers", snapshots_dir)

    X = df_zscores.as_matrix()

    af = AffinityPropagation().fit(X)
    cluster_centers_indices = af.cluster_centers_indices_
    labels = af.labels_
    labels = labels + 1  # the default labels start from 0, to be consistent with ward, add 1 so that it starts from 1
    clusters_list = output_clusters(labels, df_zscores, df_all_outlier_removed, redundancy_removed_features_names, output_dir,
        snapshots_dir)
    dunn_index = dunn(clusters_list)
    print("dunn index is %f" % dunn_index)
    return len(np.unique(labels)), dunn_index


def run_ward_cluster(df_features, feature_names, num_clusters,output_dir,output_postfix):

    redundancy_removed_features_names = remove_correlated_features(df_features, feature_names, 0.95)
    print(" The %d features that are not closely correlated are %s" % (
        len(redundancy_removed_features_names), redundancy_removed_features_names))

    #num_clusters, dunn_index1 = affinity_propagation(merged, redundancy_removed_features_names, output_dir + '/ap' + postfix, swc_screenshot_folder, REMOVE_OUTLIERS)
    linkage, df_zscore = ward_cluster(df_features, redundancy_removed_features_names, num_clusters, output_dir + '/ward' + output_postfix)
    silhouette_clusternumber(linkage, df_zscore, output_dir + '/ward' + output_postfix)
    return redundancy_removed_features_names


def main():

    ######################################################################################################################
    data_DIR = "/data/mat/xiaoxiaol/data/lims2/pw_aligned_1223"
    #default_all_feature_merged_file = data_DIR + '/keith_features_23dec.csv'

    #drop outliers, edit dendrite_type, creline
    #df_features = pd.read_csv(data_DIR +'/0107_new_features.csv')
    #df_features = df_features[df_features['QC status'] != "Outlier"]
    # #parse creline info from specimen_name
    #df_features.dropnas()
    # crelines=[]
    # swc_file_names=[]
    # for i in range(df_features.shape[0]):
    #     sn=df_features['specimen_name'][i]
    #     fn = df_features['specimen_name'][i].split('/')[-1]
    #     cl=sn.split(';')[0]
    #     crelines.append(cl)
    #     swc_file_names.append(fn)
    # df_features['cre_line'] = pd.Series(crelines)
    # df_features['swc_file_name'] = pd.Series(swc_file_names)
    # df_features.to_csv(data_DIR+'/filtered_w_cre.csv')





    input_csv_file = data_DIR + '/0108/0108_features.csv'
    out_dir = data_DIR + '/0108/clustering_results/no_GMI'
    default_swc_screenshot_folder =  data_DIR + "/figures/pw_aligned_bmps"
    #######################################################################################################################

    swc_screenshot_folder = default_swc_screenshot_folder

    method = "all"
    SEL_FEATURE = "all"


    if not os.path.exists(out_dir):
        os.mkdir(out_dir)


    ########################################################
    all_feature_file = input_csv_file
    #########################################################
    meta_feature_names = np.array(['specimen_name','specimen_id','dendrite_type','cre_line','region_info','filename','swc_file_name'])


    basal_feature_names = np.array(['basal_average_bifurcation_angle_local','basal_average_bifurcation_angle_remote','basal_average_contraction','basal_average_fragmentation',
                                    'basal_max_branch_order','basal_max_euclidean_distance','basal_max_path_distance',
                                     'basal_nodes_over_branches','basal_number_of_bifurcations',
                                    'basal_number_of_branches','basal_number_of_stems','basal_number_of_tips','basal_overall_depth','basal_overall_height',
                                     'basal_overall_width','basal_total_length','bb_first_moment_x_basal','bb_first_moment_y_basal','bb_first_moment_z_basal',
                                    'kg_soma_depth',
                                     'basal_moment1','basal_moment10','basal_moment11','basal_moment12','basal_moment13','basal_moment2',
                                      'basal_moment3','basal_moment4',
                                     'basal_moment5','basal_moment6','basal_moment7','basal_moment8','basal_moment9'])
                                      #'basal_total_surface','basal_total_volume','basal_soma_surface','basal_number_of_nodes','basal_average_diameter',
                                      # 'basal_moment1','basal_moment10','basal_moment11','basal_moment12','basal_moment13','basal_moment14','basal_moment2','basal_moment3','basal_moment4',
        #'basal_moment5','basal_moment6','basal_moment7','basal_moment8','basal_moment9','basal_average_parent_daughter_ratio'


    apical_feature_names = np.array(['apical_average_bifurcation_angle_local','apical_average_bifurcation_angle_remote','apical_average_contraction',
                                     'apical_average_fragmentation','apical_max_branch_order','apical_max_euclidean_distance',
                                     'apical_max_path_distance',
                                     'apical_nodes_over_branches','apical_number_of_bifurcations','apical_number_of_branches',
                                     'apical_number_of_tips','apical_overall_depth','apical_overall_height','apical_overall_width','apical_total_length',
                                     'kg_branch_mean_from_centroid_z_apical',
                                           'kg_branch_stdev_from_centroid_z_apical',
                                           'kg_centroid_over_farthest_branch_apical',
                                           'kg_centroid_over_farthest_neurite_apical',
                                           'kg_centroid_over_radial_dist_apical',
                                           'kg_mean_over_centroid',
                                           'kg_mean_over_farthest_branch_apical',
                                           'kg_mean_over_farthest_neurite_apical',
                                           'kg_mean_over_radial_dist_apical',
                                           'kg_mean_over_stdev',
                                           'kg_num_branches_over_radial_dist_apical',
                                           'kg_num_outer_apical_branches',
                                           'kg_outer_mean_from_center_z_apical',
                                           'kg_outer_mean_over_stdev',
                                           'kg_outer_stdev_from_center_z_apical',
                                           'kg_peak_over_moment_z_apical',
                                           'kg_radial_dist_over_moment_z_apical',
                                           'kg_soma_depth'])
                                   #, 'apical_number_of_nodes'
                                    # ])#'apical_soma_surface', 'apical_total_surface','apical_total_volume','apical_average_diameter','apical_moment1','apical_moment10','apical_moment11','apical_moment12','apical_moment13','apical_moment14',
                                     # 'apical_moment2','apical_moment3','apical_moment4','apical_moment5','apical_moment6','apical_moment7','apical_moment8','apical_moment9','apical_average_parent_daughter_ratio','apical_number_of_stems?? always 1',

    bbp_feature_names = np.array(['bb_first_moment_apical','bb_first_moment_basal','bb_first_moment_dendrite','bb_first_moment_x_apical','bb_first_moment_x_basal',
                                 'bb_first_moment_x_dendrite','bb_first_moment_y_apical','bb_first_moment_y_basal','bb_first_moment_y_dendrite','bb_first_moment_z_apical',
                                 'bb_first_moment_z_basal','bb_first_moment_z_dendrite','bb_max_branch_order_apical','bb_max_branch_order_basal','bb_max_branch_order_dendrite',
                                 'bb_max_path_length_apical','bb_max_path_length_basal','bb_max_path_length_dendrite','bb_max_radial_distance_apical','bb_max_radial_distance_basal',
                                 'bb_max_radial_distance_dendrite','bb_mean_trunk_diameter_apical','bb_mean_trunk_diameter_basal','bb_mean_trunk_diameter_dendrite',
                                 'bb_number_branches_apical','bb_number_branches_basal','bb_number_branches_dendrite','bb_number_neurites_apical','bb_number_neurites_basal',
                                 'bb_number_neurites_dendrite','bb_second_moment_apical','bb_second_moment_basal','bb_second_moment_dendrite','bb_second_moment_x_apical',
                                 'bb_second_moment_x_basal','bb_second_moment_x_dendrite','bb_second_moment_y_apical','bb_second_moment_y_basal','bb_second_moment_y_dendrite',
                                 'bb_second_moment_z_apical','bb_second_moment_z_basal','bb_second_moment_z_dendrite','bb_total_length_apical','bb_total_length_basal',
                                 'bb_total_length_dendrite'])
                                 #'bb_total_surface_area_apical','bb_total_volume_basal','bb_total_volume_apical','bb_total_volume_dendrite','bb_total_surface_area_basal','bb_total_surface_area_dendrite'




    #selected_features = ['max_euclidean_distance', 'num_stems', 'num_bifurcations', 'average_contraction',
    #'parent_daughter_ratio']

    #tmp = np.append(meta_feature_names, basal_feature_names)
    all_dendritic_feature_names =  np.append(basal_feature_names, apical_feature_names)  #bbp_feature_names
    spiny_feature_names =  apical_feature_names
    aspiny_feature_names = basal_feature_names

    df_features = pd.read_csv(all_feature_file)

    print df_features.columns
    df_features[all_dendritic_feature_names]= df_features[all_dendritic_feature_names].astype(float)
    print "There are %d neurons in this dataset" % df_features.shape[0]
    print "Dendrite types: ", np.unique(df_features['dendrite_type'])


    # df_features_all = df_features[np.append(meta_feature_names,all_dendritic_feature_names)]
    # df_features_all.to_csv(data_DIR+'/0108/all_dendrite_features.csv')

    df_groups = df_features.groupby(['dendrite_type'])

    df_spiny = df_groups.get_group('spiny')
    # df_w_spiny = df_spiny[np.append(meta_feature_names,spiny_feature_names)]
    # df_w_spiny.to_csv(data_DIR +'/0108/spiny_features.csv', index=False)


    df_aspiny =  pd.concat([df_groups.get_group('aspiny'),df_groups.get_group('sparsely spiny')],axis=0)
    # df_w_aspiny = df_aspiny[np.append(meta_feature_names,aspiny_feature_names)]
    # df_w_aspiny.to_csv(data_DIR +'/0108/aspiny_features.csv', index=False)

    print "There are %d neurons are aspiny " % df_aspiny.shape[0]

    print "There are %d neurons are spiny\n\n" % df_spiny.shape[0]


    feature_names = all_dendritic_feature_names
    method = "ward"

    REMOVE_OUTLIERS = 0
    postfix = "_" + SEL_FEATURE
    postfix += "_ol_clipped"

    #run_ward_cluster(df_features, feature_names, num_clusters,output_postfix):


    # num_clusters, dunn_index1 = affinity_propagation(df_aspiny, aspiny_feature_names,
    #                                                      out_dir + '/ap_aspiny' + postfix,
    #                                                     None, REMOVE_OUTLIERS)
    # print "spiny ap:"
    # print num_clusters
    #
    # num_clusters, dunn_index1 = affinity_propagation(df_spiny, spiny_feature_names,
    #                                                      out_dir + '/ap_spiny' + postfix,
    #                                                     None, REMOVE_OUTLIERS)
    # print "aspiny ap:"
    # print num_clusters
    # exit()
    redundancy_removed_features = run_ward_cluster(df_aspiny, aspiny_feature_names, num_clusters = 6 ,output_dir = out_dir,output_postfix= '_aspiny'+postfix)
    df_w_aspiny = df_aspiny[np.append(meta_feature_names,redundancy_removed_features)]
    df_w_aspiny.to_csv(data_DIR +'/0108/aspiny_selected_features.csv', index=False)
    # #


    # df_spiny.fillna(0,inplace=True)
    # redundancy_removed_features = run_ward_cluster(df_spiny, spiny_feature_names, num_clusters = 9 ,output_dir = out_dir, output_postfix='_spiny'+ postfix)
    # df_w_spiny = df_spiny[np.append(meta_feature_names,redundancy_removed_features)]
    # df_w_spiny.to_csv(data_DIR +'/0108/spiny_selected_features.csv', index=False)
    #

if __name__ == "__main__":
        main()
