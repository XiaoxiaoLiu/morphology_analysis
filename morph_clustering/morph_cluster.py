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
import numpy.random
from itertools import cycle
import gap_analysis as gap


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
def distance_matrix(df_all, feature_names, out_distanceMatrix_file, sample_size = 0,REMOVE_OUTLIER=0):
    feature_array = df_all[feature_names].astype(float)
    distanceMatrix = []
    normalized = zscore(feature_array)


    print "feature matrix size is:"
    print normalized.shape
    #normalized = normalizeFeatures(feature_array)


    if not REMOVE_OUTLIER:  # only clp
            normalized[normalized < -ZSCORE_OUTLIER_THRESHOLD] = -ZSCORE_OUTLIER_THRESHOLD
            normalized[normalized > ZSCORE_OUTLIER_THRESHOLD] = ZSCORE_OUTLIER_THRESHOLD
    if sample_size>0:
         random_ids =  numpy.random.random_integers(0,normalized.shape[0]-1, sample_size)
         normalized = normalized[random_ids]

    for i in range(len(normalized)):
        queryFeature = normalized[i]  # each row is a feature vector

        #scores = np.exp(-np.sum(abs(normalized - queryFeature) ** 2, 1) / 100)  #similarity
        scores = np.sum(np.abs(normalized - queryFeature) ** 2, 1)  # distance
        distanceMatrix.append(scores)
    df_dist = pd.DataFrame(distanceMatrix)
    print "distance similarity matrix size is:"
    print df_dist.shape
    df_dist.to_csv(out_distanceMatrix_file, index=False)
    print("score sim matrix is saved to : " + out_distanceMatrix_file + "\n")
    return df_dist










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
def heatmap_plot_distancematrix(df_distanceMatrix, output_dir, title=None):
    pl.figure()



    g = sns.heatmap(df_distanceMatrix, method='ward', metric='euclidean', linewidths=0.0, xticklabels=False,
                       yticklabels=False)
    if title:
        pl.title(title)


    pl.title('Similarities')

    filename = output_dir + '/similarity_heatmap.png'
    pl.savefig(filename, dpi=300)
    print("save similarity matrix heatmap figure to :" + filename)
    pl.close()



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
    print "heatmap plot:ivscc"
    sns.set_context("talk", font_scale=1.4)
    # Create a custom palette for dendrite_type colors
    dendrite_types = [np.nan, 'aspiny', 'sparsely spiny', 'spiny']
    # dendrite_type_pal = sns.color_palette("coolwarm", len(dendrite_types))
    dendrite_type_pal = sns.color_palette(["gray","black","purple","red"])
    dendrite_type_lut = dict(zip(dendrite_types, dendrite_type_pal))
    dendrite_type_colors = df_all['dendrite_type'].map(dendrite_type_lut)


    # Create a custom palette for creline colors
    cre_lines = np.unique(df_all['cre_line'])
    #print cre_lines
    #cre_lines = ['Pvalb-IRES-Cre','Sst-IRES-Cre','Gad2-IRES-Cre', 'Htr3a-Cre_NO152',
    #             'Nr5a1-Cre', 'Ntsr1-Cre','Rbp4-Cre_KL100' ,'Rorb-IRES2-Cre-D', 'Scnn1a-Tg2-Cre',
    #             'Scnn1a-Tg3-Cre','Slc17a6-IRES-Cre','Cux2-CreERT2']

    cre_line_pal = sns.color_palette("BrBG", len(cre_lines))

    cre_line_lut = dict(zip(cre_lines, cre_line_pal))  # map creline type to color
    cre_line_colors = df_all['cre_line'].map(cre_line_lut)


    layers = np.unique(df_all['layer'])
    layer_pal = sns.light_palette("black", len(layers))
    layer_lut = dict(zip(layers, layer_pal))
    layer_colors = df_all['layer'].map(layer_lut)

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
    print "linkage genearted"
    data = df_zscore_features.transpose()
    row_linkage = hierarchy.linkage(data, method='ward', metric='euclidean')
    feature_order = hierarchy.leaves_list(row_linkage)


    matchIndex = [data.index[x] for x in feature_order]

    data = data.reindex(matchIndex)

    print "plot heatmap"
    g = sns.clustermap(data, row_cluster = False, col_linkage=linkage, method='ward', metric='euclidean',
                       linewidths = 0.0,col_colors = [cre_line_colors,layer_colors],
                       cmap = sns.cubehelix_palette(light=1, as_cmap=True),figsize=(60,15))

    pl.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=90 )
    pl.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    pl.subplots_adjust(left=0.1, bottom=0.5, right=0.8, top=0.95)  # !!!!!

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

    for label in layers:
        g.ax_row_dendrogram.bar(0, 0, color=layer_lut[label], label=label, linewidth=0)
        g.ax_row_dendrogram.legend(loc=location, ncol=1,borderpad=0.0)

    #
    # for label in types:
    #      g.ax_row_dendrogram.bar(0, 0, color=type_lut[label], label=label,linewidth=0)
    #      g.ax_row_dendrogram.legend(loc=location, ncol=num_cols,borderpad=0.0)
    #
    #
    # g.ax_row_dendrogram.bar(0, 0, color = "white", label=" ", linewidth=0)
    # g.ax_row_dendrogram.legend(loc=location, ncol=num_cols, borderpad=0.0)


    # for label in dendrite_types:
    #     g.ax_row_dendrogram.bar(0, 0, color = dendrite_type_lut[label], label=label, linewidth=0)
    #     g.ax_row_dendrogram.legend(loc=location, ncol= num_cols, borderpad=0.0)


    filename = output_dir + '/zscore_feature_heatmap.png'
    pl.savefig(filename, dpi=300)
    #pl.show()
    print("save zscore matrix heatmap figure to :" + filename)
    pl.close()
    return linkage


def heatmap_plot_zscore_bigneuron(df_zscore_features, df_all, output_dir, title=None):

    print "heatmap plot:bigneuron"

    #taiwan
    metric ='nt_type'
    mtypes = np.unique(df_all[metric])
    print mtypes
    mtypes_pal = sns.color_palette("hls", len(mtypes))

    mtypes_lut = dict(zip(mtypes, mtypes_pal))
    mtypes_colors = df_all[metric].map(mtypes_lut)



    linkage = hierarchy.linkage(df_zscore_features, method='ward', metric='euclidean')

    data = df_zscore_features.transpose()
    row_linkage = hierarchy.linkage(data, method='ward', metric='euclidean')
    feature_order = hierarchy.leaves_list(row_linkage)

    #print data.index
    matchIndex = [data.index[x] for x in feature_order]
    #print matchIndex
    data = data.reindex(matchIndex)

    pl.figure()
    g = sns.clustermap(data, row_cluster = False, col_linkage=linkage, method='ward', metric='euclidean',
                       linewidths = 0.0,col_colors = [mtypes_colors],
                       cmap = sns.cubehelix_palette(light=1, as_cmap=True),figsize=(40,10))

    pl.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    pl.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
    #g.ax_heatmap.set_xticklabels([])
    pl.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.95)  # !!!!!

    if title:
        pl.title(title)


    location ="best"
    num_cols=1
    # Legend for row and col colors

    for label in mtypes:
         g.ax_row_dendrogram.bar(0, 0, color=mtypes_lut[label], label=label, linewidth=0.0)
         g.ax_row_dendrogram.legend(loc=location, ncol=num_cols,borderpad=0)

    filename = output_dir + '/zscore_feature_heatmap.png'
    pl.savefig(filename, dpi=300)
    #pl.show()
    print("save zscore matrix heatmap figure to :" + filename)
    pl.close()
    print "done clustering and heatmap plotting"
    return linkage



def heatmap_plot_zscore_bbp(df_zscore_features, df_all, output_dir, title=None):

    print "heatmap plot:bbp"
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

    subset_features_names = feature_names#.tolist()
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
    pl.subplots_adjust(left=0.5, right=0.95, top=0.95, bottom=0.3)
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

    df_all_modified = df_all  # with original feature value
    df_outliers = pd.DataFrame()

    if num_outliers > 0:
        if not REMOVE_OUTLIER:  # just clip
            normalized[normalized < -zscore_threshold] = -zscore_threshold
            normalized[normalized > zscore_threshold] = zscore_threshold

        else:
            outliers_l = np.nonzero(normalized < -zscore_threshold)
            outliers_h = np.nonzero(normalized > zscore_threshold)
            outlier_index = np.unique((np.append(outliers_l[0], outliers_h[0])))

            # remove outlier rows
            df_all_modified = df_all_modified.drop(df_all_modified.index[outlier_index])
            normalized = np.delete(normalized, outlier_index, 0)

            # re-zscoring and clipping
            # m_featureArray = df_all_modified[feature_names].astype(float)
            # normalized = zscore(m_featureArray)
            # normalized[normalized < -zscore_threshold] = -zscore_threshold
            # normalized[normalized > zscore_threshold] = zscore_threshold

            print("Removed %d outlier neurons" % len(outlier_index))

            df_outliers = df_all.iloc[outlier_index]


    df_z = pd.DataFrame(normalized)
    df_z.columns = feature_names
    df_z.index = df_all_modified['swc_file_name']

    if out_file:
        df_z.to_csv(out_file, index=True)
        df_all_modified.to_csv(out_file[0:-10]+'feature_modified.csv', index=True)
        print("save to " + out_file )


    if (df_z.shape[0] != df_all_modified.shape[0]):
        print ("error:  the sample size of the zscore and the original table does not match!")

    df_outliers.to_csv(out_file[0:-10]+'zscore_outliers.csv', index=True)

    return df_z, df_all_modified, df_outliers


#############################################################################################
def output_single_cluster_results(df_cluster, output_dir, output_prefix,  swc_path = None):
    csv_file = output_dir + '/' + output_prefix + '.csv'
    df_cluster.to_csv(csv_file, index=False)

    ano_file = output_dir + '/' + output_prefix + '.ano'
    generateLinkerFileFromDF(df_cluster, ano_file, False, swc_path)
    # copy bmp vaa3d snapshots images over

    return


def output_clusters(assign_ids, df_zscores, df_all, feature_names, output_dir, swc_path=None):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    df_assign_id = pd.DataFrame(columns=['specimen_name', 'specimen_id', 'cluster_id'])

    if 'specimen_name' in df_all.columns:
         df_assign_id['specimen_name'] = df_all['specimen_name']
    else:
        print "missing specimen_name:",df_all.columns
    if 'specimen_id' in df_all.columns:
         df_assign_id['specimen_id'] = df_all['specimen_id']
    else:
        print "missing specimen_id:",df_all.columns
    df_assign_id['cluster_id'] = assign_ids
    df_assign_id.to_csv(output_dir + "/cluster_id.csv", index=False)

    clusters = np.unique(assign_ids)
    num_cluster = len(clusters)

    cluster_list = []  # for dunn index calculation
    print("There are %d clusters in total" % num_cluster)

    df_cluster = pd.DataFrame()
    df_zscore_cluster = pd.DataFrame()

     #pick the examplar and save into ano file
    print "Save cluster examplar to  exmplar.ano."
    with open(output_dir+'/examplar.ano', 'w') as outf:
      for i in clusters:
        ids = np.nonzero(assign_ids == i)[0]  # starting from  0
        df_cluster = df_all.iloc[ids]
        #print("  %d neurons in cluster %d" % (df_cluster.shape[0], i))
        output_single_cluster_results(df_cluster, output_dir, "/cluster_" + str(i), swc_path)

        df_zscore_cluster = df_zscores.iloc[ids]
        csv_file2 = output_dir + '/cluster_zscore_' + str(i) + '.csv'
        df_zscore_cluster.to_csv(csv_file2, index=False)

        cluster_list.append(df_zscore_cluster.values)

        afile=df_cluster.iloc[len(df_cluster)/2]['swc_file_name']
        #print afile
        filename = swc_path + '/'+afile
        line = 'SWCFILE=' + filename + '\n'
        outf.write(line)

    ## pick the cluster specific feature and plot histogram
    cluster_specific_features(df_all, assign_ids, feature_names, output_dir+'/pvalues.csv')


    outf.close()


    return cluster_list


####### ward  hierachichal clustering  ###########

def fancy_dendrogram(*args, **kwargs):
    max_d = kwargs.pop('max_d', None)
    if max_d and 'color_threshold' not in kwargs:
        kwargs['color_threshold'] = max_d
    annotate_above = kwargs.pop('annotate_above', 0)

    ddata = hierarchy.dendrogram(*args, **kwargs)

    if not kwargs.get('no_plot', False):
        pl.title('Hierarchical Clustering Dendrogram (truncated)')
        pl.xlabel('sample index or (cluster size)')
        pl.ylabel('distance')
        for i, d, c in zip(ddata['icoord'], ddata['dcoord'], ddata['color_list']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            if y > annotate_above:
                pl.plot(x, y, 'o', c=c)
                pl.annotate("%.3g" % y, (x, y), xytext=(0, -5),
                             textcoords='offset points',
                             va='top', ha='center')
        if max_d:
            pl.axhline(y=max_d, c='k')
    return ddata

def truncate_dendrogram(linkage, max_cluster_num,output_dir, max_d=0 ):

    pl.figure()
    pl.title('Hierarchical Clustering Dendrogram (truncated)')
    pl.xlabel('sample index')
    pl.ylabel('distance')

    fancy_dendrogram(
        linkage,
        truncate_mode='lastp',
        p=max_cluster_num,
        leaf_rotation=90.,
        leaf_font_size=12.,
        show_contracted=True,
        annotate_above=10,
        max_d=max_d,  # plot a horizontal cut-off line
    )

    pl.savefig(output_dir+'/dendrogram.png',dpi=300)
    #pl.show()
    return



def ward_cluster(df_all, feature_names, max_cluster_num, output_dir, swc_path= None, RemoveOutliers = 0,
                 datasetType='ivscc', plot_heatmap =1):
  print("\n\n\n  ***************  ward computation, max_cluster = %d  *************:" % max_cluster_num)

  if not os.path.exists(output_dir):
    os.mkdir(output_dir)

  ##### zscores  featuer plots
  df_zscores, df_all_outlier_removed, df_outliers = get_zscore_features(df_all, feature_names,
      output_dir + '/zscore.csv', RemoveOutliers)
  if (df_outliers.shape[0] > 0 ):
    output_single_cluster_results(df_outliers, output_dir, "outliers", swc_path)

  if plot_heatmap:
      if datasetType =='ivscc':
          link = heatmap_plot_zscore_ivscc(df_zscores, df_all_outlier_removed, output_dir, "feature zscores")
      if datasetType =='bbp':
          link = heatmap_plot_zscore_bbp(df_zscores, df_all_outlier_removed, output_dir, "feature zscores")
      if datasetType =='bigneuron':
          link = heatmap_plot_zscore_bigneuron(df_zscores, df_all_outlier_removed, output_dir, "feature zscores")
  else:
      link = hierarchy.linkage(df_zscores, method='ward', metric='euclidean')


  assignments = hierarchy.fcluster(link, max_cluster_num, criterion="maxclust")
  output_clusters(assignments, df_zscores, df_all_outlier_removed, feature_names, output_dir, swc_path)

  truncate_dendrogram(link,max_cluster_num,output_dir,0)

  return link,df_zscores





def silhouette_clusternumber(linkage,df_zscores,low=1, high=5,output_dir ="."):
    #Silhouette analysis for determining the number of clusters

    print("Silhouettee analysis:")
    scores=[]
    step = 1
    if step < (high-low)/20:
         step = (high-low)/20
    for n_clusters in range(low,high,step):
         assignments = hierarchy.fcluster(linkage, n_clusters, criterion="maxclust")
         silhouette_avg = silhouette_score(df_zscores, assignments,metric='sqeuclidean')
         print("For n_clusters =", n_clusters,"The average silhouette coefficient is :", silhouette_avg)
         scores.append(silhouette_avg)
    # plot sihouettee and cut
    pl.figure()
    pl.plot(range(low,high),scores,"*-")
    pl.xlabel("cluster number")
    pl.ylabel("average sihouettee coefficient")
    pl.savefig(output_dir+'/sihouettee_clusternumber.png')
    #pl.show()
    pl.close()
    return


def dunnindex_clusternumber(linkage,df_zscores, low=1, high=5,output_dir ="."):
     index_list=[]
     for n_clusters in range(low,high):

         assignments = hierarchy.fcluster(linkage, n_clusters, criterion="maxclust")
         df_assign_id = pd.DataFrame()

         df_assign_id['cluster_id'] = assignments
         clusters = np.unique(assignments)
         cluster_list = []  # for dunn index calculation

         for i in clusters:
            ids = np.nonzero(assignments == i)[0]  # starting from  0
            df_zscore_cluster = df_zscores.iloc[ids]
            cluster_list.append(df_zscore_cluster.values)

         dunn_index = dunn(cluster_list)
         print n_clusters, ":", dunn_index
         index_list.append(dunn_index)
     pl.figure()
     pl.plot(range(low,high),index_list,"*-")
     pl.xlabel("cluster number")
     pl.ylabel("dunn index")
     pl.savefig(output_dir+'/dunnindex_clusternumber.png')
     #pl.show()
     return





def affinity_propagation(df_all, feature_names, output_dir, swc_path=None, RemoveOutliers=0):
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
    df_zscores, df_all_outlier_removed, df_outliers = get_zscore_features(df_all, redundancy_removed_features_names,
    output_dir+'/zscore.csv', RemoveOutliers)
    if (df_outliers.shape[0] > 0 ):
      output_single_cluster_results(df_outliers, output_dir, "outliers", swc_path)
    print "There are ", df_zscores.shape[0],"neurons entering the clustering"

    X = df_zscores.as_matrix()

    af = AffinityPropagation().fit(X)
    cluster_centers_indices = af.cluster_centers_indices_
    labels = af.labels_
    print("Silhouette Coefficient: %0.3f"
      % silhouette_score(df_zscores, labels, metric='sqeuclidean'))
    labels = labels + 1  # the default labels start from 0, to be consistent with ward, add 1 so that it starts from 1
    clusters_list = output_clusters(labels, df_zscores, df_all_outlier_removed, redundancy_removed_features_names, output_dir,
        swc_path)


    #dunn_index = dunn(clusters_list)
    #print("dunn index is %f" % dunn_index)
    return len(np.unique(labels))



import pickle
def run_ward_cluster(df_features, feature_names, num_clusters,output_dir,output_postfix,experiment_type='ivscc', low=2, high=5,
                     plot_heatmap=1,RemoveOutliers=0, swc_path=None, rerun=1):
    #experiment type: ivscc, bbp, bigneuron
    redundancy_removed_features_names = remove_correlated_features(df_features, feature_names, 0.95)
    print(" The %d features that are not closely correlated are %s" % (
        len(redundancy_removed_features_names), redundancy_removed_features_names))
  
    out_dir = output_dir+'/ward' + output_postfix
    if os.path.exists(out_dir + "/linkage.pik") and rerun ==0:
          pickle.dump([linkage], f, -1)
          print "load existing linkage file from disk:"
          linkage = pickle.load( open( out_dir+"/linkage.pik", "rb" ) )
    else:
	  linkage, df_zscore = ward_cluster(df_features, redundancy_removed_features_names, num_clusters, out_dir, swc_path, RemoveOutliers, experiment_type, plot_heatmap)
          with open(output_dir+'/linkage.pik', 'wb') as f:
               pickle.dump([linkage], f, -1)
               print "Save linkage to :", output_dir+'/linkage.pik'

    print "silhouette analysis"
    silhouette_clusternumber(linkage, df_zscore, low,high,output_dir + '/ward' + output_postfix)

    print "\n\n\n\ngap analysis"
    gap.gap_analysis(df_zscore.values,low,high,10,output_dir + '/ward' + output_postfix)


    #print "dunn index:"
    #dunnindex_clusternumber(linkage, df_zscore, low,high,output_dir + '/ward' + output_postfix)

    ### similarity plots
    #visualize heatmap using ward on similarity matrix
    #( randomly select 500 from the dataset)
    #sample_size = 500
    #df_simMatrix = distance_matrix(df_features, redundancy_removed_features_names, output_dir + "/"+str(sample_size)+"_morph_features_similarity_matrix.csv",sample_size , 1)

    return redundancy_removed_features_names




def run_affinity_propagation(df_features, feature_names,output_dir,output_postfix, swc_path=None,outlier_removal=0):
    redundancy_removed_features_names = remove_correlated_features(df_features, feature_names, 0.95)
    print(" The %d features that are not closely correlated are %s" % (
        len(redundancy_removed_features_names), redundancy_removed_features_names))

    num_clusters  = affinity_propagation(df_features, redundancy_removed_features_names, output_dir + '/ap' + output_postfix, swc_path, outlier_removal)
    print "ap clusters:",num_clusters
    return redundancy_removed_features_names
