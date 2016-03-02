__author__ = 'xiaoxiaoliu'

import pandas as pd
import numpy as np

import numpy as np
from feature_names import *
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

def cluster_specific_features(df_all, assign_ids, feature_names, output_csv_fn):
    #student t to get cluster specific features

    labels=[]
    clusters = np.unique(assign_ids)
    num_cluster = len(clusters)
    df_pvalues =  pd.DataFrame(index = feature_names, columns = clusters)

    for cluster_id in clusters:
        print cluster_id
        ids_a = np.nonzero(assign_ids == cluster_id)[0]  # starting from  0
        ids_b = np.nonzero(assign_ids != cluster_id)[0]  # starting from  0
        labels.append("C"+str(cluster_id) + "("+ str(len(ids_a))+")" )
        for feature in feature_names:
            a = df_all.iloc[ids_a][feature]
            b = df_all.iloc[ids_b][feature]

            t_stats,pval = stats.ttest_ind(a,b,equal_var=False)
            print pval
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
    pl.subplots_adjust(left=0.3, right=0.9, top=0.9, bottom=0.1)
    pl.title('-log10(P value)')
    filename = output_csv_fn + '.png'
    pl.savefig(filename, dpi=300)
    pl.show()
    pl.close()


    return df_pvalues





def generateLinkerFileFromCSV(result_dir, csvfile, column_name=None, strip_path=True, fpath='.'):
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
    print np.unique(types)
    for atype in np.unique(types):
        idxs = np.nonzero(types == atype)[0]
        swc_files = df['swc_file_name']
        with open(result_dir + '/' + atype + '.ano', 'w') as outf:
            for afile in swc_files[idxs]:
                filename = afile
                if strip_path:
                    filename = afile.split('/')[-1]
                line = 'SWCFILE=' + fpath +"/"+filename + '\n'
                outf.write(line)
            outf.close()
    return


#data_DIR = "/data/mat/xiaoxiaol/data/BBP/raw_converted_swcs"
data_DIR="/data/mat/xiaoxiaol/data/BBP/BBP_data_uncurated_SWC"

if 0: # figure out the meta info
    #build map
    ref=pd.read_csv(data_DIR+"/meta_refer.csv")
    keys = ref['shape_type']
    values1 = ref['dendrite_type']
    values2 = ref['detailed_type_name']
    print values2
        #print np.unique(df_check_table.algorithm)
    checkDendritetype = dict(zip(keys, values1))
    checkDetailedName = dict(zip(keys, values2))



    df_meta=pd.read_csv(data_DIR+"/morph_type_upper.csv")
    dendrite_types=[]
    detailed_type_names=[]
    for i in range(df_meta.shape[0]):
        st= df_meta.iloc[i]['shape_type']
        dy=np.nan
        if checkDendritetype.has_key(st):
            dy =checkDendritetype[st]
        dendrite_types.append(dy)

        dt=st
        if checkDetailedName.has_key(st):
              dt=checkDetailedName[st]
        detailed_type_names.append(dt)

    df_meta['dendrite_type']=pd.Series(dendrite_types)
    df_meta['detailed_type_name']=pd.Series(detailed_type_names)
    df_meta.to_csv(data_DIR+'/rich_meta_type.csv')


####  process meta info
if 0:
    csv_file = data_DIR+"/rich_meta_type.csv"
    df_bbp = pd.read_csv(csv_file)
    for i in range (df_bbp.shape[0]):
        df_bbp.ix[i,'layer'] = df_bbp.iloc[i]['m-type'].split('_')[0]


    #generateLinkerFileFromCSV(data_DIR, csv_file, column_name='m-type',strip_path=False,
     #                         fpath=data_DIR+'/sorted')
## process feature csv
if 0:
     df_f =pd.read_csv (data_DIR+"/sorted/morph_features.csv")

     df_f= df_f.drop(['specimen_name', 'specimen_id', 'dendrite_type', 'region_info','ignore'],axis=1)

     for i in range (df_f.shape[0]):
        tmp = df_f.iloc[i]['filename'].split('/')[-1]
        df_f.ix[i,'swc_file_name'] = (tmp.split('.')[0]).upper() +'.swc'


     df_f.to_csv(data_DIR+'/morph_features_fixed.csv', index=False)
     print df_f.shape
     print df_bbp.shape
     df_2=pd.merge(df_f, df_bbp, on="swc_file_name")
     print df_2.shape
     df_2.to_csv(data_DIR+'/morph_features_with_meta_tags.csv', index=False)
     col_names =[]
     col_names.extend(meta_names)
     col_names.extend(ALL_FEATURES)
     df_out= df_2[col_names]
     df_out.to_csv(data_DIR+'/morph_features_with_meta_tags_all_final.csv', index=False)

     dfg=df_out.groupby(['dendrite_type'])

     df_e =dfg.get_group('excitatory')
     df_e = df_e[EXCITATORY_cols]
     df_e.to_csv(data_DIR+'/excitatory_features.csv', index=False)

     df_i =dfg.get_group('inhibitory')
     df_i=df_i[INHIBITORY_cols]
     df_i.to_csv(data_DIR+'/inhibitory_features.csv', index=False)

if 1:
    df_e=pd.read_csv(data_DIR+'/excitatory_features.csv')
    df_i=pd.read_csv(data_DIR+'/inhibitory_features.csv')
    ###  feature specific
    #excitatory
    assign_ids=np.array(df_e['m-type'])

    cluster_specific_features(df_e, assign_ids, APICAL_features, data_DIR+'/exci_pvalues.csv')
    print done

    assign_ids=np.array(df_i['m-type'])
    cluster_specific_features(df_i, assign_ids, AXON_features, data_DIR+'/inh_pvalues.csv')