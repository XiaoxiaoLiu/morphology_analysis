import numpy as np
import matplotlib.pylab as pl
import matplotlib.pyplot as pplt

from scipy import stats

import pandas as pd
from sklearn.lda import LDA
from sklearn.decomposition import PCA
#from sklearn.lda import PCA
import scipy.cluster.hierarchy as hac
from mpl_toolkits.mplot3d import Axes3D


from os import sys, path
p = path.dirname(path.dirname(path.abspath(__file__)))
sys.path.append(p)
sys.path.append(p+'/utilities')

import lims_orca_utils as lu

# program path on this machine
#===================================================================
WORK_PATH="/Users/xiaoxiaoliu/work"
#WORK_PATH="/home/xiaoxiaol/work"

MRMR= WORK_PATH+"/src/mrmr_c_src/mrmr"
V3D="qv3d"

# data dir
data_DIR= WORK_PATH+"/data/lims2/neuron_recon_2"

list_csv_file =  data_DIR+'/list.csv'
feature_csv_file = data_DIR+"/glFeatures.csv"
data_linker_file =  data_DIR+'/original/mylinker.ano'
preprocessed_data_linker_file = data_DIR+'/preprocessed/mylinker.ano'
feature_file = data_DIR + '/preprocessed/prep_features.nfb'
feature_names= np.array(['num_nodes', 'soma_surface', 'num_stems','num_bifurcations', 'num_branches', 'num_of_tips',  'overall_width', 'overall_height',  'overall_depth', 'average_diameter',    'total_length', 'total_surface', 'total_volume', 'max_euclidean_distance',       'max_path_distance', 'max_branch_order',  'average_contraction', 'average fragmentation', 'parent_daughter_ratio', 'bifurcation_angle_local', 'bifurcation_angle_remote'])
selected_features=['max_euclidean_distance','num_stems','num_bifurcations','average_contraction','parent_daughter_ratio']


#===================================================================
def zscore(features):
    zscores = stats.zscore(features,0)
    return zscores


def np(features):
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


def plotScoreMatrix(featureArray, fig_title):
    scoreMatrix =[]
    normalized = normalizeFeatures(featureArray)

    for i in range(len(normalized)):
        queryFeature = normalized[i] # each row is a feature vecto
        scores = np.exp(-np.sum(abs(normalized-queryFeature)**2,1)/100)
        #scores = np.sum(abs(normalized-queryFeature)**2,1)
        scoreMatrix.append(scores)

    pl.figure()
    pl.imshow(scoreMatrix,interpolation='none')
    pl.colorbar()
    pl.title(fig_title)
    pl.xlabel('neuron ID')
    pl.ylabel('neuron ID')
    pl.show()
    pl.savefig(data_DIR+'/'+fig_title+'.png')


def compareHist(data1, data2,_title,tag1='data1', tag2='data2'):
    pl.figure()
    pl.show()
    pl.hist(data1, normed=True, alpha=0.5, color='b')
    pl.hist(data2, normed=True, alpha=0.5, color='r')

    # Fit a normal distribution to the data:
    mu1, std1 = stats.norm.fit(data1)
    xmin, xmax = pl.xlim()
    x = np.linspace(xmin, xmax, 100)
    p = stats.norm.pdf(x, mu1, std1)
    pl.plot(x, p, 'k', linewidth=2, color='b')

    # Fit a normal distribution to the data:
    mu2, std2 = stats.norm.fit(data2)
    xmin, xmax = pl.xlim()
    x = np.linspace(xmin, xmax, 100)
    p = stats.norm.pdf(x, mu2, std2)
    pl.plot(x, p, 'k', linewidth=2, color='r')

    pl.title(_title)
    pl.savefig(data_DIR + '/'+ _title + '.png',bbox_inches='tight')

    pl.close()
    return



def twoFeatureScatterPlotMultiClusters(raw_features, featureNames,featureId1,featureId2,  label, label1, label2, class_id_list, fig_title):
    class_f1=[]
    class_f2=[]
    colors=['blue','red','green','yellow','black','cyan','magenta','gray']
    fig = pl.figure()
    pl.show()
    i=0
    handles=[]
    for  LABEL in class_id_list:
	    pos = (np.array(np.where(label== LABEL)) ).flatten() # inhibitory
	    class_f1 = raw_features[pos, featureId1-1]
	    class_f2 = raw_features[pos, featureId2-1]
	    pl.scatter( class_f1, class_f2 ,s=50, c=colors[i] )
	    i=i+1



    pl.title(fig_title)
    pl.xlabel(featureNames[featureId1])
    pl.ylabel(featureNames[featureId2])
    fig.legend(class_id_list, 'upper left')
    pl.savefig(data_DIR + '/'+ fig_title + '.png',bbox_inches='tight')

    #pl.close()
    return




def twoFeatureScatterPlot(raw_features, label, featureNames, featureId1,featureId2, fig_title):
    pos = (np.array(np.where(label==1)) -1 ).flatten() # inhibitory
    inh1 = raw_features[pos, featureId1-1]
    inh2 = raw_features[pos, featureId2-1]

    pos = (np.array(np.where(label==2)) -1 ).flatten()
    exc1 = raw_features[pos, featureId1-1]
    exc2 = raw_features[pos, featureId2-1]

    pl.figure()
    pl.show()
    pl.scatter(inh1,inh2,s=50,c=u'b')
    pl.scatter(exc1,exc2,s=50,c=u'r')
    pl.title(fig_title)
    pl.xlabel(featureNames[featureId1])
    pl.ylabel(featureNames[featureId2])
    pl.savefig(data_DIR + '/'+ fig_title + '.png',bbox_inches='tight')
    #pl.close()
    return


#==============================================================================

def lda(X,y):
    lda = LDA(n_components=3)
    X_r2 = lda.fit(X,y).transform(X)

    pl.figure()
    for c, i, target_name in zip("gbr", [0,1, 2], ['others','inhibitory','excitatory']):
        pl.scatter(X_r2[y == i, 0], X_r2[y == i, 1], c=c, label=target_name)
        pl.legend()
        pl.title('LDA')


def pca(X):
    pca = PCA(n_components=2)
    X_r = pca.fit(X).transform(X)

    pv_tag = df_complete.IsPVALB
    idxs = np.nonzero(pv_tag== 'Y')[0]
    no_inds = np.nonzero(pv_tag== 'N')[0]


    pl.figure()
    for c, i, target_name in zip("gbr", [0,1, 2], ['pvalb','non-pvalb','spiny']):
        pl.scatter(X_r[y == i, 0], X_r[y == i, 1], c=c, label=target_name)
        pl.legend()
        pl.title('PCA')



def selectFeatures_MRMR(featureArray, threshold=0, number_of_features=5, selection_method='MID'):
    #write out feature array into a csv file, then execute MRMR
    csvfile = data_DIR+"/zscore_glFeatures.csv"

    np.savetxt(csvfile, featureArray, delimiter=",")
    # call MRMR
    # cmd = MRMR +  " -i "+ csvfile + " -t "+ threshold + " -n " + number_of_features
    # print cmd
    #os.system(cmd)
    return

def readDBFeatures(feature_file):
    # TODO: detect nan values
    glf_featureList = []  # each row is a feature vector
    gmi_featureList = []
    with open (feature_file,'r') as  f:
        for fn_line in f: # ignore the SWCFILE=* line
            line_globalFeature = (f.next()).strip()
            glf = map(float,line_globalFeature.split('\t'))
            glf_featureList.append(glf)

            line_GMI = (f.next()).strip()
            gmi = map(float,line_GMI.split('\t'))
            gmi_featureList.append(gmi)

    return  np.array(glf_featureList), np.array(gmi_featureList)



def concatCSVs(csv1, csv2, outcsv):
    df1 = pd.read_csv(csv1)
    df2 = pd.read_csv(csv2)

    #out_df = pd.merge(df1, df2)
    out_df = pd.concat([df1,df2], axis=1)
    out_df.to_csv(outcsv, index=False)
    return

def saveScoreMatrix(featureArray,scoreMatrix_file, REMOVE_OUTLIER=1):
    scoreMatrix =[]
    normalized = zscore(featureArray)
    normalized[normalized < -3]  =-3
    normalized[normalized > 3] = 3

    for i in range(len(normalized)):
        queryFeature = normalized[i] # each row is a feature vecto
        #scores = np.exp(-np.sum(abs(normalized-queryFeature)**2,1)/100)
        scores = np.sum(abs(normalized-queryFeature)**2,1)
        scoreMatrix.append(scores)

    df = pd.DataFrame(scoreMatrix)
    df.to_csv(scoreMatrix_file)

def  generateCompleFeatureCSV(new_csv_file):
    # attacheh specimen id, nrrd id,
    glfFeatures, gmiFeatures = readDBFeatures(feature_file)
    feature_names=['num_nodes',	'soma_surface',	'num_stems','num_bifurcations',	'num_branches',	'num_of_tips',	'overall_width',	'overall_height',	'overall_depth', 'average_diameter',	'total_length',	'total_surface',	'total_volume',	'max_euclidean_distance',	'max_path_distance',	'max_branch_order',	'average_contraction', 'average fragmentation',	'parent_daughter_ratio',	'bifurcation_angle_local',	'bifurcation_angle_remote']


    #attach  specimen_name
    list_df = pd.read_csv(list_csv_file)
    sp_ids = list_df.specimen_id
    sp_names = np.empty( (len(sp_ids),),dtype='object')
    for i in range(len(sp_ids)):
        sp_names[i]= lu.get_specimen_name_from_lims(str(sp_ids[i]))

    df = pd.DataFrame(glfFeatures, columns = feature_names)
    df['specimen_name'] = pd.Series(sp_names, index=df.index)
    #df = df['specimen_name',feature_names]
    df.to_csv(feature_csv_file, index=False)
    concatCSVs(list_csv_file,feature_csv_file, new_csv_file)
    return


def generateLinkerFileFromCSV(result_dir, csvfile, column_name):
	df = pd.read_csv(csvfile)
	types = df[column_name]
	for atype in np.unique(types):
		idxs = np.nonzero(types==atype)[0]
		swc_files = df['orca_path']
		with open(result_dir+'/'+atype+'.ano','w') as outf:
		    for afile in swc_files[idxs]:
			   line='SWCFILE='+afile+'\n'
			   outf.write(line)
		    outf.close()



##################################################################################################


feature_csv_with_id_file = data_DIR+"/glFeatures_withid_pvalb.csv"
#generateCompleFeatureCSV(new_csv_file)


df_complete = pd.read_csv(feature_csv_with_id_file)
mycolumns = ['specimen_id','specimen_name','nrid','orca_path','IsPVALB']
mycolumns.extend(feature_names)
df_complete = df_complete.reindex(columns=mycolumns)
print df_complete.columns

glFeatures = (df_complete.values[:,5:26]).astype(float)
zscores = zscore(glFeatures)



pv_tag = df_complete.IsPVALB
idxs = np.nonzero(pv_tag== 'Y')[0]
pv_features = glFeatures[idxs,:]
no_inds = np.nonzero(pv_tag== 'N')[0]
nonpv_features = glFeatures[no_inds,:]

for cur_feature in feature_names:
	fea_idx = np.nonzero(feature_names == cur_feature)[0][0]
    	print "plot ", fea_idx
    	compareHist(pv_features[:,fea_idx], nonpv_features[:,fea_idx], cur_feature)

# merge all info
df_type = pd.read_csv(data_DIR+'/custom_report-IVSCC_classification-April_2015.csv')
merged = pd.merge(df_complete,df_type,how='inner',on=['specimen_name'])
merged.to_csv(data_DIR+'/merged.csv')

glFeatures = merged.values[:,5:26].astype(float)
scoreMatrix_csv_file = data_DIR+"/scoreMatrix.csv"
saveScoreMatrix(glFeatures,scoreMatrix_csv_file,1)

pv_tag = merged.IsPVALB
idxs = np.nonzero(pv_tag== 'Y')[0]
pv_features = glFeatures[idxs,:]
no_inds = np.nonzero(pv_tag== 'N')[0]
nonpv_features = glFeatures[no_inds,:]


glFeatures_n= normalizeFeatures(glFeatures)
idx1=np.nonzero(glFeatures_n>3)
glFeatures_n[idx1]=3
idx2=np.nonzero(glFeatures_n<-3)
glFeatures_n[idx2]=-3



#### pca
pca = PCA(n_components=2)
scores =pca.fit(glFeatures_n).transform(glFeatures_n)

pl.scatter(scores[idxs,0],scores[idxs,1],c='y',label='Pvlab',s=50)
pl.scatter(scores[no_inds,0],scores[no_inds,1],c='r',label='other aspiny',s=50)

spiny_tag=merged.dendrite_type_tag
spiny_idxs = np.nonzero(spiny_tag== 'spiny')[0]
aspiny_idxs =np.nonzero(spiny_tag== 'aspiny')[0]
pl.scatter(scores[spiny_idxs,0],scores[spiny_idxs,1],c='b',label='spiny',s=50)
#pl.scatter(scores[aspiny_idxs,0],scores[aspiny_idxs,1],c='c',label='aspiny')
pl.xlabel('PC1')
pl.ylabel('PC2')
pl.legend()




##################  scatter plot for top two features
pl.figure()
pl.scatter(glFeatures_n[idxs,13],glFeatures_n[idxs,15],c='y',label='Pvlab',s=50 )
pl.scatter(glFeatures_n[no_inds,13],glFeatures_n[no_inds,15],c='r',label='other aspiny',s=50)
pl.scatter(glFeatures_n[spiny_idxs,13],glFeatures_n[spiny_idxs,15],c='b',label='spiny',s=50 )
#pl.scatter(glFeatures_n[aspiny_idxs,13],glFeatures_n[aspiny_idxs,15],c='c',label='aspiny',s=50 )
pl.xlabel('Max Euclidean Path Dist')
pl.ylabel('Max Branch Order')
pl.legend()



#normalized = normalizeFeatures(glFeatures)
#selectFeatures_MRMR(zscores)  # save zscores into a csv file
#plotFeatureVector(glFeatures, "normalized global features")
#plotScoreMatrix(glFeatures, "Similarity Matrix based on global features")
#plotScoreMatrix(gmiFeatures, "Similarity Matrix based on GMI features")


########################  Layer?
pca = PCA(n_components=5)
pv_features=glFeatures_n[idxs,:]
pvalb_all=merged.loc[idxs,:]
pv_scores =pca.fit(pv_features).transform(pv_features)

layer_group = pvalb_all.layer
l23 = np.nonzero(layer_group== 'L2/3')[0]
l4 =np.nonzero(layer_group== 'L4')[0]
l5 =np.nonzero(layer_group== 'L5')[0]
l6 =np.nonzero(layer_group== 'L6')[0]
pl.figure()
pl.scatter(pv_scores[l23,0],pv_scores[l23,1],c='b',label='L2/3',s=50)
pl.scatter(pv_scores[l4,0],pv_scores[l4,1],c='y',label='L4',s=50)
pl.scatter(pv_scores[l5,0],pv_scores[l5,1],c='c',label='L5',s=50)
pl.scatter(pv_scores[l6,0],pv_scores[l6,1],c='r',label='L6',s=50)


pl.xlabel('PC1')
pl.ylabel('PC2')
pl.legend()

############  region
pca = PCA(n_components=5)
pv_features=glFeatures_n[idxs,:]
pvalb_all=merged.loc[idxs,:]
pv_scores =pca.fit(pv_features).transform(pv_features)

region_group = pvalb_all.region
region1 = np.nonzero(region_group== 'AUDd')[0]
region2 =np.nonzero(region_group== 'VISl')[0]
region3 =np.nonzero(region_group== 'VISp')[0]
region4 =np.nonzero(region_group== 'VISal')[0]
pl.figure()
pl.scatter(pv_scores[region1,0],pv_scores[region1,1],c='b',label='AUDd',s=50)
pl.scatter(pv_scores[region2,0],pv_scores[region2,1],c='y',label='VIS1',s=50)
pl.scatter(pv_scores[region3,0],pv_scores[region3,1],c='c',label='VISp',s=50)
pl.scatter(pv_scores[region4,0],pv_scores[region4,1],c='r',label='VISa1',s=50)


pl.xlabel('PC1')
pl.ylabel('PC2')
pl.legend()

############  cell shape
pca = PCA(n_components=5)
pv_features=glFeatures_n[idxs,:]
pvalb_all=merged.loc[idxs,:]
pv_scores =pca.fit(pv_features).transform(pv_features)

region_group = pvalb_all.cell_shape
region1 = np.nonzero(region_group== 'asymm')[0]
region2 =np.nonzero(region_group== 'symm')[0]
region3 =np.nonzero(region_group== 'bitufted')[0]

pl.figure()
pl.scatter(pv_scores[region1,0],pv_scores[region1,1],c='b',label='asymm',s=50)
pl.scatter(pv_scores[region2,0],pv_scores[region2,1],c='y',label='symm',s=50)
pl.scatter(pv_scores[region3,0],pv_scores[region3,1],c='c',label='bifuted',s=50)



pl.xlabel('PC1')
pl.ylabel('PC2')
pl.legend()

