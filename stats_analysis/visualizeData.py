import numpy as np
import matplotlib.pylab as pl
import matplotlib.pyplot as pplt

from scipy import stats

import pandas as pd
from sklearn.lda import LDA
from sklearn.decomposition import PCA

from os import sys, path
p = path.dirname(path.dirname(path.abspath(__file__)))
sys.path.append(p)
sys.path.append(p+'/utilities')

import lims_orca_utils as lu

# program path on this machine
#===================================================================
#WORK_PATH="/Users/xiaoxiaoliu/work"
WORK_PATH="/home/xiaoxiaol/work"

MRMR= WORK_PATH+"/src/mrmr_c_src/mrmr"
V3D="qv3d"

########################################## data dir
#data_DIR= WORK_PATH+"/data/lims2/neuron_recon_2"
data_DIR= "/home/xiaoxiaol/work/data/lims2/nr_june_25_filter_aligned"
LIST_CSV_FILE =  data_DIR+'/list.csv'
#########################################################


data_linker_file =  data_DIR+'/original/mylinker.ano'
preprocessed_data_linker_file = data_DIR+'/preprocessed/mylinker.ano'
FEATURE_FILE = data_DIR + '/preprocessed/prep_features.nfb'
gl_feature_names= np.array(['num_nodes', 'soma_surface', 'num_stems','num_bifurcations', 'num_branches', 'num_of_tips',  'overall_width', 'overall_height',  'overall_depth', 'average_diameter',    'total_length', 'total_surface', 'total_volume', 'max_euclidean_distance',       'max_path_distance', 'max_branch_order',  'average_contraction', 'average fragmentation', 'parent_daughter_ratio', 'bifurcation_angle_local', 'bifurcation_angle_remote'])
gmi_feature_names = np.array(['moment1', 'moment2', 'moment3','moment4', 'moment5', 'moment6',  'moment7', 'moment8',  'moment9', 'moment10',    'moment11', 'moment12', 'moment13', 'avgR'])

selected_features=['max_euclidean_distance','num_stems','num_bifurcations','average_contraction','parent_daughter_ratio']


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

def readDBFeatures(FEATURE_FILE):
    # TODO: detect nan values
    glf_featureList = []  # each row is a feature vector
    gmi_featureList = []
    with open (FEATURE_FILE,'r') as  f:
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

    # remove outliers!!!
    normalized[normalized < -3]  =-3
    normalized[normalized > 3] = 3

    for i in range(len(normalized)):
        queryFeature = normalized[i] # each row is a feature vecto
        #scores = np.exp(-np.sum(abs(normalized-queryFeature)**2,1)/100)
        scores = np.sum(np.abs(normalized-queryFeature)**2,1)
        scoreMatrix.append(scores)

    df = pd.DataFrame(scoreMatrix)
    df.to_csv(scoreMatrix_file)

def  generateALLFeatureCSV(new_csv_file):
    feature_csv_file= data_DIR+"/allFeatures.csv"  # mid result

    glFeatures, gmiFeatures = readDBFeatures(FEATURE_FILE)

    # attacheh specimen id, nrrd id to the tabe
    list_df = pd.read_csv(LIST_CSV_FILE)
    sp_ids = list_df.specimen_id
    sp_names = np.empty( (len(sp_ids),),dtype='object')
    for i in range(len(sp_ids)):
        sp_names[i]= lu.get_specimen_name_from_lims(str(sp_ids[i]))

    allFeatures = np.append(glFeatures,gmiFeatures,1)
    allColums = np.append(gl_feature_names,gmi_feature_names,0) 
    df = pd.DataFrame(allFeatures,  columns = allColums )
    df['specimen_name'] = pd.Series(sp_names, index=df.index)

    df.to_csv(feature_csv_file, index=False)

    concatCSVs(LIST_CSV_FILE,feature_csv_file, new_csv_file)
    print 'output all feature csv file to :',new_csv_file
    return

def  generateGlFeatureCSV(new_csv_file):
    gl_feature_csv_file	= data_DIR+"/glFeatures.csv"  # mid result
    # attacheh specimen id, nrrd id,
    glfFeatures, gmiFeatures = readDBFeatures(FEATURE_FILE)


    #attach  specimen_name
    list_df = pd.read_csv(LIST_CSV_FILE)
    sp_ids = list_df.specimen_id
    sp_names = np.empty( (len(sp_ids),),dtype='object')
    for i in range(len(sp_ids)):
        sp_names[i]= lu.get_specimen_name_from_lims(str(sp_ids[i]))

    df = pd.DataFrame(glfFeatures, columns = gl_feature_names)
    df['specimen_name'] = pd.Series(sp_names, index=df.index)
    #df = df['specimen_name',feature_names]
    df.to_csv(gl_feature_csv_file, index=False)

    concatCSVs(LIST_CSV_FILE,gl_feature_csv_file, new_csv_file)
    return

def  generateGmiFeatureCSV(new_csv_file):
    gmi_feature_csv_file= data_DIR+"/gmiFeatures.csv" # mid result

    # attacheh specimen id, nrrd id,
    glFeatures, gmiFeatures = readDBFeatures(FEATURE_FILE)

    #attach  specimen_name
    list_df = pd.read_csv(LIST_CSV_FILE)
    sp_ids = list_df.specimen_id
    sp_names = np.empty( (len(sp_ids),),dtype='object')
    for i in range(len(sp_ids)):
        sp_names[i]= lu.get_specimen_name_from_lims(str(sp_ids[i]))

    df = pd.DataFrame(gmiFeatures, columns = gmi_feature_names)
    df['specimen_name'] = pd.Series(sp_names, index=df.index)
    #df = df['specimen_name',feature_names]
    df.to_csv(gmi_feature_csv_file, index=False)

    concatCSVs(LIST_CSV_FILE,gmi_feature_csv_file, new_csv_file)
    return


def generateLinkerFileFromCSV(result_dir, csvfile, column_name):
	df = pd.read_csv(csvfile)
	types = df[column_name]
	for atype in np.unique(types):
		idxs = np.nonzero(types==atype)[0]
		swc_files = df['orca_path']
		with open(result_dir+'/'+atype+'.ano','w') as outf:
        	   for afile in swc_files[idxs]:
                       filename = afile.split('/')[-1]
                       line='SWCFILE='+filename+'\n'
                       outf.write(line)
                   outf.close()


def generateFeatureMergedCSV(out_featureId_filename, ):
	### only for gl features
	feature_csv_with_id_file = data_DIR+"/glFeatures_withid_pvalb.csv"
	generateGlFeatureCSV(out_featureId_filename)

	df_complete = pd.read_csv(feature_csv_with_id_file)
	mycolumns = np.array(['specimen_id','specimen_name','nrid','orca_path','IsPVALB'])
	mycolumns = np.append(mycolumns,gl_feature_names,0)
	df_complete = df_complete.reindex(columns=mycolumns)
	print df_complete.columns


	# merge all info
	df_type = pd.read_csv(data_DIR+'/custom_report-IVSCC_classification-April_2015.csv')
	merged = pd.merge(df_complete,df_type,how='inner',on=['specimen_name'])
	merged.to_csv(data_DIR+'/merged_glFeatures.csv')



	glFeatures = merged.values[:,5:26].astype(float)
	scoreMatrix_csv_file = data_DIR+"/gl_scoreMatrix.csv"
	saveScoreMatrix(glFeatures,scoreMatrix_csv_file,1)

							###########  gmi features
	gmi_feature_csv_with_id_file = data_DIR+"/gmiFeatures_withid.csv"
	generateGmiFeatureCSV(gmi_feature_csv_with_id_file)

	df_complete = pd.read_csv(gmi_feature_csv_with_id_file)
	mycolumns = ['specimen_id','specimen_name','nrid','orca_path','IsPVALB']
	mycolumns.extend(gmi_feature_names)
	df_complete = df_complete.reindex(columns=mycolumns)
	print df_complete.columns



	# merge all info
	df_type = pd.read_csv(data_DIR+'/custom_report-IVSCC_classification-April_2015.csv')
	merged = pd.merge(df_complete,df_type,how='inner',on=['specimen_name'])
	merged.to_csv(data_DIR+'/merged_gmiFeatures.csv')



	gmiFeatures = merged.values[:,5:19].astype(float)
	scoreMatrix_csv_file = data_DIR+"/gmi_scoreMatrix.csv"
	saveScoreMatrix(gmiFeatures,scoreMatrix_csv_file,1)


##################################################################################################
all_feature_csv_with_id_file = data_DIR+"/allFeatures_withid.csv"
generateALLFeatureCSV( data_DIR+"/allFeatures_withid.csv")

df_complete = pd.read_csv(all_feature_csv_with_id_file)
mycolumns = np.array(['specimen_id','specimen_name','orca_path'])
mycolumns = np.append(mycolumns,gl_feature_names,0)
mycolumns = np.append(mycolumns,gmi_feature_names,0)
df_complete = df_complete.reindex(columns=mycolumns)
print df_complete.columns



# merge all info
df_type = pd.read_csv(data_DIR+'/../custom_report-IVSCC_classification-April_2015.csv')
merged = pd.merge(df_complete,df_type,how='inner',on=['specimen_name'])
merged.to_csv(data_DIR+'/merged_allFeatures.csv')


generateLinkerFileFromCSV(data_DIR+'/original',data_DIR +'/merged_allFeatures.csv','cre_line')










